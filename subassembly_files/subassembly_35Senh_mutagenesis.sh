#!/bin/sh

# required tools: bioawk (https://github.com/lh3/bioawk), cutadapt (https://cutadapt.readthedocs.io/en/stable/), TrimGalore (https://github.com/FelixKrueger/TrimGalore), and PANDAseq (https://github.com/neufeld/pandaseq)

# this script was used for linking the 35S enhancer single base variants to barcodes
# the plasmid library was sequenced with 150bp paired end reads
# read 1 (R1) was used to sequence the enhancer variant
# read 2 (R3) was used to sequence the barcode (15bp) and also contained the 5' UTR, the minimal promoter, and and the 3' end (50bp) of the enhancer variant
# index 1 (R2, 8bp) was used to sequence the 5' end of the barcodes

SAMPLE="base name of sequencing files"


### assemble reads ###
if [ ! -d panda_logs ]; then
	mkdir panda_logs
fi

pandaseq -F -f ${SAMPLE}_R1_001.fastq.gz -r ${SAMPLE}_R3_001.fastq.gz -g panda_logs/log_${SAMPLE%_*}.txt > assembled_${SAMPLE%_*}.fq

### trim reads  (removes the trainling promoter, UTR, and barcode)###
trim_galore --length 30 --quality 0 --adapter AGGAGCTGGCAAGACCCTTCCTCTATATAAGGAAGTTCATTTCATTTG assembled_${SAMPLE%_*}.fq

### extract read name and sequence ###
bioawk -c fastx -t '{print substr($name, 1, length($name) - 9), $seq}' assembled_${SAMPLE%_*}_trimmed.fq | sort > variants_${SAMPLE%_*}.tsv

### trim read 2 to the barcode only ###
trim_galore --length 0 --quality 0 --adapter CATGGTTTAGTATACTAGAATTCCAGCGTGTCCTCTC ${SAMPLE}_R3_001.fastq.gz
gzip ${SAMPLE}_R3_001_trimmed.fq

### assemble barcodes ###
pandaseq -F -f ${SAMPLE}_R2_001.fastq.gz -r ${SAMPLE}_R3_001_trimmed.fq.gz -g panda_log_${SAMPLE%_*}_barcodes.txt > barcodes_${SAMPLE%_*}.fq

### extract barcodes of correct length ###
bioawk -c fastx -t '{if(length($seq) == 15) print substr($name, 1, length($name) - 9), $seq}' barcodes_${SAMPLE%_*}.fq | sort > barcodes_${SAMPLE%_*}.tsv

### join barcode and variant data (has to have a tab character between quotes) ###
join -t '	' --nocheck-order barcodes_${SAMPLE%_*}.tsv variants_${SAMPLE%_*}.tsv > joined_${SAMPLE%_*}.tsv

### combine and count unique reads (min 10 reads) ###
awk '{print $2, $3}' joined_${SAMPLE%_*}.tsv | sort | uniq -c | awk -v OFS='\t' 'BEGIN{print "barcode", "sequence", "assembly.count"} {print $2, $3, $1}' > count_${SAMPLE%_*}.tsv
awk 'BEGIN{getline; print $0} {if($3 >= 10) print $0}' count_${SAMPLE%_*}.tsv > subassembly_${SAMPLE%_*}.tsv