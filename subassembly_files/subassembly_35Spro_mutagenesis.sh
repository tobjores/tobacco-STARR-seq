#!/bin/sh

# required tools: bioawk (https://github.com/lh3/bioawk), cutadapt (https://cutadapt.readthedocs.io/en/stable/), TrimGalore (https://github.com/FelixKrueger/TrimGalore), and PANDAseq (https://github.com/neufeld/pandaseq)

# this script was used for linking the 35S minimal promoter single base variants to barcodes
# the plasmid library was sequenced with 74/86bp paired end reads
# read 1 was used to sequence the promoter variant
# read 2 was used to sequence the barcode (15bp) and also contained the 5' UTR, and 3' end (40bp) of the promoter variant

SAMPLE="base name of sequencing files"

### assemble reads ###
if [ ! -d panda_logs ]; then
	mkdir panda_logs
fi

pandaseq -F -f ${SAMPLE}_R1_001.fastq.gz -r ${SAMPLE}_R2_001.fastq.gz -g panda_logs/log_${SAMPLE%_*}.txt > assembled_${SAMPLE%_*}.fq

### trim reads  (removes the trainling promoter, UTR, and barcode)###
trim_galore --length 45 --quality 0 --adapter ACACGCTGGAATTCTAGTATACTAAACCATG assembled_${SAMPLE%_*}.fq

### extract read name and sequence ###
bioawk -c fastx -t '{print substr($name, 1, length($name) - 9), $seq}' assembled_${SAMPLE%_*}_trimmed.fq | sort > variants_${SAMPLE%_*}.tsv

### trim read 2 to the barcode only ###
trim_galore --length 0 --quality 0 --adapter CATGGTTTAGTATACTAGAATTCCAGCGTGTCCTCTC ${SAMPLE}_R2_001.fastq.gz

### extract barcodes of correct length ###
bioawk -c fastx -t '{if(length($seq) == 15 && $seq !~ /N/) print $name, revcomp($seq)}' ${SAMPLE}_R2_001_trimmed.fq | sort > barcodes_${SAMPLE%_*}.tsv

### join barcode and variant data (has to have a tab character between quotes) ###
join -t '	' --nocheck-order barcodes_${SAMPLE%_*}.tsv variants_${SAMPLE%_*}.tsv > joined_${SAMPLE%_*}.tsv

### combine and count unique reads (min 5 reads) ###
awk '{print $2, $3}' joined_${SAMPLE%_*}.tsv | sort | uniq -c | awk -v OFS='\t' 'BEGIN{print "barcode", "sequence", "assembly.count"} {print $2, $3, $1}' > count_${SAMPLE%_*}.tsv
awk 'BEGIN{getline; print $0} {if($3 >= 5) print $0}' count_${SAMPLE%_*}.tsv > subassembly_${SAMPLE%_*}.tsv