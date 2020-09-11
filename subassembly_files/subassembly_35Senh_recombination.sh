#!/bin/sh

# required tools: bioawk (https://github.com/lh3/bioawk), cutadapt (https://cutadapt.readthedocs.io/en/stable/), TrimGalore (https://github.com/FelixKrueger/TrimGalore), and PANDAseq (https://github.com/neufeld/pandaseq)

# this script was used for linking the recombined 35S enhancer fragments to barcodes
# the plasmid library was sequenced with 155bp paired end reads
# read 1 was used to sequence the fragments
# read 2 was used to sequence the barcode (15bp) and also contained the 5' UTR, the minimal promoter, and and the 3' end (55bp) of the recombined fragments

SAMPLE="base name of sequencing files"

### split out barcode ###
bioawk -c fastx -t '{print $name, revcomp(substr($seq, 1, 15))}' ${SAMPLE}_R2_001.fastq.gz | sort > barcodes_${SAMPLE%_*}.tsv

### assemble reads ###
if [ ! -d panda_logs ]; then
	mkdir panda_logs
fi

pandaseq -F -f ${SAMPLE}_R1_001.fastq.gz -r ${SAMPLE}_R2_001.fastq.gz -g panda_logs/log_${SAMPLE%_*}.txt > assembled_${SAMPLE%_*}.fq

### trim reads  (removes the trainling promoter, UTR, and barcode)###
trim_galore --quality 0 --length 0 --dont_gzip --adapter GTGATGAGGAGCTGGCAAGACCCTT assembled_${SAMPLE%_*}.fq

### extract read name and sequence ###
bioawk -c fastx -t '{print substr($name, 1, length($name) - 9), $seq}' assembled_${SAMPLE%_*}_trimmed.fq | sort > reads_${SAMPLE%_*}.tsv

### merge barcode and fragment files (there has to be a tab between the quotes) ###
join -t '	' --nocheck-order barcodes_${SAMPLE%_*}.tsv reads_${SAMPLE%_*}.tsv > joined_${SAMPLE%_*}.tsv

### combine and count unique reads (min 5 reads) ###
awk '{print $2, $3}' joined_${SAMPLE%_*}.tsv | grep -v 'N' |  sort | uniq -c | awk -v OFS='\t' 'BEGIN{print "barcode", "sequence", "assembly.count"} {if($1 >= 5) print $2, $3, $1}' > subassembly_${SAMPLE%_*}.tsv

### extract empty construct reads ###
zcat ${SAMPLE}_R1_001.fastq.gz | grep -B 1 -A 2 '^AGGAGCTGGCAAGACCCTTCCTCTATATAAG' | grep -v '\-\-' > empty_${SAMPLE%_*}.fq

### extract read name and sequence ###
bioawk -c fastx -t '{print $name, $seq}' empty_${SAMPLE%_*}.fq | sort > empty_reads_${SAMPLE%_*}.tsv

### merge barcode and empty construct reads (tab between quotes) ###
join -t '	' --nocheck-order barcodes_${SAMPLE%_*}.tsv empty_reads_${SAMPLE%_*}.tsv > empty_joined_${SAMPLE%_*}.tsv

### combine and count unique reads (min 5 reads) ###
awk '{print $2}' empty_joined_${SAMPLE%_*}.tsv | grep -v 'N' |  sort | uniq -c | awk -v OFS='\t' 'BEGIN{print "barcode", "sequence", "assembly.count"} {if($1 >= 5) print $2, "empty", $1}' > subassembly_empty_${SAMPLE%_*}.tsv