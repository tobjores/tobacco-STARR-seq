#!/bin/sh

# required tools: fastx-toolkit (http://hannonlab.cshl.edu/fastx_toolkit/), bioawk (https://github.com/lh3/bioawk), cutadapt (https://cutadapt.readthedocs.io/en/stable/), and TrimGalore (https://github.com/FelixKrueger/TrimGalore)

# this script was used for linking the fragments of the pZS11_4enh plasmid to barcodes
# the plasmid library was sequenced with 155bp paired end reads
# read 1 was used to sequence the fragments
# read 2 was used to sequence the barcode (15bp) and also contained the 3' UTR, the minimal promoter, and and the 3' end (36bp) of the recombined fragments

SAMPLE="base name of sequencing files"


### reverse-complement read 2 ###
zcat ${SAMPLE}_R2_001.fastq.gz | fastx_reverse_complement -z > ${SAMPLE}_R2_revcom.fastq.gz

### split out barcode ###
bioawk -c fastx -t '{print $name, substr($seq, length($seq) - 14, 15)}' ${SAMPLE}_R2_revcom.fastq.gz | sort > barcodes_${SAMPLE%_*}.tsv

### trim reads to remove nextera adapter trailing sequences (also removes the promoter, UTR, and barcode from read 2) ###
trim_galore --quality 0 --length 0 --dont_gzip --nextera --paired ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_revcom.fastq.gz

### align reads to reference sequence ###
if [ ! -d bowtie2-index ]; then
	mkdir bowtie2-index
fi

bowtie2-build pZS11_4enh.fa bowtie2-index/pZS11_4enh

bowtie2 --xeq -x bowtie2-index/pZS11_4enh --no-discordant --no-mixed -X 1000 --ff -S ${SAMPLE%_*}.sam -1 ${SAMPLE}_R1_001_val_1.fq -2 ${SAMPLE}_R2_revcom_val_2.fq

### generate .bed files for fragments ###
bioawk -c sam -t '{if(and($flag, 2) && and($flag, 64) && $cigar !~ /I|D/) {if($tlen > 0) print $rname, $pos - 1, $pos + $tlen - 1, $qname, $mapq, "+"; else print $rname, $pnext - 1, $pnext - $tlen - 1, $qname, $mapq, "-";}}' ${SAMPLE%_*}.sam > ${SAMPLE%_*}.bed

### convert .bed file ###
awk -v OFS='\t' '{print $4, $2, $3, $6}' ${SAMPLE%_*}.bed | sort > frags_${SAMPLE%_*}.tsv

### merge barcode and fragment files (has to have a tab character between quotes) ###
join -t '	' --nocheck-order barcodes_${SAMPLE%_*}.tsv frags_${SAMPLE%_*}.tsv > joined_${SAMPLE%_*}.tsv

### combine and count unique fragments (min 5 reads) ###
awk '{print $2, $3, $4, $5}' joined_${SAMPLE%_*}.tsv | sort | uniq -c | awk -v OFS='\t' 'BEGIN{print "barcode", "start", "stop", "strand", "assembly.count"} {if($1 >= 5) print $2, $3, $4, $5, $1}' > subassembly_${SAMPLE%_*}.tsv