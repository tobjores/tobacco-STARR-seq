#!/bin/sh

# required tools: cutadapt (https://cutadapt.readthedocs.io/en/stable/), TrimGalore (https://github.com/FelixKrueger/TrimGalore), and PANDAseq (https://github.com/neufeld/pandaseq)

SAMPLE="base name of sequencing files (e.g. Fig1B_Rep1_cDNA)"

# set correct adapter sequences and barcode length
ADAPTER1="GTGAGCAAGGGCGAGGAGCTGTTCACCGGG"
ADAPTER2="CATGGTTTAGTATACTAGAATTCCAGCGTG"
LENGTH=15

# trim barcode reads
trim_galore --length ${LENGTH} --quality 0 --adapter ${ADAPTER1} --adapter2 ${ADAPTER2} --paired ${SAMPLE}_R1.fq.gz ${SAMPLE}_R2.fq.gz

# assemble reads
if [ ! -d panda_logs ]; then
  mkdir panda_logs
fi

pandaseq -f ${SAMPLE}_R1_val_1.fq.gz -r ${SAMPLE}_R2_val_2.fq.gz -g panda_logs/log_${SAMPLE}.txt > barcodes_${SAMPLE}.fa

# filter for barcodes of correct length; combine and count unique barcodes
awk -v bclen=${LENGTH} '{if (NR % 2 == 0 && length($0) == bclen) print}' barcodes_${SAMPLE}.fa | sort | uniq -c > barcodes_${SAMPLE}.count
