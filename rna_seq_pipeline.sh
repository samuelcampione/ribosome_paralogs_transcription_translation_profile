#!/bin/bash

########################################################################
# Bioinformatics pipeline for RNA-seq data:
# - (1) download raw reads from SRR inputs
# - (2) trim adapters
# - (3) align and output quantification of non-rRNA reads
# Notes: 
# - requires stdin file with one SRR per line
########################################################################



# check if the input argument is provided
if [ $# -eq 0 ]; then
    echo "No file path provided. Exiting."
    exit 1
fi

# get the directory name from input file path
input_dir=$(dirname "$1") 
cd $input_dir

# loop through each line of input file with 1 srr per line
while IFS= read -r srr; do 
    
    echo "RNA-seq Analysis of ${srr}"



    ########################
    # Step 1 
    ########################
    echo "1. Aquire raw data"

    mkdir -p ${srr}_rna
    cd ${srr}_rna || exit
    
    if [ ! -f "${srr}_1.fastq" ]; then
        prefetch --progress ${srr}
        fasterq-dump --split-files --threads 8 --progress ${srr}
    else
        echo "${srr}_1.fastq and ${srr}_2.fastq already exists. Skipping download."
    fi

    if [ ! -f "${srr}_1_fastqc.html" ]; then
        echo "Running quality check on pairs ${srr}"
        fastqc -q -t 8 ${srr}_1.fastq 
        fastqc -q -t 8 ${srr}_2.fastq
    else
        echo "Skipping quality check on ${srr} pairs."
    fi
    
    

    # ########################
    # # Step 2
    # ########################
    # echo "2. Trim adapters"

    # # trim and filter using default parameters
    # trimmomatic PE \
    #     -phred33 \
    #     -threads 8 \
    #     ${srr}_1.fastq ${srr}_2.fastq \
    #     ${srr}_1_paired.fastq ${srr}_1_unpaired.fastq \
    #     ${srr}_2_paired.fastq ${srr}_2_unpaired.fastq \
    #     ILLUMINACLIP:/Users/scampione/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10 \
    #     LEADING:3 TRAILING:3 \
    #     SLIDINGWINDOW:4:15 MINLEN:36
        

    # # quality check on trimmed reads
    # fastqc ${srr}_1_paired.fastq ${srr}_2_paired.fastq


    # ########################
    # # Step 3
    # ########################
    # echo "3. Alignment and quantification"
    
    # # align against S Cer R64 (ver 4.1) genome assembly
    # echo "Aligning against s_cerevisiae_R64 4.1 using bowtie 2"
    # bowtie2 -x /Users/scampione/data/s_cer_r64_4_1_bowtie2_index/s_cer_r64_4_1 \
    #     -1 ${srr}_1_paired.fastq \
    #     -2 ${srr}_2_paired.fastq \
    #     -S ${srr}_aligned_reads.sam



    # # convert SAM to BAM, sort by name (so pairs are together), and filter
    # samtools view -bS ${srr}_aligned_reads.sam > ${srr}_aligned_reads.bam
    # samtools sort -n ${srr}_aligned_reads.bam -o ${srr}_sorted_aligned_reads.bam
    # samtools view -b -F 4 ${srr}_sorted_aligned_reads.bam > ${srr}_filtered_aligned_reads.bam


    # # quantification
    # Mode: intersection-strict, min qual: 10
    # htseq-count -f bam -m intersection-strict \
    #             -s reverse \
    #             -q \
    #             -i gene_name \
    #             ${srr}_filtered_aligned_reads.bam \
    #             /Users/scampione/data/Blevins_Tavella_etal_Scer_transcriptome.gtf \
    #             > ${srr}_gene_counts.txt

    # Mode: intersection-strict, min qual: 0
    htseq-count -f bam -m intersection-strict --minaqual 0 \
                -s reverse \
                -q \
                -i gene_name \
                ${srr}_filtered_aligned_reads.bam \
                /Users/scampione/data/Blevins_Tavella_etal_Scer_transcriptome.gtf \
                > quality_thresh_0_${srr}_gene_counts.txt

    # Mode: Union, Min qual: 0
    # htseq-count -f bam -m union --minaqual 0 \
    #             -s reverse \
    #             -q \
    #             -i gene_name \
    #             ${srr}_filtered_aligned_reads.bam \
    #             /Users/scampione/data/Blevins_Tavella_etal_Scer_transcriptome.gtf \
    #             > union_thesh_0_${srr}_gene_counts.txt
    

    # filtered_read_count=$(samtools view -c ${srr}_filtered_aligned_reads.bam)
    # echo "Filtered read count for ${srr}: ${filtered_read_count}"


    echo "End of analysis of ${srr}"
    
    cd ..

done < "$1"