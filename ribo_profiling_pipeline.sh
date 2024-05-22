#!/bin/bash

########################################################################
# Bioinformatics pipeline for Ribo-seq data:
# - (1) download raw reads from SRR inputs
# - (2) filter and trim
# - (4) align and output quantification of reads
# Notes: 
# - requires standard input file with one SRR per line
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
    
    echo "Ribo-seq Analysis of ${srr}"



    ########################
    # Step 1 
    ########################
    echo "1. Aquire raw data"

    mkdir -p ${srr}_ribo
    cd ${srr}_ribo || exit
    
    if [ ! -f "${srr}_1.fastq" ]; then
        prefetch --progress ${srr}
        fasterq-dump --split-files --threads 8 --progress ${srr}
    else
        echo "${srr}_1.fastq already exists. Skipping download."
    fi

    if [ ! -f "${srr}_1_fastqc.html" ]; then
        echo "Running quality check on pairs ${srr}_1.fastq"
        fastqc -q -t 8 ${srr}_1.fastq
    else
        echo "Skipping quality check on ${srr}_1.fastq"
    fi



    #######################
    # Step 2
    #######################
    echo "2. Trim adapters"

    # Trim adapter and length (Trim 5nt from the 5' and 4nt from the 3' end)
    cutadapt --report minimal \
             --cores 8 \
             -a file:/Users/scampione/data/unique_sequencing_adapters.fasta \
             -u 5 -u -4 \
             --minimum-length 25 \
             -o trimmed_adapter_${srr}.fastq \
             ${srr}_1.fastq

    # cutadapt --report minimal \
    #          --cores 8 \
    #          -a file:/Users/scampione/data/ribo_seq_adapters.fasta \
    #          -u 5 -u -4 \
    #          --minimum-length 25 \
    #          -o trimmed_adapter_${srr}.fastq \
    #          ${srr}_1.fastq



    # Trim polyG 
    cutadapt --report minimal \
             -a GGGGGG \
             --minimum-length 25 \
             -o trimmed_${srr}.fastq trimmed_adapter_${srr}.fastq


    # Remove reads with unknown bases NNNNNNNNN
    fastp -i trimmed_${srr}.fastq -o filtered_${srr}.fastq --n_base_limit 0

    echo "Running quality check report on trimmed filtered reads"    
    fastqc -t 8 filtered_${srr}.fastq



    # Ribosomal RNA removal (output is non-rRNA reads)
    
    echo "Aligning against rRNA_Saccharomyces_cerevisiae_index"
    bowtie2 -x /Users/scampione/data/rRNA_Saccharomyces_cerevisiae_index/rRNA.Saccharomyces_cerevisiae \
            -p 8 \
            -U filtered_${srr}.fastq \
            -S rrna_mapped_${srr}.sam \
            --un non_rrna_reads_${srr}.fastq


    # Don't use reads under 25 nt
    seqtk seq -L 25 non_rrna_reads_${srr}.fastq > filtered_non_rrna_reads_${srr}.fastq

    echo "Running quality check report on non_rrna_reads_${srr}.fastq"
    fastqc -q -t 8 filtered_non_rrna_reads_${srr}.fastq



    ########################
    # Step 3
    ########################
    echo "4. Alignment and quantification of non rRNA reads"

    echo "Aligning against s_cerevisiae_R64 using bowtie 2"

    # Build updated reference genome RF64-4-1
    # bowtie2-build /Users/scampione/data/s_cer_r64_4_1_bowtie2_index/GCF_000146045.2_R64_genomic.fna \
    #               /Users/scampione/data/s_cer_r64_4_1_bowtie2_index/s_cer_r64_4_1

    bowtie2 -x /Users/scampione/data/s_cer_r64_4_1_bowtie2_index/s_cer_r64_4_1 \
            -U filtered_non_rrna_reads_${srr}.fastq \
            -S aligned_reads_${srr}.sam


    # Filter Unmapped Reads and convert SAM to BAM
    samtools view -bS aligned_reads_${srr}.sam > aligned_reads_${srr}.bam
    samtools sort aligned_reads_${srr}.bam -o sorted_aligned_reads_${srr}.bam
    samtools index sorted_aligned_reads_${srr}.bam
    samtools view -b -F 4 sorted_aligned_reads_${srr}.bam > filtered_aligned_reads_${srr}.bam


    echo "Quantification:"

    echo "Using Blevins_Tavella annotation file stranded=yes"
    htseq-count -f bam -m intersection-strict \
                -s yes \
                -q \
                -i gene_name \
                filtered_aligned_reads_${srr}.bam \
                /Users/scampione/data/Blevins_Tavella_etal_Scer_transcriptome.gtf \
                > yes_stranded_Blevins_Tavella_gtf_counts_${srr}_full.txt


                

    echo "End of analysis of ${srr}"

    cd ..

done < "$1"
