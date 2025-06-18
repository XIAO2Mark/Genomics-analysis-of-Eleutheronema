#!/bin/bash
# Genome Assembly Pipeline for E. rhadinum
# This script outlines the hybrid assembly strategy combining MGI-seq 2000 short reads,
# PacBio Sequel II long reads, and Hi-C technology

# Set variables
WORKING_DIR="/path/to/working/directory"
SHORT_READ_DIR="${WORKING_DIR}/short_reads"
LONG_READ_DIR="${WORKING_DIR}/pacbio_reads"
HIC_READ_DIR="${WORKING_DIR}/hic_reads"
GENOME_SIZE=800000000  # Estimated genome size in bp
THREADS=32
OUTPUT_DIR="${WORKING_DIR}/assembly"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Step 1: Quality control of raw reads
echo "Step 1: Quality control of raw reads"

## Short reads QC
mkdir -p ${OUTPUT_DIR}/qc_reports
fastp --in1 ${SHORT_READ_DIR}/R1.fq.gz --in2 ${SHORT_READ_DIR}/R2.fq.gz \
      --out1 ${OUTPUT_DIR}/clean_R1.fq.gz --out2 ${OUTPUT_DIR}/clean_R2.fq.gz \
      --html ${OUTPUT_DIR}/qc_reports/short_reads_qc.html \
      --json ${OUTPUT_DIR}/qc_reports/short_reads_qc.json \
      --thread ${THREADS}

## Long reads QC
NanoPlot --fastq ${LONG_READ_DIR}/*.fastq.gz \
         --outdir ${OUTPUT_DIR}/qc_reports/long_reads \
         --threads ${THREADS}

# Step 2: De novo assembly using PacBio long reads
echo "Step 2: De novo assembly using PacBio long reads"
mkdir -p ${OUTPUT_DIR}/hifiasm
cd ${OUTPUT_DIR}/hifiasm

# Using Hifiasm for PacBio HiFi reads assembly
hifiasm -o erhadinum.asm -t ${THREADS} ${LONG_READ_DIR}/*.fastq.gz

# Convert the output to FASTA format
awk '/^S/{print ">"$2;print $3}' erhadinum.asm.p_ctg.gfa > erhadinum.asm.p_ctg.fa

# Step 3: Polishing the assembly with short reads
echo "Step 3: Polishing the assembly with short reads"
mkdir -p ${OUTPUT_DIR}/polishing

# Map short reads to the assembly
bwa index ${OUTPUT_DIR}/hifiasm/erhadinum.asm.p_ctg.fa
bwa mem -t ${THREADS} ${OUTPUT_DIR}/hifiasm/erhadinum.asm.p_ctg.fa \
        ${OUTPUT_DIR}/clean_R1.fq.gz ${OUTPUT_DIR}/clean_R2.fq.gz | \
        samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/polishing/short_reads.sorted.bam -
samtools index ${OUTPUT_DIR}/polishing/short_reads.sorted.bam

# Polish with Pilon
pilon --genome ${OUTPUT_DIR}/hifiasm/erhadinum.asm.p_ctg.fa \
      --frags ${OUTPUT_DIR}/polishing/short_reads.sorted.bam \
      --output erhadinum.polished \
      --outdir ${OUTPUT_DIR}/polishing \
      --changes --threads ${THREADS}

# Step 4: Scaffolding with Hi-C data
echo "Step 4: Scaffolding with Hi-C data"
mkdir -p ${OUTPUT_DIR}/scaffolding

# Map Hi-C reads to polished assembly
bwa index ${OUTPUT_DIR}/polishing/erhadinum.polished.fasta
bwa mem -t ${THREADS} ${OUTPUT_DIR}/polishing/erhadinum.polished.fasta \
        ${HIC_READ_DIR}/HiC_R1.fastq.gz | samtools sort -n -@ ${THREADS} \
        -o ${OUTPUT_DIR}/scaffolding/hic_R1.sorted.bam -
bwa mem -t ${THREADS} ${OUTPUT_DIR}/polishing/erhadinum.polished.fasta \
        ${HIC_READ_DIR}/HiC_R2.fastq.gz | samtools sort -n -@ ${THREADS} \
        -o ${OUTPUT_DIR}/scaffolding/hic_R2.sorted.bam -

# Process Hi-C reads with HiC-Pro
HiC-Pro -i ${OUTPUT_DIR}/scaffolding/hic_R1.sorted.bam,${OUTPUT_DIR}/scaffolding/hic_R2.sorted.bam \
        -o ${OUTPUT_DIR}/scaffolding/hicpro_results \
        -g ${OUTPUT_DIR}/polishing/erhadinum.polished.fasta

# Perform 3D-DNA scaffolding
3d-dna ${OUTPUT_DIR}/polishing/erhadinum.polished.fasta \
       ${OUTPUT_DIR}/scaffolding/hicpro_results/hic_contacts.bed \
       -o ${OUTPUT_DIR}/scaffolding/final_assembly

# Step 5: Assembly evaluation
echo "Step 5: Assembly evaluation"
mkdir -p ${OUTPUT_DIR}/evaluation

# Calculate basic assembly statistics
assembly-stats ${OUTPUT_DIR}/scaffolding/final_assembly.fasta > ${OUTPUT_DIR}/evaluation/assembly_stats.txt

# Assess assembly completeness with BUSCO
busco -i ${OUTPUT_DIR}/scaffolding/final_assembly.fasta \
      -o busco_results \
      -l actinopterygii_odb10 \
      -m genome \
      -c ${THREADS} \
      --out_path ${OUTPUT_DIR}/evaluation/

# Map reads back to final assembly to assess mapping rates
## Map short reads
bwa index ${OUTPUT_DIR}/scaffolding/final_assembly.fasta
bwa mem -t ${THREADS} ${OUTPUT_DIR}/scaffolding/final_assembly.fasta \
        ${OUTPUT_DIR}/clean_R1.fq.gz ${OUTPUT_DIR}/clean_R2.fq.gz | \
        samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/evaluation/short_reads_to_final.bam -
samtools index ${OUTPUT_DIR}/evaluation/short_reads_to_final.bam
samtools flagstat ${OUTPUT_DIR}/evaluation/short_reads_to_final.bam > ${OUTPUT_DIR}/evaluation/short_reads_mapping_stats.txt

## Map long reads
minimap2 -ax map-pb -t ${THREADS} ${OUTPUT_DIR}/scaffolding/final_assembly.fasta \
         ${LONG_READ_DIR}/*.fastq.gz | \
         samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/evaluation/long_reads_to_final.bam -
samtools index ${OUTPUT_DIR}/evaluation/long_reads_to_final.bam
samtools flagstat ${OUTPUT_DIR}/evaluation/long_reads_to_final.bam > ${OUTPUT_DIR}/evaluation/long_reads_mapping_stats.txt

echo "Genome assembly pipeline completed!"