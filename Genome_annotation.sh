#!/bin/bash
# Genome Annotation Pipeline for E. rhadinum

WORKING_DIR="/Home/jiexiao/annotation"
GENOME_FASTA="${WORKING_DIR}/assembly/scaffolding/final_assembly.fasta"
RNA_SEQ_DIR="${WORKING_DIR}/RNA_data"
PROTEIN_DB="${WORKING_DIR}/databases/protein_db"
OUTPUT_DIR="${WORKING_DIR}/annotation"
THREADS=32

mkdir -p ${OUTPUT_DIR}

#Repeat identification and masking
mkdir -p ${OUTPUT_DIR}/repeats

# De novo repeat identification
BuildDatabase -name erhadinum_db -engine ncbi ${GENOME_FASTA}
RepeatModeler -database erhadinum_db -engine ncbi -pa ${THREADS} -LTRStruct
cat RepeatModeler-*/consensi.fa.classified ${WORKING_DIR}/databases/RepBase/RepBase_fish.fa > ${OUTPUT_DIR}/repeats/repeat_library.fa
RepeatMasker -pa ${THREADS} -lib ${OUTPUT_DIR}/repeats/repeat_library.fa -xsmall -gff -html -dir ${OUTPUT_DIR}/repeats/ ${GENOME_FASTA}

# Generate statistics
perl ${WORKING_DIR}/tools/RepeatMasker/util/buildSummary.pl -species fish ${OUTPUT_DIR}/repeats/final_assembly.fasta.out > ${OUTPUT_DIR}/repeats/repeat_summary.txt

# RNA-Seq for gene prediction
mkdir -p ${OUTPUT_DIR}/transcriptome
hisat2-build ${GENOME_FASTA} ${OUTPUT_DIR}/transcriptome/genome_idx
hisat2 -p ${THREADS} -x ${OUTPUT_DIR}/transcriptome/genome_idx \
       -1 ${RNA_SEQ_DIR}/RNA_R1.fastq.gz -2 ${RNA_SEQ_DIR}/RNA_R2.fastq.gz \
       | samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/transcriptome/rnaseq.sorted.bam -
samtools index ${OUTPUT_DIR}/transcriptome/rnaseq.sorted.bam

stringtie ${OUTPUT_DIR}/transcriptome/rnaseq.sorted.bam \
          -o ${OUTPUT_DIR}/transcriptome/stringtie_assembly.gtf \
          -p ${THREADS}

# Convert to FASTA for MAKER input
gffread -w ${OUTPUT_DIR}/transcriptome/transcripts.fa \
        -g ${GENOME_FASTA} ${OUTPUT_DIR}/transcriptome/stringtie_assembly.gtf

#Protein-coding gene prediction using MAKER
mkdir -p ${OUTPUT_DIR}/maker

# Create MAKER control files
cd ${OUTPUT_DIR}/maker
maker -CTL
mpiexec -n ${THREADS} maker -base erhadinum_maker

# Process MAKER output
cd ${OUTPUT_DIR}/maker/erhadinum_maker.maker.output
gff3_merge -d erhadinum_maker_master_datastore_index.log -o erhadinum_genes.gff
fasta_merge -d erhadinum_maker_master_datastore_index.log -o ${OUTPUT_DIR}/maker/erhadinum_proteins.fa

# Extract gene and transcript sequences
./gff3_to_fasta ${GENOME_FASTA} erhadinum_genes.gff -t cds -o ${OUTPUT_DIR}/maker/erhadinum_cds.fa
./gff3_to_fasta ${GENOME_FASTA} erhadinum_genes.gff -t mrna -o ${OUTPUT_DIR}/maker/erhadinum_transcripts.fa

#Non-coding RNA prediction
mkdir -p ${OUTPUT_DIR}/ncrna

# tRNA prediction with tRNAscan-SE
tRNAscan-SE -o ${OUTPUT_DIR}/ncrna/trna.out \
            -f ${OUTPUT_DIR}/ncrna/trna.ss \
            -m ${OUTPUT_DIR}/ncrna/trna_stats.txt \
            ${GENOME_FASTA}

# rRNA prediction with RNAmmer
rnammer -S euk -m lsu,ssu,tsu -gff ${OUTPUT_DIR}/ncrna/rrna.gff \
        -h ${OUTPUT_DIR}/ncrna/rrna.hmmreport \
        -f ${OUTPUT_DIR}/ncrna/rrna.fasta \
        ${GENOME_FASTA}
        
cmscan --cpu ${THREADS} --tblout ${OUTPUT_DIR}/ncrna/rfam.tblout \
       --fmt 2 --clanin ${WORKING_DIR}/databases/Rfam/Rfam.clanin \
       ${WORKING_DIR}/databases/Rfam/Rfam.cm ${GENOME_FASTA}

./cmsscan2gff.pl ${OUTPUT_DIR}/ncrna/rfam.tblout > ${OUTPUT_DIR}/ncrna/rfam.gff

# Step 5: Synteny analysis with E. tetradactylum
echo "Step 5: Synteny analysis"
mkdir -p ${OUTPUT_DIR}/synteny

# Format databases for BLAST
makeblastdb -in ${OUTPUT_DIR}/maker/erhadinum_proteins.fa \
            -dbtype prot -out ${OUTPUT_DIR}/synteny/erhadinum_proteins_db

makeblastdb -in ${WORKING_DIR}/E_tetradactylum/proteins.fa \
            -dbtype prot -out ${OUTPUT_DIR}/synteny/etetradactylum_proteins_db

# BLASTP search
blastp -query ${OUTPUT_DIR}/maker/erhadinum_proteins.fa \
       -db ${OUTPUT_DIR}/synteny/etetradactylum_proteins_db \
       -evalue 1e-5 -outfmt 6 -num_threads ${THREADS} \
       -out ${OUTPUT_DIR}/synteny/erhadinum_vs_etetradactylum.blastp

# Run MCScanX for synteny analysis
# Prepare input files for MCScanX
cat ${OUTPUT_DIR}/synteny/erhadinum_vs_etetradactylum.blastp > ${OUTPUT_DIR}/synteny/mcscanx_input.blast

# Create GFF for MCScanX
awk '{print $1"\t"$2"\t"$3"\t"$4}' ${OUTPUT_DIR}/maker/erhadinum_genes.gff > ${OUTPUT_DIR}/synteny/erhadinum.gff
awk '{print $1"\t"$2"\t"$3"\t"$4}' ${WORKING_DIR}/E_tetradactylum/genes.gff > ${OUTPUT_DIR}/synteny/etetradactylum.gff
cat ${OUTPUT_DIR}/synteny/erhadinum.gff ${OUTPUT_DIR}/synteny/etetradactylum.gff > ${OUTPUT_DIR}/synteny/combined.gff

# Run MCScanX
cd ${OUTPUT_DIR}/synteny
${WORKING_DIR}/tools/MCScanX/MCScanX mcscanx_input

# Generate dot plot visualization
${WORKING_DIR}/tools/MCScanX/downstream_analyses/dot_plotter mcscanx_input.collinearity
