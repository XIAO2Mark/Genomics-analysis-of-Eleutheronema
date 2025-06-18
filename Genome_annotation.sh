#!/bin/bash
# Genome Annotation Pipeline for E. rhadinum
# This script details the annotation process, including repeat identification and gene prediction

# Set variables
WORKING_DIR="/path/to/working/directory"
GENOME_FASTA="${WORKING_DIR}/assembly/scaffolding/final_assembly.fasta"
RNA_SEQ_DIR="${WORKING_DIR}/transcriptome_data"
PROTEIN_DB="${WORKING_DIR}/databases/protein_db"
OUTPUT_DIR="${WORKING_DIR}/annotation"
THREADS=32

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Step 1: Repeat identification and masking
echo "Step 1: Repeat identification and masking"
mkdir -p ${OUTPUT_DIR}/repeats

# De novo repeat identification
## Build repeat library using RepeatModeler
BuildDatabase -name erhadinum_db -engine ncbi ${GENOME_FASTA}
RepeatModeler -database erhadinum_db -engine ncbi -pa ${THREADS} -LTRStruct

# Combine with known repeats from RepBase
cat RepeatModeler-*/consensi.fa.classified ${WORKING_DIR}/databases/RepBase/RepBase_fish.fa > ${OUTPUT_DIR}/repeats/repeat_library.fa

# Mask repeats using RepeatMasker
RepeatMasker -pa ${THREADS} -lib ${OUTPUT_DIR}/repeats/repeat_library.fa \
             -xsmall -gff -html -dir ${OUTPUT_DIR}/repeats/ ${GENOME_FASTA}

# Generate statistics
perl ${WORKING_DIR}/tools/RepeatMasker/util/buildSummary.pl -species fish ${OUTPUT_DIR}/repeats/final_assembly.fasta.out > ${OUTPUT_DIR}/repeats/repeat_summary.txt

# Step 2: Transcriptome assembly for gene prediction support
echo "Step 2: Transcriptome assembly"
mkdir -p ${OUTPUT_DIR}/transcriptome

# Map RNA-seq reads to the genome
hisat2-build ${GENOME_FASTA} ${OUTPUT_DIR}/transcriptome/genome_idx
hisat2 -p ${THREADS} -x ${OUTPUT_DIR}/transcriptome/genome_idx \
       -1 ${RNA_SEQ_DIR}/RNA_R1.fastq.gz -2 ${RNA_SEQ_DIR}/RNA_R2.fastq.gz \
       | samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/transcriptome/rnaseq.sorted.bam -
samtools index ${OUTPUT_DIR}/transcriptome/rnaseq.sorted.bam

# Assemble transcripts with StringTie
stringtie ${OUTPUT_DIR}/transcriptome/rnaseq.sorted.bam \
          -o ${OUTPUT_DIR}/transcriptome/stringtie_assembly.gtf \
          -p ${THREADS}

# Convert to FASTA for MAKER input
gffread -w ${OUTPUT_DIR}/transcriptome/transcripts.fa \
        -g ${GENOME_FASTA} ${OUTPUT_DIR}/transcriptome/stringtie_assembly.gtf

# Step 3: Protein-coding gene prediction using MAKER
echo "Step 3: Protein-coding gene prediction using MAKER"
mkdir -p ${OUTPUT_DIR}/maker

# Create MAKER control files
cd ${OUTPUT_DIR}/maker
maker -CTL

# Edit the maker_opts.ctl file with the correct paths and parameters
cat > maker_opts.ctl << EOF
#-----Genome (these are always required)
genome=${GENOME_FASTA}
organism_type=eukaryotic
#-----Re-annotation Using MAKER Derived GFF3
maker_gff=
est_pass=0
altest_pass=0
protein_pass=0
rm_pass=0
#-----EST Evidence (for best results provide a file for at least one)
est=${OUTPUT_DIR}/transcriptome/transcripts.fa
altest=
est_gff=
altest_gff=
#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=${PROTEIN_DB}/zebrafish_proteins.fa,${PROTEIN_DB}/uniprot_fish.fa
protein_gff=
#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=fish
rmlib=${OUTPUT_DIR}/repeats/repeat_library.fa
repeat_protein=${WORKING_DIR}/tools/RepeatMasker/Libraries/RepeatPeps.lib
rm_gff=${OUTPUT_DIR}/repeats/final_assembly.fasta.out.gff
softmask=1
#-----Gene Prediction
snaphmm=${WORKING_DIR}/tools/snap/fish.hmm
gmhmm=${WORKING_DIR}/tools/GeneMark/gmhmm_fish.mod
augustus_species=zebrafish
fgenesh_par_file=${WORKING_DIR}/tools/fgenesh/fish.par
pred_gff=
model_gff=
run_evm=1
est2genome=1
protein2genome=1
trna=1
#-----MAKER Behavior Options
max_dna_len=300000
min_contig=1000
keep_preds=1
split_hit=10000
single_exon=1
correct_est_fusion=1
EOF

# Run MAKER annotation pipeline
mpiexec -n ${THREADS} maker -base erhadinum_maker

# Process MAKER output
cd ${OUTPUT_DIR}/maker/erhadinum_maker.maker.output
gff3_merge -d erhadinum_maker_master_datastore_index.log -o erhadinum_genes.gff
fasta_merge -d erhadinum_maker_master_datastore_index.log -o ${OUTPUT_DIR}/maker/erhadinum_proteins.fa

# Extract gene and transcript sequences
${WORKING_DIR}/tools/maker/bin/gff3_to_fasta ${GENOME_FASTA} erhadinum_genes.gff \
    -t cds -o ${OUTPUT_DIR}/maker/erhadinum_cds.fa
${WORKING_DIR}/tools/maker/bin/gff3_to_fasta ${GENOME_FASTA} erhadinum_genes.gff \
    -t mrna -o ${OUTPUT_DIR}/maker/erhadinum_transcripts.fa

# Step 4: Non-coding RNA prediction
echo "Step 4: Non-coding RNA prediction"
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

# Other ncRNA prediction with Infernal and Rfam
cmscan --cpu ${THREADS} --tblout ${OUTPUT_DIR}/ncrna/rfam.tblout \
       --fmt 2 --clanin ${WORKING_DIR}/databases/Rfam/Rfam.clanin \
       ${WORKING_DIR}/databases/Rfam/Rfam.cm ${GENOME_FASTA}

# Convert to GFF3 format
${WORKING_DIR}/tools/infernal/cmsscan2gff.pl ${OUTPUT_DIR}/ncrna/rfam.tblout > ${OUTPUT_DIR}/ncrna/rfam.gff

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

echo "Genome annotation pipeline completed!"