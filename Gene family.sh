#!/bin/bash
# Gene Family Analysis Pipeline

WORKING_DIR="/Home/jiexiao/gene_family"
PROTEIN_DIR="${WORKING_DIR}/proteins"
OUTPUT_DIR="${WORKING_DIR}/gene_family_analysis"
THREADS=32

mkdir -p ${OUTPUT_DIR}

# Step 1: Prepare protein sequences from 17 fish species
mkdir -p ${OUTPUT_DIR}/protein_sequences

# List of species to be analyzed
SPECIES=(
    "E_rhadinum"
    "E_tetradactylum"
    "D_rerio"
    "B_splendens"
    "L_calcarifer"
    "C_melampygus"
    "E_naucrates"
    "S_dumerili"
    "C_semilaevis"
    "O_latipes"
    "T_jaculatrix"
    "M_salmoides"
    "E_lanceolatus"
    "H_stenolepis"
    "P_olivaceus"
    "C_lumpus"
    "L_oculatus"
)

# Combine all protein sequences for OrthoFinder
cat ${OUTPUT_DIR}/protein_sequences/*_longest.fa > ${OUTPUT_DIR}/protein_sequences/all_species_proteins.fa

# All-against-all protein comparison using Diamond
mkdir -p ${OUTPUT_DIR}/diamond
diamond makedb --in ${OUTPUT_DIR}/protein_sequences/all_species_proteins.fa \
              --db ${OUTPUT_DIR}/diamond/all_species_proteins

diamond blastp --db ${OUTPUT_DIR}/diamond/all_species_proteins \
              --query ${OUTPUT_DIR}/protein_sequences/all_species_proteins.fa \
              --outfmt 6 --evalue 1e-5 --threads ${THREADS} \
              --out ${OUTPUT_DIR}/diamond/all_vs_all.blastp

# Gene family clustering with OrthoFinder
mkdir -p ${OUTPUT_DIR}/orthofinder
# Run OrthoFinder on the protein sequences
orthofinder -f ${OUTPUT_DIR}/protein_sequences/ \
           -t ${THREADS} \
           -o ${OUTPUT_DIR}/orthofinder \
           -S diamond


# Build maximum likelihood tree with RAxML
raxmlHPC-PTHREADS -T ${THREADS} -f a -x 12345 -p 12345 -N 1000 \
                  -m PROTGAMMAJTT -s ${OUTPUT_DIR}/phylogeny/concatenated.phy \
                  -n species_tree -w ${OUTPUT_DIR}/phylogeny/

# Divergence time estimation
mkdir -p ${OUTPUT_DIR}/timetree
cd ${OUTPUT_DIR}/timetree
mcmctree mcmctree.ctl

# Run r8s
r8s -b -f ${OUTPUT_DIR}/timetree/r8s_input.txt > ${OUTPUT_DIR}/timetree/r8s_output.log

# Gene family expansion/contraction analysis with CAFE
mkdir -p ${OUTPUT_DIR}/cafe
cd ${OUTPUT_DIR}/cafe
cafe cafe.ctl

