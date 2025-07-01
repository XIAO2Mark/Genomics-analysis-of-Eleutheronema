#!/bin/bash
# Olfactory Receptor Analysis Pipeline

WORKING_DIR="/Home/jiexiao/ORy"
OUTPUT_DIR="${WORKING_DIR}/olfactory_receptor_analysis"
GENOME_DIR="${WORKING_DIR}/genome_assemblies"
ZEBRAFISH_OR_DB="${WORKING_DIR}/databases/zebrafish_or_genes.fa"
THREADS=32

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Collect genomes from all target species
mkdir -p ${OUTPUT_DIR}/genomes

SPECIES=(
    "E_rhadinum"
    "E_tetradactylum"
    "D_rerio"
    ...
)

# Create symbolic links to genome files
for species in "${SPECIES[@]}"; do
    if [ -f "${GENOME_DIR}/${species}.fa" ]; then
        ln -sf "${GENOME_DIR}/${species}.fa" "${OUTPUT_DIR}/genomes/${species}.fa"
    else
        echo "Warning: Genome file for ${species} not found!"
    fi
done

# Identify putative OR genes using homology-based search with TBLASTN
mkdir -p ${OUTPUT_DIR}/tblastn_results

for species in "${SPECIES[@]}"; do
    if [ -f "${OUTPUT_DIR}/genomes/${species}.fa" ]; then
        echo "Processing ${species}..."
        
        makeblastdb -in "${OUTPUT_DIR}/genomes/${species}.fa" \
                    -dbtype nucl \
                    -out "${OUTPUT_DIR}/genomes/${species}_db"
        
        tblastn -query "${ZEBRAFISH_OR_DB}" \
                -db "${OUTPUT_DIR}/genomes/${species}_db" \
                -evalue 1e-10 \
                -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sframe" \
                -num_threads ${THREADS} \
                -out "${OUTPUT_DIR}/tblastn_results/${species}_tblastn.out"
    fi
done


# Multiple sequence alignment of OR proteins
mkdir -p ${OUTPUT_DIR}/alignments

# Combine all OR protein sequences
cat ${OUTPUT_DIR}/protein_sequences/*_OR_proteins.fa > ${OUTPUT_DIR}/all_OR_proteins.fa

# Perform multiple sequence alignment with MAFFT
mafft --auto --thread ${THREADS} ${OUTPUT_DIR}/all_OR_proteins.fa > ${OUTPUT_DIR}/alignments/all_OR_aligned.fa

# Phylogenetic tree construction
mkdir -p ${OUTPUT_DIR}/phylogeny

# Build maximum likelihood tree with IQ-TREE
iqtree -s ${OUTPUT_DIR}/alignments/all_OR_aligned.fa \
       -m MFP -bb 1000 -nt ${THREADS} \
       -pre ${OUTPUT_DIR}/phylogeny/OR_tree

