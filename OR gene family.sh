#!/bin/bash
# Olfactory Receptor Analysis Pipeline
# This script details the identification, classification, and phylogenetic analysis of olfactory receptor genes

# Set variables
WORKING_DIR="/path/to/working/directory"
OUTPUT_DIR="${WORKING_DIR}/olfactory_receptor_analysis"
GENOME_DIR="${WORKING_DIR}/genome_assemblies"
ZEBRAFISH_OR_DB="${WORKING_DIR}/databases/zebrafish_or_genes.fa"
THREADS=32

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Step 1: Collect genomes from all target species
echo "Step 1: Collecting genome sequences"
mkdir -p ${OUTPUT_DIR}/genomes

# List of fish species to analyze
SPECIES=(
    "E_rhadinum"
    "E_tetradactylum"
    # Add the other 47 fish genomes here
    "D_rerio"
    # ... additional species ...
)

# Create symbolic links to genome files
for species in "${SPECIES[@]}"; do
    if [ -f "${GENOME_DIR}/${species}.fa" ]; then
        ln -sf "${GENOME_DIR}/${species}.fa" "${OUTPUT_DIR}/genomes/${species}.fa"
    else
        echo "Warning: Genome file for ${species} not found!"
    fi
done

# Step 2: Identify putative OR genes using homology-based search with TBLASTN
echo "Step 2: Identifying putative OR genes"
mkdir -p ${OUTPUT_DIR}/tblastn_results

# Run TBLASTN for each species
for species in "${SPECIES[@]}"; do
    if [ -f "${OUTPUT_DIR}/genomes/${species}.fa" ]; then
        echo "Processing ${species}..."
        
        # Create BLAST database for the genome
        makeblastdb -in "${OUTPUT_DIR}/genomes/${species}.fa" \
                    -dbtype nucl \
                    -out "${OUTPUT_DIR}/genomes/${species}_db"
        
        # Run TBLASTN search
        tblastn -query "${ZEBRAFISH_OR_DB}" \
                -db "${OUTPUT_DIR}/genomes/${species}_db" \
                -evalue 1e-10 \
                -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sframe" \
                -num_threads ${THREADS} \
                -out "${OUTPUT_DIR}/tblastn_results/${species}_tblastn.out"
    fi
done

# Step 3: Extract and curate OR gene sequences
echo "Step 3: Extracting and curating OR sequences"
mkdir -p ${OUTPUT_DIR}/extracted_sequences

# Extract OR gene sequences for each species
for species in "${SPECIES[@]}"; do
    if [ -f "${OUTPUT_DIR}/tblastn_results/${species}_tblastn.out" ]; then
        echo "Extracting OR sequences for ${species}..."
        
        # Extract regions with high sequence similarity to OR genes
        python - << EOF
import sys
from Bio import SeqIO
from Bio.Seq import Seq

# Read TBLASTN results
blast_results = []
with open("${OUTPUT_DIR}/tblastn_results/${species}_tblastn.out", 'r') as f:
    for line in f:
        fields = line.strip().split('\t')
        query_id = fields[0]
        subject_id = fields[1]
        pident = float(fields[2])
        sstart = int(fields[8])
        send = int(fields[9])
        evalue = float(fields[10])
        sframe = int(fields[12])
        
        # Filter by E-value and percent identity
        if evalue <= 1e-10 and pident >= 50:
            # Ensure sstart < send
            if sstart > send:
                sstart, send = send, sstart
            
            # Add padding around the match
            sstart = max(1, sstart - 100)
            send = send + 100
            
            blast_results.append({
                'query_id': query_id,
                'subject_id': subject_id,
                'sstart': sstart,
                'send': send,
                'sframe': sframe
            })

# Read genome sequence
genome_dict = {record.id: record for record in SeqIO.parse("${OUTPUT_DIR}/genomes/${species}.fa", "fasta")}

# Extract putative OR gene sequences
count = 0
with open("${OUTPUT_DIR}/extracted_sequences/${species}_putative_ORs.fa", 'w') as out_f:
    for result in blast_results:
        subj_id = result['subject_id']
        if subj_id in genome_dict:
            seq_record = genome_dict[subj_id]
            
            # Extract the region
            start = result['sstart'] - 1  # 0-based indexing
            end = result['send']
            extracted_seq = seq_record.seq[start:end]
            
            # Adjust for frame if necessary
            frame = result['sframe']
            if frame < 0:
                extracted_seq = extracted_seq.reverse_complement()
            
            # Write to output file
            out_f.write(f">{species}_OR_{count}|{subj_id}:{start+1}-{end}|frame={frame}\n")
            out_f.write(f"{extracted_seq}\n")
            count += 1

print(f"Extracted {count} putative OR gene sequences for {species}.")
EOF
    fi
done

# Step 4: Identify coding sequences and translate to protein
echo "Step 4: Identifying coding sequences and translating to protein"
mkdir -p ${OUTPUT_DIR}/coding_sequences ${OUTPUT_DIR}/protein_sequences

# Process each species
for species in "${SPECIES[@]}"; do
    if [ -f "${OUTPUT_DIR}/extracted_sequences/${species}_putative_ORs.fa" ]; then
        echo "Identifying coding sequences for ${species}..."
        
        # Identify open reading frames and translate to protein
        python - << EOF
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def find_orf(seq, min_length=300):
    # Look for start codon (ATG) and a downstream stop codon
    # Returns the longest ORF that exceeds min_length
    starts = []
    for i in range(0, len(seq) - 2, 3):
        if seq[i:i+3] == 'ATG':
            starts.append(i)
    
    best_orf = None
    max_length = 0
    
    for start in starts:
        for i in range(start, len(seq) - 2, 3):
            if seq[i:i+3] in ['TAA', 'TAG', 'TGA']:  # Stop codon
                orf_length = i - start + 3
                if orf_length > max_length and orf_length >= min_length:
                    max_length = orf_length
                    best_orf = (start, i + 3)
                break
    
    return best_orf

# Process putative OR sequences
count = 0
valid_count = 0
with open("${OUTPUT_DIR}/coding_sequences/${species}_OR_cds.fa", 'w') as cds_out, \
     open("${OUTPUT_DIR}/protein_sequences/${species}_OR_proteins.fa", 'w') as prot_out:
    
    for record in SeqIO.parse("${OUTPUT_DIR}/extracted_sequences/${species}_putative_ORs.fa", "fasta"):
        count += 1
        
        # Try all reading frames
        best_orf = None
        max_length = 0
        best_frame = 0
        best_prot = None
        
        for frame in range(3):
            seq = record.seq[frame:]
            orf = find_orf(str(seq))
            if orf and (orf[1] - orf[0]) > max_length:
                max_length = orf[1] - orf[0]
                best_orf = orf
                best_frame = frame
                
                # Translate the ORF
                cds = seq[orf[0]:orf[1]]
                prot = cds.translate(to_stop=True)
                if len(prot) >= 100:  # Minimum protein length
                    best_prot = prot
        
        # Check reverse complement as well
        rc_seq = record.seq.reverse_complement()
        for frame in range(3):
            seq = rc_seq[frame:]
            orf = find_orf(str(seq))
            if orf and (orf[1] - orf[0]) > max_length:
                max_length = orf[1] - orf[0]
                best_orf = orf
                best_frame = frame + 3  # Indicate reverse frame
                
                # Translate the ORF
                cds = seq[orf[0]:orf[1]]
                prot = cds.translate(to_stop=True)
                if len(prot) >= 100:  # Minimum protein length
                    best_prot = prot
        
        # If valid ORF found, write to output
        if best_orf and best_prot:
            valid_count += 1
            record_id = record.id.split('|')[0]
            
            if best_frame < 3:  # Forward frame
                seq = record.seq[best_frame:]
                cds = seq[best_orf[0]:best_orf[1]]
            else:  # Reverse frame
                seq = rc_seq[best_frame - 3:]
                cds = seq[best_orf[0]:best_orf[1]]
            
            # Write CDS and protein sequences
            cds_out.write(f">{record_id}_CDS\n{cds}\n")
            prot_out.write(f">{record_id}_PROT\n{best_prot}\n")

print(f"Processed {count} sequences for ${species}, found {valid_count} valid OR coding sequences.")
EOF
    fi
done

# Step 5: Multiple sequence alignment of OR proteins
echo "Step 5: Multiple sequence alignment of OR proteins"
mkdir -p ${OUTPUT_DIR}/alignments

# Combine all OR protein sequences
cat ${OUTPUT_DIR}/protein_sequences/*_OR_proteins.fa > ${OUTPUT_DIR}/all_OR_proteins.fa

# Perform multiple sequence alignment with MAFFT
mafft --auto --thread ${THREADS} ${OUTPUT_DIR}/all_OR_proteins.fa > ${OUTPUT_DIR}/alignments/all_OR_aligned.fa

# Step 6: Phylogenetic tree construction
echo "Step 6: Constructing phylogenetic tree"
mkdir -p ${OUTPUT_DIR}/phylogeny

# Build maximum likelihood tree with IQ-TREE
iqtree -s ${OUTPUT_DIR}/alignments/all_OR_aligned.fa \
       -m MFP -bb 1000 -nt ${THREADS} \
       -pre ${OUTPUT_DIR}/phylogeny/OR_tree

# Step 7: Classification of OR genes into subfamilies
echo "Step 7: Classifying OR genes into subfamilies"
mkdir -p ${OUTPUT_DIR}/classification

# Use known zebrafish OR genes for classification
python - << EOF
import sys
from Bio import SeqIO
from Bio import AlignIO
import pandas as pd
import numpy as np
from ete3 import Tree

# Load the phylogenetic tree
tree = Tree("${OUTPUT_DIR}/phylogeny/OR_tree.treefile", format=1)

# Load known zebrafish OR classification
zebrafish_or_classification = {}
with open("${WORKING_DIR}/databases/zebrafish_or_classification.tsv", 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            gene_id = parts[0]
            subfamily = parts[1]
            zebrafish_or_classification[gene_id] = subfamily

# Function to find the closest zebrafish OR in the tree
def find_closest_zebrafish_or(node, classified_genes):
    # If this is a leaf node
    if node.is_leaf():
        if node.name in classified_genes:
            return node.name, 0
        return None, float('inf')
    
    # Check children
    best_match = None
    min_dist = float('inf')
    
    for child in node.children:
        match, dist = find_closest_zebrafish_or(child, classified_genes)
        if match and dist < min_dist:
            best_match = match
            min_dist = dist
    
    # Add distance to this node
    if best_match:
        return best_match, min_dist + 1
    
    return None, float('inf')

# Classify OR genes based on tree topology and sequence identity
classifications = {}
for record in SeqIO.parse("${OUTPUT_DIR}/protein_sequences/*_OR_proteins.fa", "fasta"):
    gene_id = record.id
    
    # Skip if it's a zebrafish OR (already classified)
    if gene_id.startswith("D_rerio"):
        if gene_id in zebrafish_or_classification:
            classifications[gene_id] = zebrafish_or_classification[gene_id]
        continue
    
    # Find the node in the tree
    node = None
    for leaf in tree.get_leaves():
        if leaf.name == gene_id:
            node = leaf
            break
    
    if node:
        # Find closest classified zebrafish OR
        closest_zf, _ = find_closest_zebrafish_or(node, zebrafish_or_classification)
        if closest_zf:
            classifications[gene_id] = zebrafish_or_classification[closest_zf]
        else:
            # If no close zebrafish OR found, use sequence identity
            max_identity = 0
            best_subfamily = None
            
            for zf_gene, zf_subfamily in zebrafish_or_classification.items():
                # Calculate sequence identity (would require pairwise alignment)
                # This is a placeholder - in a real script, you would calculate actual identity
                identity = 0.5  # Placeholder value
                
                if identity > max_identity:
                    max_identity = identity
                    best_subfamily = zf_subfamily
            
            if best_subfamily:
                classifications[gene_id] = best_subfamily
            else:
                classifications[gene_id] = "Unclassified"
    else:
        classifications[gene_id] = "Unclassified"

# Write classifications to file
with open("${OUTPUT_DIR}/classification/OR_subfamily_classification.tsv", 'w') as f:
    f.write("Gene_ID\tSubfamily\n")
    for gene_id, subfamily in classifications.items():
        f.write(f"{gene_id}\t{subfamily}\n")

# Count OR genes by species and subfamily
counts = {}
for gene_id, subfamily in classifications.items():
    species = gene_id.split('_')[0]
    if species not in counts:
        counts[species] = {}
    
    if subfamily not in counts[species]:
        counts[species][subfamily] = 0
    
    counts[species][subfamily] += 1

# Write counts to file
with open("${OUTPUT_DIR}/classification/OR_counts_by_species_subfamily.tsv", 'w') as f:
    # Get all unique subfamilies
    all_subfamilies = sorted(set(subfamily for species_data in counts.values() for subfamily in species_data))
    
    # Write header
    f.write("Species\t" + "\t".join(all_subfamilies) + "\tTotal\n")
    
    # Write counts for each species
    for species in sorted(counts.keys()):
        row = [species]
        total = 0
        
        for subfamily in all_subfamilies:
            count = counts[species].get(subfamily, 0)
            row.append(str(count))
            total += count
        
        row.append(str(total))
        f.write("\t".join(row) + "\n")
EOF

