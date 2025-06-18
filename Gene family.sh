#!/bin/bash
# Gene Family Analysis Pipeline
# This script details the gene family expansion/contraction analysis across 17 fish species

# Set variables
WORKING_DIR="/path/to/working/directory"
PROTEIN_DIR="${WORKING_DIR}/proteins"
OUTPUT_DIR="${WORKING_DIR}/gene_family_analysis"
THREADS=32

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Step 1: Prepare protein sequences from 17 fish species
echo "Step 1: Preparing protein sequences"
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

# Get longest isoform for each gene
for species in "${SPECIES[@]}"; do
    echo "Processing ${species}..."
    
    # Extract the longest transcript for each gene
    python - << EOF
import sys
from Bio import SeqIO

def select_longest_transcripts(input_file, output_file):
    gene_to_transcript = {}
    
    # Read all sequences and track the longest transcript per gene
    for record in SeqIO.parse(input_file, "fasta"):
        # Parse gene ID from header (adjust pattern based on your data format)
        parts = record.id.split('|')
        if len(parts) > 1:
            gene_id = parts[0]
        else:
            gene_id = record.id.split('.')[0]  # Default fallback
        
        seq_len = len(record.seq)
        
        if gene_id not in gene_to_transcript or seq_len > gene_to_transcript[gene_id][1]:
            gene_to_transcript[gene_id] = (record, seq_len)
    
    # Write the longest transcripts to output file
    with open(output_file, 'w') as out_handle:
        for gene_id, (record, _) in gene_to_transcript.items():
            SeqIO.write(record, out_handle, "fasta")
    
    print(f"Processed {len(gene_to_transcript)} genes.")

select_longest_transcripts("${PROTEIN_DIR}/${species}.fa", "${OUTPUT_DIR}/protein_sequences/${species}_longest.fa")
EOF

done

# Combine all protein sequences for OrthoFinder
cat ${OUTPUT_DIR}/protein_sequences/*_longest.fa > ${OUTPUT_DIR}/protein_sequences/all_species_proteins.fa

# Step 2: All-against-all protein comparison using Diamond
echo "Step 2: All-against-all protein comparison"
mkdir -p ${OUTPUT_DIR}/diamond

# Create Diamond database
diamond makedb --in ${OUTPUT_DIR}/protein_sequences/all_species_proteins.fa \
              --db ${OUTPUT_DIR}/diamond/all_species_proteins

# Run Diamond for all-against-all comparison
diamond blastp --db ${OUTPUT_DIR}/diamond/all_species_proteins \
              --query ${OUTPUT_DIR}/protein_sequences/all_species_proteins.fa \
              --outfmt 6 --evalue 1e-5 --threads ${THREADS} \
              --out ${OUTPUT_DIR}/diamond/all_vs_all.blastp

# Step 3: Gene family clustering with OrthoFinder
echo "Step 3: Gene family clustering with OrthoFinder"
mkdir -p ${OUTPUT_DIR}/orthofinder

# Run OrthoFinder on the protein sequences
orthofinder -f ${OUTPUT_DIR}/protein_sequences/ \
           -t ${THREADS} \
           -o ${OUTPUT_DIR}/orthofinder \
           -S diamond

# Step 4: Extract single-copy orthologs for phylogenetic analysis
echo "Step 4: Extracting single-copy orthologs"
mkdir -p ${OUTPUT_DIR}/single_copy_orthologs

# Extract single-copy orthologs from OrthoFinder results
python - << EOF
import os
import sys
from Bio import SeqIO

# Path to OrthoFinder results
orthogroups_file = "${OUTPUT_DIR}/orthofinder/Results_*/Orthogroups/Orthogroups.tsv"
single_copy_file = "${OUTPUT_DIR}/orthofinder/Results_*/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"

# Read single-copy orthologs
with open(single_copy_file, 'r') as f:
    single_copy_ogs = set(line.strip() for line in f)

# Read orthogroups
og_to_genes = {}
with open(orthogroups_file, 'r') as f:
    header = f.readline().strip().split('\t')
    species_cols = {header[i]: i for i in range(1, len(header))}
    
    for line in f:
        parts = line.strip().split('\t')
        og_id = parts[0]
        if og_id in single_copy_ogs:
            og_to_genes[og_id] = {}
            for species, idx in species_cols.items():
                if idx < len(parts) and parts[idx]:
                    og_to_genes[og_id][species] = parts[idx]

# Extract sequences for each single-copy ortholog
for og_id, genes_by_species in og_to_genes.items():
    output_file = f"${OUTPUT_DIR}/single_copy_orthologs/{og_id}.fa"
    with open(output_file, 'w') as out_f:
        for species, gene_id in genes_by_species.items():
            species_file = f"${OUTPUT_DIR}/protein_sequences/{species}_longest.fa"
            for record in SeqIO.parse(species_file, "fasta"):
                if record.id == gene_id or record.id.startswith(gene_id + ".") or record.id.startswith(gene_id + "|"):
                    out_f.write(f">{species}\n{record.seq}\n")
                    break

print(f"Extracted sequences for {len(og_to_genes)} single-copy orthologs.")
EOF

# Step 5: Multiple sequence alignment and phylogenetic tree construction
echo "Step 5: Building phylogenetic tree"
mkdir -p ${OUTPUT_DIR}/alignments ${OUTPUT_DIR}/phylogeny

# Align each ortholog group
for og_file in ${OUTPUT_DIR}/single_copy_orthologs/*.fa; do
    og_id=$(basename ${og_file} .fa)
    echo "Aligning ${og_id}..."
    mafft --auto --thread ${THREADS} ${og_file} > ${OUTPUT_DIR}/alignments/${og_id}.aln
done

# Concatenate all alignments
python - << EOF
from Bio import AlignIO
import os

# Get all alignment files
aln_files = sorted([f for f in os.listdir("${OUTPUT_DIR}/alignments") if f.endswith(".aln")])

# Get all species from the first alignment
first_aln = AlignIO.read("${OUTPUT_DIR}/alignments/" + aln_files[0], "fasta")
species = [record.id for record in first_aln]

# Initialize concatenated sequences
concat_seqs = {sp: "" for sp in species}

# Concatenate alignments
for aln_file in aln_files:
    alignment = AlignIO.read("${OUTPUT_DIR}/alignments/" + aln_file, "fasta")
    for record in alignment:
        concat_seqs[record.id] += str(record.seq)

# Write concatenated alignment
with open("${OUTPUT_DIR}/phylogeny/concatenated.phy", "w") as out_f:
    out_f.write(f"{len(species)} {len(list(concat_seqs.values())[0])}\n")
    for sp, seq in concat_seqs.items():
        out_f.write(f"{sp} {seq}\n")

print(f"Concatenated {len(aln_files)} alignments into a matrix of {len(species)} species and {len(list(concat_seqs.values())[0])} sites.")
EOF

# Build maximum likelihood tree with RAxML
raxmlHPC-PTHREADS -T ${THREADS} -f a -x 12345 -p 12345 -N 1000 \
                  -m PROTGAMMAJTT -s ${OUTPUT_DIR}/phylogeny/concatenated.phy \
                  -n species_tree -w ${OUTPUT_DIR}/phylogeny/

# Step 6: Divergence time estimation
echo "Step 6: Divergence time estimation"
mkdir -p ${OUTPUT_DIR}/timetree

# Create MCMCTree control file
cat > ${OUTPUT_DIR}/timetree/mcmctree.ctl << EOF
seed = -1
seqfile = ${OUTPUT_DIR}/phylogeny/concatenated.phy
treefile = ${OUTPUT_DIR}/phylogeny/RAxML_bestTree.species_tree
outfile = ${OUTPUT_DIR}/timetree/mcmctree_output

ndata = 1
seqtype = 2    * 0: nucleotides; 1:codons; 2:AAs
usedata = 1    * 0: no data; 1:seq; 2:use approximate likelihood
clock = 2      * 1: global clock; 2: independent rates; 3: correlated rates
model = 0      * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:T92, 6:TN93, 7:REV
alpha = 0.5    * alpha parameter for gamma rates at sites
ncatG = 5      * No. categories in discrete gamma

cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?

BDparas = 1 1 0.1    * birth, death, sampling

kappa_gamma = 6 2     * gamma prior for kappa
alpha_gamma = 1 1     * gamma prior for alpha

rgene_gamma = 2 2     * gamma prior for overall rate
sigma2_gamma = 1 10   * gamma prior for sigma^2     (for clock=2 or 3)

finetune = 1: 0.1 0.1 0.1 0.01 0.5 0.5  * times, rates, mixing, paras, RateParas, FossilErr

print = 1
burnin = 20000
sampfreq = 100
nsample = 20000
EOF

# Add calibration points based on TimeTree (modify as needed)
# Format calibration points in the tree file according to MCMCTree format
# Example: ((A, B)@0.8, C);  where 0.8 is the calibration time

# Run MCMCTree for divergence time estimation
cd ${OUTPUT_DIR}/timetree
mcmctree mcmctree.ctl

# Additionally use R8S for cross-validation
# Prepare r8s input file
cat > ${OUTPUT_DIR}/timetree/r8s_input.txt << EOF
#NEXUS
BEGIN TAXA;
DIMENSIONS NTAX=17;
TAXLABELS
E_rhadinum E_tetradactylum D_rerio B_splendens L_calcarifer C_melampygus E_naucrates S_dumerili C_semilaevis O_latipes T_jaculatrix M_salmoides E_lanceolatus H_stenolepis P_olivaceus C_lumpus L_oculatus
;
END;
BEGIN TREES;
TREE tree1 = [RAxML tree in Newick format with calibration points]
END;
BEGIN r8s;
blformat lengths=persite nsites=1000 ultrametric=no;
collapse;
fixage taxon=L_oculatus age=350;  # Example calibration based on TimeTree
divtime method=PL algorithm=TN crossv=yes fossilconstrained=yes;
showage;
END;
EOF

# Run r8s
r8s -b -f ${OUTPUT_DIR}/timetree/r8s_input.txt > ${OUTPUT_DIR}/timetree/r8s_output.log

# Step 7: Gene family expansion/contraction analysis with CAFE
echo "Step 7: Gene family expansion/contraction analysis"
mkdir -p ${OUTPUT_DIR}/cafe

# Create CAFE input file from OrthoFinder results
python - << EOF
import os

# Path to OrthoFinder results
orthogroups_file = "${OUTPUT_DIR}/orthofinder/Results_*/Orthogroups/Orthogroups.tsv"
tree_file = "${OUTPUT_DIR}/timetree/mcmctree_output.tre"  # Using MCMCTree output

# Read the tree file (simplified for this script)
with open(tree_file, 'r') as f:
    tree_line = f.readlines()[-1].strip()

# Read orthogroups and count genes per species
with open(orthogroups_file, 'r') as f:
    header = f.readline().strip().split('\t')
    species_cols = {header[i]: i for i in range(1, len(header))}
    
    # Open output file for CAFE
    with open("${OUTPUT_DIR}/cafe/cafe_input.txt", 'w') as out_f:
        # Write header
        out_f.write("Desc\tFamily ID\t" + "\t".join([sp for sp in header[1:]]) + "\n")
        
        # Write gene counts per family
        family_id = 1
        for line in f:
            parts = line.strip().split('\t')
            og_id = parts[0]
            counts = []
            
            for species in header[1:]:
                idx = species_cols[species]
                if idx < len(parts) and parts[idx]:
                    # Count genes in this species for this orthogroup
                    gene_count = len(parts[idx].split(', '))
                else:
                    gene_count = 0
                counts.append(str(gene_count))
            
            out_f.write(f"{og_id}\tFamily{family_id}\t" + "\t".join(counts) + "\n")
            family_id += 1
EOF

# Create CAFE control file
cat > ${OUTPUT_DIR}/cafe/cafe.ctl << EOF
#CAFE input file
tree (((((((((E_rhadinum:_,E_tetradactylum:_):_,L_calcarifer:_):_,(((C_melampygus:_,E_naucrates:_):_,S_dumerili:_):_,T_jaculatrix:_):_):_,M_salmoides:_):_,E_lanceolatus:_):_,B_splendens:_):_,((C_semilaevis:_,P_olivaceus:_):_,H_stenolepis:_):_):_,C_lumpus:_):_,O_latipes:_):_,D_rerio:_,L_oculatus:_);

# Lambda (family-wide rate of gene gain/loss)
lambda -s

# Input data file
load -i ${OUTPUT_DIR}/cafe/cafe_input.txt -p 0.05 -t 10

# Report results
report ${OUTPUT_DIR}/cafe/cafe_results
EOF

# Update the tree in the control file with the actual divergence times
python - << EOF
import re

# Read the time tree file
with open("${OUTPUT_DIR}/timetree/mcmctree_output.tre", 'r') as f:
    time_tree = f.read().strip()

# Read CAFE control file
with open("${OUTPUT_DIR}/cafe/cafe.ctl", 'r') as f:
    cafe_ctl = f.read()

# Extract the tree structure from CAFE file
tree_pattern = r'tree \(.*?\);'
tree_match = re.search(tree_pattern, cafe_ctl, re.DOTALL)
if tree_match:
    cafe_tree = tree_match.group(0)
    
    # Replace placeholder branch lengths with divergence times
    time_tree = time_tree.replace(';', '')  # Remove semicolon
    species = ["E_rhadinum", "E_tetradactylum", "L_calcarifer", "C_melampygus", 
               "E_naucrates", "S_dumerili", "T_jaculatrix", "M_salmoides", 
               "E_lanceolatus", "B_splendens", "C_semilaevis", "P_olivaceus", 
               "H_stenolepis", "C_lumpus", "O_latipes", "D_rerio", "L_oculatus"]
    
    # Parse divergence times and update CAFE tree
    updated_tree = cafe_tree
    for sp in species:
        for match in re.finditer(f'{sp}:_', updated_tree):
            # Find the corresponding time in time_tree
            sp_pattern = f'{sp}:(\\d+\\.\\d+)'
            sp_time_match = re.search(sp_pattern, time_tree)
            if sp_time_match:
                time_value = sp_time_match.group(1)
                updated_tree = updated_tree.replace(f'{sp}:_', f'{sp}:{time_value}')
    
    # Update internal branches
    internal_pattern = r'):_'
    for match in re.finditer(internal_pattern, updated_tree):
        # We need a more sophisticated approach here to match internal nodes
        # This is a simplified placeholder
        updated_tree = updated_tree.replace('):_', '):0.1')
    
    # Replace tree in CAFE control file
    updated_cafe_ctl = cafe_ctl.replace(cafe_tree, updated_tree)
    
    # Write updated CAFE control file
    with open("${OUTPUT_DIR}/cafe/cafe.ctl", 'w') as f:
        f.write(updated_cafe_ctl)
EOF

# Run CAFE analysis
cd ${OUTPUT_DIR}/cafe
cafe cafe.ctl

# Step 8: Analyze CAFE results and perform enrichment analysis
echo "Step 8: Analyzing gene family evolution results"
mkdir -p ${OUTPUT_DIR}/cafe/enrichment

# Extract significantly expanded/contracted gene families
python - << EOF
import os
import re

# Parse CAFE results
results_file = "${OUTPUT_DIR}/cafe/cafe_results.cafe"
expanded_families = {}
contracted_families = {}

with open(results_file, 'r') as f:
    for line in f:
        if line.startswith('Final'):
            parts = line.strip().split('\t')
            family_id = parts[1]
            pvalue = float(parts[2])
            
            if pvalue <= 0.05:  # Significant change
                # Extract species-specific expansions/contractions
                for species_info in parts[3:]:
                    if '<' in species_info or '>' in species_info:
                        species, change = species_info.split('(')
                        species = species.strip()
                        
                        # Extract direction and magnitude
                        if '+' in change:  # Expansion
                            if species not in expanded_families:
                                expanded_families[species] = []
                            expanded_families[species].append(family_id)
                        elif '-' in change:  # Contraction
                            if species not in contracted_families:
                                contracted_families[species] = []
                            contracted_families[species].append(family_id)

# Write expanded/contracted families for each species
for species, families in expanded_families.items():
    with open(f"${OUTPUT_DIR}/cafe/enrichment/{species}_expanded.txt", 'w') as f:
        for family in families:
            f.write(f"{family}\n")

for species, families in contracted_families.items():
    with open(f"${OUTPUT_DIR}/cafe/enrichment/{species}_contracted.txt", 'w') as f:
        for family in families:
            f.write(f"{family}\n")

# Focus on E. rhadinum
er_expanded = expanded_families.get('E_rhadinum', [])
er_contracted = contracted_families.get('E_rhadinum', [])

print(f"E. rhadinum has {len(er_expanded)} significantly expanded gene families and {len(er_contracted)} significantly contracted gene families.")
EOF

# Perform GO and KEGG enrichment analysis for E. rhadinum
mkdir -p ${OUTPUT_DIR}/cafe/enrichment/E_rhadinum

# Extract gene IDs for the significantly changed families
python - << EOF
import os

# Path to OrthoFinder results and gene family lists
orthogroups_file = "${OUTPUT_DIR}/orthofinder/Results_*/Orthogroups/Orthogroups.tsv"
expanded_file = "${OUTPUT_DIR}/cafe/enrichment/E_rhadinum_expanded.txt"
contracted_file = "${OUTPUT_DIR}/cafe/enrichment/E_rhadinum_contracted.txt"

# Read expanded and contracted family IDs
expanded_families = []
if os.path.exists(expanded_file):
    with open(expanded_file, 'r') as f:
        expanded_families = [line.strip() for line in f]

contracted_families = []
if os.path.exists(contracted_file):
    with open(contracted_file, 'r') as f:
        contracted_families = [line.strip() for line in f]

# Map family IDs to OrthoFinder OG IDs
family_to_og = {}
with open("${OUTPUT_DIR}/cafe/cafe_input.txt", 'r') as f:
    next(f)  # Skip header
    for line in f:
        parts = line.strip().split('\t')
        og_id = parts[0]
        family_id = parts[1]
        family_to_og[family_id] = og_id

# Extract genes from expanded and contracted families
expanded_genes = []
contracted_genes = []

with open(orthogroups_file, 'r') as f:
    header = f.readline().strip().split('\t')
    er_idx = header.index('E_rhadinum')
    
    for line in f:
        parts = line.strip().split('\t')
        og_id = parts[0]
        
        # Check if this OG belongs to a significantly changed family
        for family_id, og in family_to_og.items():
            if og == og_id:
                if family_id in expanded_families and er_idx < len(parts):
                    genes = parts[er_idx].split(', ')
                    expanded_genes.extend(genes)
                elif family_id in contracted_families and er_idx < len(parts):
                    genes = parts[er_idx].split(', ')
                    contracted_genes.extend(genes)
                break

# Write gene lists for enrichment analysis
with open("${OUTPUT_DIR}/cafe/enrichment/E_rhadinum/expanded_genes.txt", 'w') as f:
    for gene in expanded_genes:
        f.write(f"{gene}\n")

with open("${OUTPUT_DIR}/cafe/enrichment/E_rhadinum/contracted_genes.txt", 'w') as f:
    for gene in contracted_genes:
        f.write(f"{gene}\n")
EOF

# Run GO enrichment analysis
Rscript - << EOF
library(clusterProfiler)
library(org.Dr.eg.db)  # Using zebrafish as reference (adjust as needed)

# Read gene lists
expanded_genes <- readLines("${OUTPUT_DIR}/cafe/enrichment/E_rhadinum/expanded_genes.txt")
contracted_genes <- readLines("${OUTPUT_DIR}/cafe/enrichment/E_rhadinum/contracted_genes.txt")

# Read all genes as background
all_genes <- unique(readLines("${OUTPUT_DIR}/maker/erhadinum_genes.list"))

# Map to gene symbols (this would need adjustment for your actual gene IDs)
# This is an example assuming zebrafish homology mapping exists
gene_mapping <- read.delim("${WORKING_DIR}/gene_mappings/erhadinum_to_zebrafish.tsv", 
                           header=TRUE, stringsAsFactors=FALSE)

expanded_entrez <- gene_mapping$zebrafish_entrez[gene_mapping$erhadinum_gene %in% expanded_genes]
contracted_entrez <- gene_mapping$zebrafish_entrez[gene_mapping$erhadinum_gene %in% contracted_genes]
all_entrez <- gene_mapping$zebrafish_entrez[gene_mapping$erhadinum_gene %in% all_genes]

# GO enrichment for expanded genes
if(length(expanded_entrez) > 0) {
  ego_expanded <- enrichGO(gene = expanded_entrez,
                         universe = all_entrez,
                         OrgDb = org.Dr.eg.db,
                         ont = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
  
  write.csv(as.data.frame(ego_expanded), 
            "${OUTPUT_DIR}/cafe/enrichment/E_rhadinum/expanded_GO_enrichment.csv")
}

# GO enrichment for contracted genes
if(length(contracted_entrez) > 0) {
  ego_contracted <- enrichGO(gene = contracted_entrez,
                           universe = all_entrez,
                           OrgDb = org.Dr.eg.db,
                           ont = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
  
  write.csv(as.data.frame(ego_contracted), 
            "${OUTPUT_DIR}/cafe/enrichment/E_rhadinum/contracted_GO_enrichment.csv")
}

# KEGG pathway enrichment for expanded genes
if(length(expanded_entrez) > 0) {
  ekegg_expanded <- enrichKEGG(gene = expanded_entrez,
                             organism = 'dre',  # zebrafish KEGG organism code
                             universe = all_entrez,
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05)
  
  write.csv(as.data.frame(ekegg_expanded), 
            "${OUTPUT_DIR}/cafe/enrichment/E_rhadinum/expanded_KEGG_enrichment.csv")
}

# KEGG pathway enrichment for contracted genes
if(length(contracted_entrez) > 0) {
  ekegg_contracted <- enrichKEGG(gene = contracted_entrez,
                               organism = 'dre',  # zebrafish KEGG organism code
                               universe = all_entrez,
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05)
  
  write.csv(as.data.frame(ekegg_contracted), 
            "${OUTPUT_DIR}/cafe/enrichment/E_rhadinum/contracted_KEGG_enrichment.csv")
}
EOF

echo "Gene family analysis pipeline completed!"
