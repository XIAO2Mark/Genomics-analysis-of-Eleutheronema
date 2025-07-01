#!/bin/bash

REFERENCE_GENOME="E_rhadinum_reference.fasta"
OUTGROUP_SRA="SRR11793813"
THREADS=16
GENERATION_TIME=4.5
MUTATION_RATE=1.9e-8

mkdir -p {raw_data,aligned,variants,results/{admixture,phylogeny,pca,ld_decay,psmc,fst,introgression}}

# DATA PREPARATION AND ALIGNMENT

fastq-dump --split-files ${OUTGROUP_SRA} -O raw_data/

# Index reference genome
bwa index ${REFERENCE_GENOME}
samtools faidx ${REFERENCE_GENOME}
gatk CreateSequenceDictionary -R ${REFERENCE_GENOME}

# Align reads using BWA
for sample in raw_data/*_1.fastq; do
    base=$(basename ${sample} _1.fastq)
    echo "Processing sample: ${base}"
    
    # BWA alignment with default settings
    bwa mem -t ${THREADS} ${REFERENCE_GENOME} \
        raw_data/${base}_1.fastq raw_data/${base}_2.fastq | \
        samtools sort -@ ${THREADS} -o aligned/${base}.sorted.bam
    
    # Index BAM files
    samtools index aligned/${base}.sorted.bam
    
    # Add read groups (required for GATK)
    gatk AddOrReplaceReadGroups \
        -I aligned/${base}.sorted.bam \
        -O aligned/${base}.rg.bam \
        -RGID ${base} \
        -RGLB lib1 \
        -RGPL illumina \
        -RGPU unit1 \
        -RGSM ${base}
    
    samtools index aligned/${base}.rg.bam
done

# VARIANT CALLING AND FILTERING
# Call variants using GATK HaplotypeCaller
for sample in aligned/*.rg.bam; do
    base=$(basename ${sample} .rg.bam)
    echo "Calling variants for: ${base}"
    
    gatk HaplotypeCaller \
        -R ${REFERENCE_GENOME} \
        -I ${sample} \
        -O variants/${base}.g.vcf.gz \
        -ERC GVCF
done

# Joint genotyping
gatk GenomicsDBImport \
    -V variants/*.g.vcf.gz \
    --genomicsdb-workspace-path variants/genomicsdb \
    --intervals ${REFERENCE_GENOME%.*}.intervals

# Joint calling
gatk GenotypeGVCFs \
    -R ${REFERENCE_GENOME} \
    -V gendb://variants/genomicsdb \
    -O variants/joint_calls.vcf.gz

# Filter variants
echo "Filtering variants..."
gatk VariantFiltration \
    -R ${REFERENCE_GENOME} \
    -V variants/joint_calls.vcf.gz \
    -O variants/filtered_variants.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

# Select only PASS variants
gatk SelectVariants \
    -R ${REFERENCE_GENOME} \
    -V variants/filtered_variants.vcf.gz \
    -O variants/final_snps.vcf.gz \
    --exclude-filtered

# POPULATION STRUCTURE ANALYSIS

# Convert VCF to PLINK format for ADMIXTURE
vcftools --gzvcf variants/final_snps.vcf.gz \
    --plink --out results/admixture/snps

# Convert to BED format
plink --file results/admixture/snps \
    --make-bed --out results/admixture/snps

# Run ADMIXTURE for K=1 to K=10
echo "Running ADMIXTURE analysis..."
cd results/admixture/
for K in {1..10}; do
    echo "Running ADMIXTURE for K=${K}"
    admixture --cv snps.bed ${K} > admixture_K${K}.log 2>&1
done

# Find optimal K (lowest CV error)
grep -h CV *.log > cv_errors.txt
cd ../../

# PHYLOGENETIC ANALYSIS
echo "Building phylogenetic tree..."
cd results/phylogeny/
treebest nj -b 1000 alignment.phy > tree_bootstrap.nwk
cd ../../

# PRINCIPAL COMPONENT ANALYSIS
# Convert VCF to GCTA format
vcftools --gzvcf variants/final_snps.vcf.gz \
    --plink --out results/pca/snps

plink --file results/pca/snps \
    --make-bed --out results/pca/snps

# Run PCA using GCTA
gcta64 --bfile results/pca/snps \
    --make-grm --out results/pca/grm

gcta64 --grm results/pca/grm \
    --pca --out results/pca/pca_results

#LINKAGE DISEQUILIBRIUM DECAY ANALYSIS
# Create population files for each population
# ErZJ (n=22), ErJH (n=12), ErZS (n=16), ErSA (n=22), EtZJ (n=20), EtPA (n=25)
# Run PopLDdecay for each population
for pop in ErZJ ErJH ErZS ErSA EtZJ EtPA; do
    echo "Calculating LD decay for population: ${pop}"
    
    # Create population-specific VCF
    vcftools --gzvcf variants/final_snps.vcf.gz \
        --keep ${pop}_samples.txt \
        --recode --recode-INFO-all \
        --out results/ld_decay/${pop}
    
    # Run PopLDdecay with default parameters
    PopLDdecay -InVCF results/ld_decay/${pop}.recode.vcf \
        -OutDir results/ld_decay/ \
        -OutPrefix ${pop}
done

#DEMOGRAPHIC HISTORY ANALYSIS
# PSMC analysis using heterozygous sites for effective population size estimation
echo "Running PSMC analysis..."
for sample in aligned/*.rg.bam; do
    base=$(basename ${sample} .rg.bam)
    echo "Running PSMC for sample: ${base}"
    
    # Generate consensus sequence focusing on heterozygous sites
    samtools mpileup -C50 -uf ${REFERENCE_GENOME} ${sample} | \
        bcftools call -c - | \
        vcfutils.pl vcf2fq -d 10 -D 100 | \
        gzip > results/psmc/${base}.fq.gz
    
    # Convert to PSMC input format
    fq2psmcfa -q20 results/psmc/${base}.fq.gz > results/psmc/${base}.psmcfa
    
    # Run PSMC with specified parameters: -N 30 -t 15 -r 5 -p 4+25*2+4+6
    psmc -N30 -t15 -r5 -p "4+25*2+4+6" \
        -o results/psmc/${base}.psmc results/psmc/${base}.psmcfa
    
    # Plot individual PSMC results
    psmc_plot.pl -u ${MUTATION_RATE} -g ${GENERATION_TIME} \
        results/psmc/${base} results/psmc/${base}.psmc
done

# SMC++ multigenome analysis for demographic history inference
echo "Running SMC++ multigenome analysis..."

# Convert VCF to SMC++ format for each population
for pop in ErZJ ErJH ErZS ErSA EtZJ EtPA; do
    echo "Converting VCF to SMC++ format for population: ${pop}"
    
    # Get sample list for this population
    samples=$(cat ${pop}_samples.txt | tr '\n' ',' | sed 's/,$//')
    
    # Convert to SMC++ format with enhanced handling of varying coverages
    smc++ vcf2smc variants/final_snps.vcf.gz \
        results/psmc/smc_${pop}.smc.gz \
        $(head -1 ${pop}_samples.txt) \
        -d ${samples}
done

# Run SMC++ estimate method with composite likelihood for demographic inference
for pop in ErZJ ErJH ErZS ErSA EtZJ EtPA; do
    echo "Running SMC++ estimate for population: ${pop}"
    
    smc++ estimate \
        --em-iterations 20 \
        --knots 32 \
        --timepoints 32 1000 \
        ${MUTATION_RATE} \
        results/psmc/smc_${pop}.smc.gz \
        -o results/psmc/smc_${pop}_output/
done

# SMC++ split time analysis between E. rhadinum and E. tetradactylum populations
echo "Calculating split times using SMC++ split function..."

# Define population pairs for split time analysis
declare -A population_pairs=(
    ["ErZJ_vs_EtZJ"]="ErZJ EtZJ"
    ["ErJH_vs_EtZJ"]="ErJH EtZJ" 
    ["ErZS_vs_EtZJ"]="ErZS EtZJ"
    ["ErSA_vs_EtPA"]="ErSA EtPA"
)

for pair in "${!population_pairs[@]}"; do
    read -r pop1 pop2 <<< "${population_pairs[$pair]}"
    echo "Calculating split time for: ${pair}"
    
    # Prepare joint SMC format for split analysis
    smc++ vcf2smc variants/final_snps.vcf.gz \
        results/psmc/split_${pair}.smc.gz \
        $(head -1 ${pop1}_samples.txt) \
        ${pop1}:$(cat ${pop1}_samples.txt | tr '\n' ',' | sed 's/,$//') \
        ${pop2}:$(cat ${pop2}_samples.txt | tr '\n' ',' | sed 's/,$//')
    
    # Run split analysis
    smc++ split \
        --em-iterations 20 \
        --timepoints 32 1000 \
        ${MUTATION_RATE} \
        results/psmc/smc_${pop1}_output/model.final.json \
        results/psmc/smc_${pop2}_output/model.final.json \
        results/psmc/split_${pair}.smc.gz \
        -o results/psmc/split_${pair}_output/
done

# Scale results with generation time and create visualizations

# Plot SMC++ results with spline cubic visualization
for pop in ErZJ ErJH ErZS ErSA EtZJ EtPA; do
    smc++ plot \
        --spline cubic \
        --generation-time ${GENERATION_TIME} \
        results/psmc/plot_${pop}.pdf \
        results/psmc/smc_${pop}_output/model.final.json
done

# Plot split time results
for pair in "${!population_pairs[@]}"; do
    echo "Plotting split time results for: ${pair}"
    
    smc++ plot \
        --spline cubic \
        --generation-time ${GENERATION_TIME} \
        results/psmc/split_plot_${pair}.pdf \
        results/psmc/split_${pair}_output/model.final.json
done

# Combine all population plots for comparison
echo "Creating combined demographic history plot..."
all_models=""
for pop in ErZJ ErJH ErZS ErSA EtZJ EtPA; do
    all_models="${all_models} results/psmc/smc_${pop}_output/model.final.json"
done

smc++ plot \
    --spline cubic \
    --generation-time ${GENERATION_TIME} \
    results/psmc/combined_demographics.pdf \
    ${all_models}

# GENOMIC DIFFERENTIATION ANALYSIS
# Calculate FST using sliding windows
echo "Calculating FST values..."
vcftools --gzvcf variants/final_snps.vcf.gz \
    --weir-fst-pop ErZJ_samples.txt \
    --weir-fst-pop EtZJ_samples.txt \
    --fst-window-size 100000 \
    --fst-window-step 10000 \
    --out results/fst/ErZJ_vs_EtZJ

# Identify highly differentiated regions (HDRs)
python3 << 'EOF'
import pandas as pd
import numpy as np
from scipy import stats
fst_data = pd.read_csv('results/fst/ErZJ_vs_EtZJ.windowed.weir.fst', sep='\t')
fst_data = fst_data.dropna(subset=['WEIGHTED_FST'])
fst_data['FST_Z'] = stats.zscore(fst_data['WEIGHTED_FST'])
threshold_99 = np.percentile(fst_data['FST_Z'], 99)
hdrs = fst_data[fst_data['FST_Z'] >= threshold_99]
from statsmodels.stats.multitest import multipletests
_, fdr_pvals, _, _ = multipletests(
    1 - stats.norm.cdf(hdrs['FST_Z']), 
    alpha=0.05, 
    method='fdr_bh'
)

hdrs_filtered = hdrs[fdr_pvals < 0.05]
hdrs_filtered.to_csv('results/fst/HDRs.csv', index=False)

print(f"Identified {len(hdrs_filtered)} highly differentiated regions")
EOF

# Calculate additional population genetic statistics for HDRs
# Nucleotide diversity (π) and Tajima's D
vcftools --gzvcf variants/final_snps.vcf.gz \
    --window-pi 10000 \
    --out results/fst/nucleotide_diversity
vcftools --gzvcf variants/final_snps.vcf.gz \
    --TajimaD 10000 \
    --out results/fst/tajimas_d

# Calculate DXY using pixy
pixy --stats pi,dxy,fst \
    --vcf variants/final_snps.vcf.gz \
    --populations population_file.txt \
    --window_size 10000 \
    --output_folder results/fst/pixy_output/

# Calculate LD in windows
vcftools --gzvcf variants/final_snps.vcf.gz \
    --geno-r2 \
    --ld-window-bp 10000 \
    --min-r2 0.1 \
    --out results/fst/ld_analysis

# Estimate recombination rate using LDhat
echo "Estimating recombination rates..."
# Convert VCF to LDhat format
python3 convert_vcf_to_ldhat.py variants/final_snps.vcf.gz results/fst/ldhat_input

# Run LDhat
interval -seq results/fst/ldhat_input.seq \
    -loc results/fst/ldhat_input.loc \
    -lk lk_files/lk_n100_t0.01 \
    -its 1000000 -samp 2000 -pen 5 \
    -prefix results/fst/ldhat_output

# INTROGRESSION ANALYSIS

# Random sampling of 10 individuals per population

# Create sample lists (10 individuals each)
for pop in EtZJ EtPA ErZJ; do
    head -10 ${pop}_samples.txt > results/introgression/${pop}_subset.txt
done

# Extract subset VCF
vcftools --gzvcf variants/final_snps.vcf.gz \
    --keep results/introgression/EtZJ_subset.txt \
    --keep results/introgression/EtPA_subset.txt \
    --keep results/introgression/ErZJ_subset.txt \
    --keep outgroup_sample.txt \
    --recode --recode-INFO-all \
    --out results/introgression/introgression_subset

# Run D-statistic and f-statistic analysis
python3 << 'EOF'
import subprocess
import os

# Convert VCF to required format for genomics_general
subprocess.run([
    "python", "/path/to/genomics_general/VCF_processing/parseVCF.py",
    "-i", "results/introgression/introgression_subset.recode.vcf",
    "-o", "results/introgression/introgression.geno.gz"
])

# Run ABBABABAwindows.py for D-statistic and fd calculation
# P1 = EtZJ (allopatric E. tetradactylum)
# P2 = EtPA (sympatric E. tetradactylum) 
# P3 = ErZJ (sympatric E. rhadinum)
# O = outgroup

with open('results/introgression/populations.txt', 'w') as f:
    # Write population assignments
    for i in range(10):
        f.write(f"EtZJ_sample_{i}\tP1\n")
        f.write(f"EtPA_sample_{i}\tP2\n")
        f.write(f"ErZJ_sample_{i}\tP3\n")
    f.write("outgroup_sample\tO\n")

# Run ABBA-BABA analysis
subprocess.run([
    "python", "/path/to/genomics_general/ABBABABAwindows.py",
    "-g", "results/introgression/introgression.geno.gz",
    "-o", "results/introgression/ABBABABA_results.csv",
    "-f", "phased",
    "-w", "10000",
    "-m", "20",
    "-s", "10000",
    "-P1", "P1", "-P2", "P2", "-P3", "P3", "-O", "O",
    "--popsFile", "results/introgression/populations.txt"
])
EOF

# Filter significant introgression regions
python3 << 'EOF'
import pandas as pd
import numpy as np

# Read ABBA-BABA results
results = pd.read_csv('results/introgression/ABBABABA_results.csv')

# Filter regions with significant introgression (P < 0.05)
significant = results[results['fd_pvalue'] < 0.05]
significant.to_csv('results/introgression/significant_introgression.csv', index=False)

print(f"Found {len(significant)} regions with significant introgression")
print(f"Fraction of genome showing introgression: {len(significant)/len(results):.4f}")
EOF

# Extract genes from introgressed regions for functional analysis
echo "Extracting genes from introgressed regions..."
bedtools intersect -a genome_annotation.gff3 \
    -b results/introgression/significant_introgression.bed \
    -wa > results/introgression/introgressed_genes.gff3
