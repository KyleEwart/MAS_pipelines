#!/usr/bin/env bash

# Variables
SAMPLE=TG_003_08a_50plex_v2
fastq=../data/5.demux/barcode05/fastq_runid_unknown_0.fastq
REFERENCE=Tiger-SNPs_sequences-templates_biallelic_longer.fas
SNP_POSITIONS=snp-positions.txt
MIN_COVERAGE=8
THREADS=2
REMOVE_FILES=true

# snippy binaries path. Hash out if not required.
export PATH=/Users/kyleewart/anaconda3/envs/snippy_env/bin:$PATH
# Activate conda environment:
eval "$(conda shell.bash hook)"
conda activate snippy_env

# Run the alignment using minimap2
#minimap2 -ax map-ont -t ${THREADS} -I 3G ${REFERENCE} ${fastq} > ${SAMPLE}.sam
bwa index ${REFERENCE}
bwa mem -x ont2d -t ${THREADS} ${REFERENCE} ${fastq} > ${SAMPLE}.sam
# Convert and sort
samtools view -bT ${REFERENCE} ${SAMPLE}.sam | samtools sort -l 0 -T temp - > ${SAMPLE}.bam

# Index the bam (Not sure if I need this)
samtools index ${SAMPLE}.bam

# Getting alignment stats
samtools flagstat ${SAMPLE}.bam > ${SAMPLE}_alignment-stats.txt

# I think this allows freebayes to run in parallel. 100 might need to be altered, and this step might not be needed.
# More info on this using 'freebayes-parallel -h'
fasta_generate_regions.py ${REFERENCE}.fai 1000 > ${REFERENCE}.txt

# Call SNPs with freebayes
freebayes-parallel ${REFERENCE}.txt ${THREADS} --min-coverage ${MIN_COVERAGE} --min-alternate-fraction 0.1 --pooled-continuous --strict-vcf -f ${REFERENCE} ${SAMPLE}.bam > ${SAMPLE}.raw.vcf

# Filtering the snps
# If this causes problems, change '(FMT/AO)/(FMT/DP)' back to '>=0' 
bcftools view --include "QUAL>=50 && FMT/DP>=$MIN_COVERAGE && (FMT/AO)/(FMT/DP)>=0.2" ${SAMPLE}.raw.vcf | vt normalize -r ${REFERENCE} - | bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > ${SAMPLE}.filt.vcf
# QUAL>=100: This filter specifies that variants with a quality score (QUAL) of 100 or higher should be included in the output. The QUAL field represents the phred-scaled quality score for the variant, indicating the confidence of the variant call.
# FMT/DP>=10: This filter ensures that variants have a depth (DP) value of 10 or greater in the sample(s). The FMT/DP field represents the read depth at the variant position, indicating the number of reads covering that position.
# (FMT/AO)/(FMT/DP)>=0: This filter calculates the allele frequency of the variant (AO) and divides it by the read depth (DP) to obtain the allele frequency ratio. The filter ensures that the allele frequency ratio is greater than or equal to zero, meaning there is at least one read supporting the variant allele. The FMT/AO field represents the alternate allele observation count.

# What snippy uses:
#bcftools view --include 'FMT/GT="1/1" && QUAL>=100 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0' MAS_alignment.raw.vcf | vt normalize -r Tiger-SNPs_sequences-templates_biallelic_longer.fas - | bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > MAS_alignment.filt.vcf

# Format to tab file:
snippy-vcf_to_tab --ref ${REFERENCE} --vcf ${SAMPLE}.filt.vcf > ${SAMPLE}.filt.tab
# /home/ubuntu/miniconda3/envs/snippy_env/bin/snippy-vcf_to_tab --gff reference/ref.gff --ref reference/ref.fa --vcf MAS_alignment.filt.vcf > MAS_alignment.filt.tab

# If the SNP positions file is provided, then this filters those specific SNPs.
# This if statement tests if the SNP_POSITIONS variable contains an empty string
# These lines also remove columns that aren't needed
if [ -z ${SNP_POSITIONS} ]; then
    tail -n +2 ${SAMPLE}.filt.tab | awk -F "\t" -v OFS='\t' '{print $1,$2,$6}' > ${SAMPLE}.filt-target.tab
else
    grep -Ff ${SNP_POSITIONS} ${SAMPLE}.filt.tab | awk -F "\t" -v OFS='\t' '{print $1,$2,$6}' > ${SAMPLE}.filt-target.tab
fi

# Add column headings to the file
echo -e "SNP\tPosition\tEvidence" | cat - ${SAMPLE}.filt-target.tab > ${SAMPLE}.snps.tsv

# Check the condition and delete intermediate files if true
if [ ${REMOVE_FILES} = true ]; then
    echo "Deleting intermediate files..."
    rm ${SAMPLE}.sam ${SAMPLE}.bam* ${SAMPLE}_alignment-stats.txt ${SAMPLE}.*.vcf ${SAMPLE}.*.tab ${REFERENCE}.*
fi
