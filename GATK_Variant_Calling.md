# The Genome Analysis Toolkit (GATK)
---
The GATK (Genome Analysis Toolkit) is one of the most used programs for genotype calling in sequencing data in model and non model organisms. 
Designed to analyze human genetic data and all its pipelines are optimized for the current purpose.
---

**Author:** Bhuwan Singh Bist

**Date:** 2025-10-16

#### Structural variants(SVs) are DNA rearrangements that involve at least 50 nucleotides. 
#### By genetic standards, these mutations are fairly large, and they are significantly abundant.
#### SVs are one of the strongest forces that direct genome evolution and diseases. 
#### They are impactful coz they encompass a number of mutational classes that can disrupt proteion-coding genes and cis-regulatory architecture through different mechanisms. 
#### Use GATK-SV if you want these Structural Variants such as:
a. Copy number variants (CNVs), including deletions and duplications
b. Insertions
c. Inversions
d. Reciprocal chromosomal translocations
e. Additional forms of complex structural variation
#### Overall, it can identify, genotype, and annotate structural variation from the above types of variants:
#### Below is the image showing how different types of strcutural variants look like and what this pipeline can detect.
---
This tutorial demonstrates a complete workflow for variant calling using GATK. The pipeline includes:
<img width="2639" height="2056" alt="image" src="https://github.com/user-attachments/assets/9ddc14d2-a8bd-4131-a78f-af4f038267f2" />
---


###### Step 1. Haplotype Caller
###### Aligned BAM file(reads mapped to reference genome) and call all possible variants(SNPs and indels) from the samples.
###### Generates Genomic VCF(GVCF) which contains both variants and non-variants.

```bash
#!/bin/bash -l
#SBATCH --job-name=HaplotypeCaller
#SBATCH --time=300:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=logs/HaplotypeCaller_%A.out
#SBATCH --error=logs/HaplotypeCaller_%A.err

# Load GATK
module load gatk-4.1.2.0 

# Input BAM and reference
BAM=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates/all_samples_merged_rmdup.bam
REFERENCE=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates/Dama_gazelle_hifiasm-ULONT_primary.fasta

# Output folder
GVCF_DIR=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates/GVCFs
mkdir -p $GVCF_DIR

# Run HaplotypeCaller
gatk HaplotypeCaller \
    -R $REFERENCE \
    -I $BAM \
    -O $GVCF_DIR/all_samples.g.vcf.gz \
    -ERC GVCF \
  
echo "✅ HaplotypeCaller finished. Output: $GVCF_DIR/all_samples.g.vcf.gz"
```





