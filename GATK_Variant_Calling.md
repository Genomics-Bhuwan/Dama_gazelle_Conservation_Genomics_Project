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
###### Since, my deduplicated bam file consists of five samples. I want to know the name of the samples.

```bash
samtools view -H all_samples_merged_rmdup.bam | grep '^@RG' | awk '{for(i=1;i<=NF;i++) if($i ~ /^SM:/) print $i}' | sed 's/SM://'
```

### Run the batch script 
### Run all samples in parallel (fast and efficient)
```bash
#!/bin/bash -l
#SBATCH --job-name=HaplotypeCaller_parallel
#SBATCH --time=300:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --partition=bigmem
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=logs/HaplotypeCaller_%A.out
#SBATCH --error=logs/HaplotypeCaller_%A.err

# -------------------------------
# Load modules
# -------------------------------
module load gatk-4.1.2.0
module load samtools-1.22.1

# -------------------------------
# Input/output paths
# -------------------------------
BAM=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates/all_samples_merged_rmdup.bam
REFERENCE=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates/Dama_gazelle_hifiasm-ULONT_primary.fasta
GVCF_DIR=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates/GVCFs

mkdir -p $GVCF_DIR

# -------------------------------
# Extract sample names
# -------------------------------
samples=$(samtools view -H $BAM | grep '^@RG' | awk '{for(i=1;i<=NF;i++) if($i ~ /^SM:/) print $i}' | sed 's/SM://')

# -------------------------------
# Run HaplotypeCaller for each sample in parallel
# -------------------------------
for sample in $samples; do
  echo "🚀 Launching HaplotypeCaller for: $sample"

  gatk --java-options "-Xmx16g" HaplotypeCaller \
      -R $REFERENCE \
      -I $BAM \
      -O $GVCF_DIR/${sample}.g.vcf.gz \
      -ERC GVCF \
      --sample-name $sample \
      --native-pair-hmm-threads 4 \
      > $GVCF_DIR/${sample}.log 2>&1 &

done

# Wait for all background processes to complete
wait

echo "✅ All parallel HaplotypeCaller jobs finished successfully."
```

#### I should do it Chromosome by Chromosome for structural variants and genomic islands.
#### List of chromosomes (contigs) with a file named chromosomes_list.txt:

```bash
grep "^>" Dama_gazelle_hifiasm-ULONT_primary.fasta | sed 's/>//' > chromosomes_list.txt

```
#### Create a text file named sample_map.txt mapping sample names to their gVCF paths.
```bash
Sample1   /path/to/GVCFs/Sample1.g.vcf.gz
Sample2   /path/to/GVCFs/Sample2.g.vcf.gz
Sample3   /path/to/GVCFs/Sample3.g.vcf.gz
Sample4   /path/to/GVCFs/Sample4.g.vcf.gz
Sample5   /path/to/GVCFs/Sample5.g.vcf.gz
```
#### Get Genomics Database by Chromosomes.
```bash
#!/bin/bash -l
#SBATCH --job-name=GenomicsDBImport
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=1-68%10    # adjust to the number of chromosomes and concurrent tasks
#SBATCH --partition=bigmem
#SBATCH --output=logs/GenomicsDBImport_%A_%a.out
#SBATCH --error=logs/GenomicsDBImport_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# Load GATK
module load gatk-4.1.2.0

# Directories
REF=/localscratch/bistbs/.../Dama_gazelle_hifiasm-ULONT_primary.fasta
SAMPLEMAP=/localscratch/bistbs/.../sample_map.txt
CHROM_LIST=/localscratch/bistbs/.../chromosomes_list.txt
OUTDIR=/localscratch/bistbs/.../GenomeImport
TMPDIR=/localscratch/bistbs/.../GenomeImport_tmp

mkdir -p $OUTDIR $TMPDIR logs

# Get chromosome for this array task
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $CHROM_LIST)

echo "Processing chromosome: $CHR"

gatk --java-options "-Xmx16g" GenomicsDBImport \
  --genomicsdb-workspace-path ${OUTDIR}/${CHR}_gvcf_db \
  --sample-name-map ${SAMPLEMAP} \
  --batch-size 50 \
  -L ${CHR} \
  --tmp-dir ${TMPDIR}

echo "✅ Finished chromosome ${CHR}"
```

#### Step 2. Joint Genotyping
```bash
#!/bin/bash -l
#SBATCH --job-name=GenotypeGVCFs
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=1-68%10
#SBATCH --partition=bigmem
#SBATCH --output=logs/GenotypeGVCFs_%A_%a.out
#SBATCH --error=logs/GenotypeGVCFs_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

module load gatk-4.1.2.0

REF=/localscratch/bistbs/.../Dama_gazelle_hifiasm-ULONT_primary.fasta
CHROM_LIST=/localscratch/bistbs/.../chromosomes_list.txt
DBDIR=/localscratch/bistbs/.../GenomeImport
OUTDIR=/localscratch/bistbs/.../JointGenotyping
TMPDIR=/localscratch/bistbs/.../JointGenotyping_tmp

mkdir -p $OUTDIR $TMPDIR logs

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $CHROM_LIST)

echo "Running GenotypeGVCFs for $CHR"

gatk --java-options "-Xmx16g" GenotypeGVCFs \
  -R ${REF} \
  -V gendb://${DBDIR}/${CHR}_gvcf_db \
  -O ${OUTDIR}/${CHR}_joint.vcf.gz \
  --tmp-dir ${TMPDIR}

echo "✅ Finished chromosome ${CHR}"

```

#### Step 3. Variant Filtration
```bash
# -------------------------------
# Load GATK module
# -------------------------------
module load gatk-4.1.2.0

# -------------------------------
# Input/output paths
# -------------------------------
RAW_VCF=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates/GVCFs/all_samples_joint.vcf.gz
FILTERED_VCF=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates/GVCFs/all_samples_joint_filtered.vcf.gz

# -------------------------------
# Run VariantFiltration
# -------------------------------
gatk VariantFiltration \
   -V $RAW_VCF \
   
   # -------------------------------
   # Filter sites with extremely high depth
   # DP < 1800: remove sites where coverage is unusually high
   # High depth can indicate collapsed repeats or duplicated regions
   # -------------------------------
   -filter "DP < 1800" --filter-name "DP1800" \
   
   # -------------------------------
   # Strand bias filters
   # FS > 60.0: Fisher Strand test; removes sites with strong strand bias
   # SOR > 3.0: Strand Odds Ratio; another measure of strand bias
   # Sites failing these filters are likely sequencing/mapping artifacts
   # -------------------------------
   -filter "FS > 60.0" --filter-name "FS60" \
   -filter "SOR > 3.0" --filter-name "SOR3" \
   
   # -------------------------------
   # Mapping quality filter
   # MQ < 40.0: remove sites where reads are poorly aligned
   # Poorly mapped reads can introduce false variants
   # -------------------------------
   -filter "MQ < 40.0" --filter-name "MQ40" \
   
   # -------------------------------
   # Mapping quality rank sum test
   # MQRankSum < -12.5: tests if alternate alleles have lower mapping quality than reference
   # Extreme negative values indicate alignment bias against alternate alleles
   # -------------------------------
   -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
   
   # -------------------------------
   # Read position rank sum test
   # ReadPosRankSum < -8.0: tests if alternate alleles occur at biased positions in reads
   # Negative values indicate alt alleles mostly at read ends, often errors
   # -------------------------------
   -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
   
   # -------------------------------
   # Output filtered VCF
   # -------------------------------
   -O $FILTERED_VCF
```

#### Step 4.A VCF Filtering for Population Genomics
```bash

module load vcf-tools
RAW_VCF=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates/GVCFs/all_samples_joint_filtered.vcf.gz
VCFTOOLS_OUT=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates/GVCFs/vcftools_filtered

# ------------------------------
# Filter with VCFtools
# ------------------------------
vcftools --gzvcf $RAW_VCF \                # Input VCF (can be gzipped)
         --minQ 30 \                        # Keep only sites with minimum quality score of 30 (high confidence genotypes)
         --remove-indels \                  # Remove insertions/deletions (INDELs), keep only SNPs.
         --recode --recode-INFO-all \       # Produce a new VCF keeping all INFO fields
         --out $VCFTOOLS_OUT                # Output prefix

# ------------------------------
# Check missingness per individual
# ------------------------------
vcftools --vcf ${VCFTOOLS_OUT}.recode.vcf \
         --missing-indv                    # Reports fraction of missing genotypes per sample

# This will produce 'out.miss' file with % missing genotypes per individual
# Individuals with very high missingness can be removed in later filtering steps
```
#### Step 4.B Keep indel if you want to do SNpeff and VEP
```bash
module load vcf-tools
RAW_VCF=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates/GVCFs/all_samples_joint_filtered.vcf.gz
VCFTOOLS_OUT=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates/GVCFs/vcftools_filtered_with_indels

# ------------------------------
# Filter with VCFtools (keep indels)
# ------------------------------
vcftools --gzvcf $RAW_VCF \       # Input VCF (gzipped)
         --minQ 30 \              # Keep sites with minimum quality score of 30
         --recode --recode-INFO-all \  # Produce new VCF with all INFO fields
         --out $VCFTOOLS_OUT      # Output prefix

# ------------------------------
# Check missingness per individual
# ------------------------------
vcftools --vcf ${VCFTOOLS_OUT}.recode.vcf \
         --missing-indv
```


