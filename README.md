# Dama gazelle Whole-Genome Resequencing (WGS) Variant Calling Pipeline Tutorial

**Author:** Bhuwan Singh Bist  

**Affiliation:** Jezkova lab  

**Date:** 2025-10-03

This tutorial demonstrates a complete workflow for WGS variant calling in the Dama gazelle. The pipeline includes:

1. Quality control (FastQC)  
2. Adapter trimming (Trim Galore + Cutadapt)  
3. Mapping reads to a haplotype-resolved reference genome  
4. BAM processing (Picard & Samtools)  
5. Variant calling and filtering (GATK & VCFtools)

> **Note:** All code blocks are for **demonstration purposes only** and are not executed in this document.

---

# Dama gazelle WGS Variant Calling Pipeline

This tutorial demonstrates the workflow for whole-genome sequencing (WGS) variant calling in the Dama gazelle. Each step includes the commands in a copyable code block.

---

## 1. Quality control using FastQC

```bash
#!/bin/bash -l
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --job-name=Fastqc

module load fastqc

THREADS=24
INDIR="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Dama_gazelle_sequences/dama_fastq"
OUTDIR="${INDIR}/1_Fastqc"

cd $INDIR
fastqc *_1.fastq *_2.fastq -o $OUTDIR -t $THREADS
```

## 2. Adapter trimming with Trim Galore

```bash
module load trimgalore-0.6.7
# Ensure cutadapt v5.1 is available

THREADS=24
INDIR="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Dama_gazelle_sequences/dama_fastq"
OUTDIR="${INDIR}/2_Adapter_trimming"

cd $INDIR

for R1 in *_1.fastq; do
    SAMPLE=${R1%_1.fastq}
    R2="${SAMPLE}_2.fastq"
    echo "Processing sample: $SAMPLE"
    trim_galore --paired --cores $THREADS --output_dir $OUTDIR $R1 $R2
done
```

## 3a. Index the reference

```bash
cd /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Dama_gazelle_sequences/dama_fastq/3_Indexing/

module load samtools-1.22.1
module load bwa-mem
module load bwa-0.7.17

REFERENCE=/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Dama_gazelle_sequences/dama_fastq/3_Indexing/Dama_gazelle_hifiasm-ULONT_primary.fasta.gz

samtools faidx $REFERENCE
bwa index $REFERENCE
```

## 3b. Align reads and process BAMs

```bash
cd /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Dama_gazelle_sequences/dama_fastq/4_aligning_with_BWA_Mem/

module load samtools-1.22.1
module load bwa-mem
module load bwa-0.7.17

INPUT_DIR=/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Dama_gazelle_sequences/dama_fastq/4_aligning_with_BWA_Mem
OUTPUT_DIR=/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Dama_gazelle_sequences/dama_fastq/4_aligning_with_BWA_Mem/aligning_out
REFERENCE=$INPUT_DIR/Dama_gazelle_hifiasm-ULONT_primary.fasta.gz

mkdir -p $OUTPUT_DIR
THREADS=24

for R1 in $INPUT_DIR/*_1_val_1.fq; do
    BASE=$(basename $R1 _1_val_1.fq)
    R2="${INPUT_DIR}/${BASE}_2_val_2.fq"
    
    echo "Processing sample: $BASE"
    
    bwa mem -t $THREADS $REFERENCE $R1 $R2 \
        | samtools view -bF 4 - > $OUTPUT_DIR/${BASE}_mapped.bam

    samtools fixmate -m $OUTPUT_DIR/${BASE}_mapped.bam $OUTPUT_DIR/${BASE}.fixmate.bam
    samtools sort $OUTPUT_DIR/${BASE}.fixmate.bam -o $OUTPUT_DIR/${BASE}.sorted.bam

    picard AddOrReplaceReadGroups \
        I=$OUTPUT_DIR/${BASE}.sorted.bam \
        O=$OUTPUT_DIR/${BASE}.rg.bam \
        RGID=$BASE RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$BASE

    picard MarkDuplicates \
        I=$OUTPUT_DIR/${BASE}.rg.bam \
        O=$OUTPUT_DIR/${BASE}.dedup.bam \
        M=$OUTPUT_DIR/${BASE}.dedup.metrics

    samtools index $OUTPUT_DIR/${BASE}.dedup.bam

    gatk --java-options "-Xmx4g" RealignerTargetCreator \
        -R $REFERENCE -I $OUTPUT_DIR/${BASE}.dedup.bam \
        -O $OUTPUT_DIR/${BASE}.intervals

    gatk --java-options "-Xmx4g" IndelRealigner \
        -R $REFERENCE -I $OUTPUT_DIR/${BASE}.dedup.bam \
        -targetIntervals $OUTPUT_DIR/${BASE}.intervals \
        -O $OUTPUT_DIR/${BASE}.realigned.bam
done
```

## 4. Variant calling (multi-sample)

```bash
BAM_LIST=$(ls $OUTPUT_DIR/*.realigned.bam | tr '\n' ' ')

gatk --java-options "-Xmx8g" HaplotypeCaller \
    -R $REFERENCE \
    -I $BAM_LIST \
    -O $OUTPUT_DIR/dama_raw.vcf
```

## 5. Variant filtering

```bash
vcftools --vcf $OUTPUT_DIR/dama_raw.vcf \
         --remove-indels \
         --minQ 30 \
         --recode --recode-INFO-all \
         --out $OUTPUT_DIR/dama_filtered

echo "Pipeline complete. Filtered VCF is at $OUTPUT_DIR/dama_filtered.recode.vcf"
```
