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

## 3. Index the reference

```bash
cd /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Dama_gazelle_sequences/dama_fastq/3_Indexing/

module load samtools-1.22.1
module load bwa-mem
module load bwa-0.7.17

REFERENCE=/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Dama_gazelle_sequences/dama_fastq/3_Indexing/Dama_gazelle_hifiasm-ULONT_primary.fasta.gz

samtools faidx $REFERENCE
bwa index $REFERENCE
```

## 4. Align FILES with BWA MEM
```bash
#!/bin/bash

# Base directories
INPUT_DIR=/scratch/bistbs_new/2_Adapter_trimming
OUTPUT_DIR=/scratch/bistbs_new/4_aligning_with_BWA_Mem_Final_1
REFERENCE=/scratch/bistbs_new/3_Indexing/Dama_gazelle_primary.fasta

# Make output directory
mkdir -p $OUTPUT_DIR

# Loop through sample folders
for SAMPLE_DIR in $INPUT_DIR/SRR*; do
    SAMPLE=$(basename $SAMPLE_DIR)
    
    # Forward and reverse reads
    R1="$SAMPLE_DIR/${SAMPLE}_1_val_1.fq"
    R2="$SAMPLE_DIR/${SAMPLE}_2_val_2.fq"

    # Decide number of threads (24 for SRR17129394, 12 for others)
    if [ "$SAMPLE" == "SRR17129394" ]; then
        THREADS=20
        MEM=128G
        TIME="300:00:00"
    else
        THREADS=8
        MEM=100G
        TIME="200:00:00"
    fi

    # Submit SLURM job
    sbatch --job-name=BWA-$SAMPLE \
           --output=$OUTPUT_DIR/${SAMPLE}_bwa.out \
           --error=$OUTPUT_DIR/${SAMPLE}_bwa.err \
           --ntasks=1 \
           --cpus-per-task=$THREADS \
           --mem=$MEM \
           --time=$TIME \
           --partition=batch \
           --mail-type=BEGIN,END \
           --mail-user=bistbs@miamioh.edu \
           <<EOF
#!/bin/bash
module load samtools-1.22.1
module load bwa-0.7.17
cd $OUTPUT_DIR

echo "Processing sample: $SAMPLE"
echo "Read1: $R1"
echo "Read2: $R2"

bwa mem -t $THREADS $REFERENCE $R1 $R2 \
    | samtools view -bF 4 - > $OUTPUT_DIR/${SAMPLE}_mapped.bam
EOF

done
```

## 5. Sort the SAM
```bash
#!/bin/bash -l
# To be submitted by: sbatch slurm_job.sh
#SBATCH --time=300:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --job-name=Sort_SAM


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load Dependencies and setup env variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load Java
module load java-20 

# Path to Picard
PICARD_JAR="/scratch/bistbs_new/4_aligning_with_BWA_Mem_Final_1/picard.jar"

# Input directory (where your BAM files are)
INPUT_DIR="/scratch/bistbs_new/4_aligning_with_BWA_Mem_Final_1"
OUTPUT_DIR="${INPUT_DIR}/Sorted_BAMs"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Loop through all BAM files and sort them
for bam in ${INPUT_DIR}/*.bam; do
    base=$(basename "${bam}" .bam)
    echo "Sorting ${bam} ..."
    
    java -Xmx12g -jar "${PICARD_JAR}" SortSam \
        I="${bam}" \
        O="${OUTPUT_DIR}/${base}_sorted.bam" \
        SORT_ORDER=coordinate \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=SILENT
done

echo "✅ All BAM files have been sorted and indexed successfully."
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
