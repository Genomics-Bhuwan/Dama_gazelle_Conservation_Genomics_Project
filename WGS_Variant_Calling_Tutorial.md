# Dama Gazelle Whole-Genome Resequencing (WGS) Variant Calling Pipeline Tutorial

**Author:** Bhuwan Singh Bist

**Affiliation:** Jezkova lab

**Date:** 2025-10-03

This tutorial demonstrates a complete workflow for WGS variant calling in the Dama gazelle. The pipeline includes:

1. Quality control (FastQC)
2. Adapter trimming (Trim Galore + Cutadapt)
3. Mapping reads to a haplotype-resolved Ruminant Telomere-to-Telomere(T2T)reference genome assembly of Dama gazelle(Nanger dama)
4. BAM processing (Picard & Samtools)
5. Variant calling and filtering (GATK & VCFtools)

> **Note:** All code blocks are for **demonstration purposes only** and are not executed in this document.

---

# Dama Gazelle WGS Variant Calling Pipeline

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

REFERENCE=/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Dama_gazelle_sequences/dama_fastq/3_Indexing/Dama_gazelle_hifiasm-ULONT_primary.fasta

samtools faidx $REFERENCE
bwa index $REFERENCE
```

## 3b. Align reads and process BAMs

```bash
# Change directory to where the T2T reference is stored.
cd /scratch/bistbs/   

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load required modules                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
module load samtools-1.22.1
module load bwa-mem
module load bwa-0.7.17 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# T2T consertium generated reference assembly#
# Indexing with samtools and BWA             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

REFERENCE=/scratch/bistbs/Dama_gazelle_hifiasm-ULONT_primary.fasta

# Index with samtools (for downstream use)
samtools faidx $REFERENCE

# Index with BWA (for mapping)
bwa index $REFERENCE
```
## 3C. Align reads and process BAMs using BWA-MEM

```bash
# Change directory to where the T2T reference is stored.
cd /scratch/bistbs/  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load required modules                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
module load samtools-1.22.1
module load bwa-mem
module load bwa-0.7.17 

# Directories
INPUT_DIR=/scratch/bistbs/
OUTPUT_DIR=/scratch/bistbs/4_Aligning_BWA/Alignment_output
REFERENCE=/scratch/bistbs/Dama_gazelle_hifiasm-ULONT_primary.fasta

# Make output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Number of threads (matches ntasks-per-node)
THREADS=24

# Loop through all paired-end files
for R1 in $INPUT_DIR/*_1_val_1.fq; do
    BASE=$(basename $R1 _1_val_1.fq)
    R2="${INPUT_DIR}/${BASE}_2_val_2.fq"
    
    echo "Processing sample: $BASE"
    echo "Read1: $R1"
    echo "Read2: $R2"
    
# Align reads with BWA-MEM and output BAM (mapped only) and Process and filter unmapped reads with samtools
#bF 4-binary format and remove reads that are 4. 4 means reads are unmapped.
    bwa mem -t $THREADS $REFERENCE $R1 $R2 \
        | samtools view -bF 4 - > $OUTPUT_DIR/${BASE}_mapped.bam
done
```    
## 3d. Sort the mapped BAM files
```bash
# Change directory to where the T2T reference is stored.
cd /scratch/bistbs/4_Aligning_BWA/Alignment_output/ 


# Number of threads (matches ntasks-per-node)
THREADS=24

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load Dependencies and Setup Environment   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Load Java (Picard requires Java)
module load java-20

# Directories
INPUT_DIR=/scratch/bistbs/4_Aligning_BWA/Alignment_output
OUTPUT_DIR=/scratch/bistbs/4_Aligning_BWA/Alignment_output/5_Sort_Bam
REFERENCE=/scratch/bistbs/Dama_gazelle_hifiasm-ULONT_primary.fasta
SCRATCH=/scratch/bistbs/tmp

# Create output and temp directories if not exist
## mkdir -p "$OUTPUT_DIR" ##It already exists.
mkdir -p "$SCRATCH"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Merge BAM Files                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

cd "$INPUT_DIR" || { echo "Error: Input directory not found."; exit 1; }

# List all BAM files to be merged
bam_files=($(ls *_mapped.bam))
if [ ${#bam_files[@]} -eq 0 ]; then
    echo "No BAM files found in $INPUT_DIR"
    exit 1
fi

# Prepare Picard input argument list
bamlist=""
for f in "${bam_files[@]}"; do
    bamlist+="I=$INPUT_DIR/$f "
done

# Define output merged BAM name
OUTPUT_BAM=$OUTPUT_DIR/All_samples_merged.bam

# Run Picard MergeSamFiles
java -Xmx60g -jar /scratch/bistbs/4_Aligning_BWA/Alignment_output/picard.jar MergeSamFiles \
    $bamlist \
    O="$OUTPUT_BAM" \
    USE_THREADING=true \
    TMP_DIR="$SCRATCH" \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT

echo "Merge complete. Output written to: $OUTPUT_BAM"

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
