## Index your downsampled BAMS.


----
```bash
#!/bin/bash -l
#SBATCH --job-name=BAM_Index
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=logs/BAM_Index_%A.out
#SBATCH --error=logs/BAM_Index_%A.err

# Load samtools
module load samtools-1.22.1

# Set working directory
WORKDIR=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates
cd $WORKDIR

# Reference genome
REF_GENOME=Dama_gazelle_hifiasm-ULONT_primary.fasta.gz

# BAM file
BAM_FILE=all_samples_merged_rmdup.bam

# Create output directory for indexing
INDEX_DIR=9_Indexing
mkdir -p $INDEX_DIR

# Index the BAM into the new folder
OUTPUT_BAI="$INDEX_DIR/${BAM_FILE}.bai"
echo "Indexing BAM: $BAM_FILE using reference: $REF_GENOME ..."
samtools index -o $OUTPUT_BAI $BAM_FILE

echo "✅ BAM indexing completed. Index saved to $OUTPUT_BAI"

```
---
