## Mark and Remove Duplicates

## This script marks and removes the PCR duplicates from the merged BAM file.
----
```bash
#!/bin/bash -l
#SBATCH --job-name=MarkDuplicates
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=120G
#SBATCH --partition=batch
#SBATCH --nodelist=mualhpcp36.mpi          # force run on node 36
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=logs/MarkDup_%A.out
#SBATCH --error=logs/MarkDup_%A.err

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Load modules and variables        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
module load java-20

PICARD_JAR="/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/picard.jar"
INPUT_BAM="/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/all_samples_merged.bam"
OUTPUT_DIR="/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates"
SCRATCH="/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/tmp_MarkDuplicates"

# Create directories if missing
mkdir -p "$OUTPUT_DIR" "$SCRATCH" logs

# Output files
MARKED_BAM="${OUTPUT_DIR}/all_samples_merged_rmdup.bam"
METRICS_FILE="${OUTPUT_DIR}/all_samples_merged_rmdup.metrics"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Run Picard MarkDuplicates Tool      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo "🔍 Running MarkDuplicates on: $INPUT_BAM ..."
java -Xmx100g -jar "$PICARD_JAR" MarkDuplicates \
    I="$INPUT_BAM" \
    O="$MARKED_BAM" \
    METRICS_FILE="$METRICS_FILE" \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR="$SCRATCH" \
    MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 \
    MAX_RECORDS_IN_RAM=50000 \
    CREATE_INDEX=true

echo "✅ Finished MarkDuplicates"
echo "Output written to: $MARKED_BAM"

```
---
