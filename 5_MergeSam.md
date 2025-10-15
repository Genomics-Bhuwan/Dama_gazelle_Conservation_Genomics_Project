# Merge RG BAMs

This script merges all BAM files with read groups (`_RG.bam`) into a single coordinate-sorted BAM using Picard `MergeSamFiles`.

```bash
#!/bin/bash -l
#SBATCH --job-name=MergeReadGroups
#SBATCH --time=300:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=100G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=logs/MergeRG_%A.out
#SBATCH --error=logs/MergeRG_%A.err
# #SBATCH --array=1-5   # <-- ❌ REMOVE THIS

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#         Load modules and variables        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
module load java-20

PICARD_JAR="/scratch/bistbs_new/4_aligning_with_BWA_Mem_Final_1/picard.jar"
INPUT_DIR="/scratch/bistbs_new/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/"
OUTPUT_DIR="/scratch/bistbs_new/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam"
SCRATCH="/scratch/bistbs_new/tmp_MergeSam"

mkdir -p "$OUTPUT_DIR" "$SCRATCH" logs

# Collect all BAMs
bam_list=(${INPUT_DIR}/*_RG.bam)

# Build input arguments for Picard
bam_args=""
for f in "${bam_list[@]}"; do
    bam_args+="I=$f "
done

MERGED_BAM="${OUTPUT_DIR}/all_samples_merged.bam"

echo "Merging ${#bam_list[@]} BAMs into $MERGED_BAM..."
java -Xmx200g -jar "$PICARD_JAR" MergeSamFiles \
    $bam_args \
    O="$MERGED_BAM" \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    TMP_DIR="$SCRATCH"

echo "✅ Finished merging BAMs"
```
  
