## Find out the coverage of your five BAM files.


----
```bash
#!/bin/bash -l
#SBATCH --job-name=BAM_Coverage
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=50G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=logs/Coverage_%A.out
#SBATCH --error=logs/Coverage_%A.err

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Load dependencies and set paths   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
module load roslin/samtools/1.10

# Path to your merged BAM
BAM="/scratch/bistbs_new/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/8_MarkDuplicates/all_samples_merged_rmdup.bam"

# Output coverage file
OUTPUT="${BAM%.bam}_coverage.gz"

# Create logs folder if it does not exist
mkdir -p logs

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          Run coverage calculation         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo "Calculating coverage for $BAM ..."
samtools coverage "$BAM" | gzip > "$OUTPUT"

echo "✅ Coverage stats saved to $OUTPUT"

```
---