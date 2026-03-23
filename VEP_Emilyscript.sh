#### Finding deleterious mutation in the  Dama gazelle.
- Usually kept three species but here I am keeping two outgroups.
- Grant’s gazellehttps://trace.ncbi.nlm.nih.gov/Traces/?run=SRR6878810)  and Thomson's gazelle (https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR6894844)
#### Step 1. Download the short-read sequences of Grant's gazelle and Thomson's gazelle.
- Map these with BWA-MEM using the reference genome assembly of Dama gazelle.
- Since, I had already mapped each outgroup individually to the ref. genome and then already added RG and dedeuplicated, now I am mapping. 
```bash
#!/bin/bash -l
#SBATCH --job-name=MergeOutgroups_Dama
#SBATCH --time=50:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=80G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/merge_%j.out
#SBATCH --error=/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/merge_%j.err

# 1. Load Java (Using your specific version)
module load java-20

# 2. Define Paths from your provided script
PICARD_JAR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/picard.jar"
IN_DIR="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Haplotype_Caller/Five_samples/rmdup_Dama"
OUT_DIR="/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# 3. Define the two specific outgroup files to merge
# Grant's Gazelle and Thompson's Gazelle
BAM1="$IN_DIR/SRR6878810_sorted_RG_rmdup10X.bam"
BAM2="$IN_DIR/SRR6894844_rmdup_10X.bam"

echo "Starting merge for Outgroups..."

# 4. Run Picard MergeSamFiles
java -Xmx30g -jar "$PICARD_JAR" MergeSamFiles \
    I="$BAM1" \
    I="$BAM2" \
    O="$OUT_DIR/gazelle_outgroup_merged.bam" \
    TMP_DIR="$OUT_DIR" \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    USE_THREADING=true

echo "Merge Complete. Output file: $OUT_DIR/gazelle_outgroup_merged.bam"
```
