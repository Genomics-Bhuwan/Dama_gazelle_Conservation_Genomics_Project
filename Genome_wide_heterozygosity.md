#### Dama gazelle Genome-wide heterozygosity

#!/bin/bash -l
#SBATCH --job-name=ClapperRail_Het
#SBATCH --time=68:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --output=clapper_rail_het_%A.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# ----------------------------
# 1. Define Tool Paths & Files
# ----------------------------
ANGSD="/software/ngsTools/1.0.2/ngsTools/ANGSD/angsd/angsd"
REALSFS="/software/ngsTools/1.0.2/ngsTools/ANGSD/angsd/misc/realSFS"

# Paths updated for Clapper Rail
REFERENCE="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/Clapper_rail/2_Adapter_trimming/GCA_028554615.1_bRalCre1.1_genomic.fna"
BAM_INPUT="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/Clapper_rail/05_Remove_duplicates/SRR23269683_downsampled_28x.bam"
OUT_PREFIX="Clapper_Rail_SRR23269683"

# Setup Output Directory (Following your path structure)
OUTPUT_DIR="/shared/jezkovt_bistbs_shared/Guam_Rail/Heterozygosity/Comparative_genomics/Heterozygosity_all/Clapper_rail"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# ----------------------------
# 2. Step 1: Generate SAF File
# ----------------------------
echo "[$(date)] Starting Step 1: Generating SAF file for Clapper Rail..."

$ANGSD -P 8 \
     -i "$BAM_INPUT" \
     -anc "$REFERENCE" \
     -ref "$REFERENCE" \
     -out "$OUT_PREFIX" \
     -dosaf 1 \
     -gl 1 \
     -C 50 \
     -minQ 20 \
     -minmapq 30 \
     -baq 1 \
     -uniqueOnly 1 \
     -remove_bads 1

# ----------------------------
# 3. Step 2: Calculate Folded SFS
# ----------------------------
echo "[$(date)] Starting Step 2: Running realSFS..."

if [ -f "${OUT_PREFIX}.saf.idx" ]; then
    $REALSFS "${OUT_PREFIX}.saf.idx" \
         -P 8 \
         -fold 1 \
         -maxIter 100 \
         > "${OUT_PREFIX}.est.ml"
else
    echo "ERROR: SAF index file not found for Clapper Rail."
    exit 1
fi

# ----------------------------
# 4. Final Calculation & Output File
# ----------------------------
echo "------------------------------------------------"
# Calculation: Heterozygous sites / (Homozygous + Heterozygous sites)
HET=$(awk '{print $2/($1+$2)}' "${OUT_PREFIX}.est.ml")

# Save result to a final output file as requested
echo -e "Species\tSample\tHeterozygosity" > "${OUT_PREFIX}_final_heterozygosity.txt"
echo -e "Clapper_Rail\tSRR23269683\t$HET" >> "${OUT_PREFIX}_final_heterozygosity.txt"

echo "Clapper Rail Heterozygosity: $HET"
echo "Results saved to: ${OUTPUT_DIR}/${OUT_PREFIX}_final_heterozygosity.txt"
echo "------------------------------------------------"
echo "[$(date)] Pipeline completed."
