#### I have 5 individuals. 
- Addra: SRR17129394, SRR17134085, SRR17134086
- Mhorr: SRR17134087, SRR17134088
  ---
SRR17129394 SRR17134087
SRR17129394 SRR17134088
SRR17134085 SRR17134087
SRR17134085 SRR17134088
SRR17134086 SRR17134087
SRR17134086 SRR17134088
---
#### Rung ANGSD for six pairs and generate six sfs
#### Step 1. and Step 2.
```bash
#!/bin/bash -l
#SBATCH --job-name=Dama_6Pairs_Fixed
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=13
#SBATCH --mem=32G
#SBATCH --partition=batch
#SBATCH --array=1-6
#SBATCH --output=logs/dama_pair_%A_%a.log

# ----------------------------
# Paths & References
# ----------------------------
OUT_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD"
REF="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Dama_gazelle_hifiasm-ULONT_primary.fasta"

cd $OUT_DIR
mkdir -p Individual_SAFs Pairwise_SFS logs

# ----------------------------
# Get Full Paths for This Task
# ----------------------------
BAM_A=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pairs.txt | awk '{print $1}')
BAM_B=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pairs.txt | awk '{print $2}')

# Extract Sample IDs for naming (e.g., SRR17129394)
ID_A=$(basename $BAM_A _10x.bam)
ID_B=$(basename $BAM_B _10x.bam)

echo "------------------------------------------------------"
echo "Job Array ID: ${SLURM_ARRAY_TASK_ID}"
echo "Pair: $ID_A vs $ID_B"
echo "BAM A: $BAM_A"
echo "BAM B: $BAM_B"
echo "------------------------------------------------------"

# ----------------------------
# STEP 1: Individual SAF Generation
# ----------------------------
for BAM in $BAM_A $BAM_B; do
    ID=$(basename $BAM _10x.bam)
    if [ ! -f Individual_SAFs/${ID}.saf.idx ]; then
        echo "Generating SAF for $ID..."
        angsd -P 13 -i $BAM -anc $REF -ref $REF -out Individual_SAFs/${ID} \
              -dosaf 1 -gl 1 -C 50 -minQ 20 -minmapq 30
    else
        echo "SAF for $ID already exists. Skipping."
    fi
done

# ----------------------------
# STEP 2: Pairwise 2DSFS Generation
# ----------------------------
echo "Generating Pairwise SFS for $ID_A vs $ID_B..."

realSFS Individual_SAFs/${ID_A}.saf.idx Individual_SAFs/${ID_B}.saf.idx -P 13 \
> Pairwise_SFS/${ID_A}_${ID_B}.real.sfs

echo "Done. Pairwise SFS saved to Pairwise_SFS/${ID_A}_${ID_B}.real.sfs"
```
#### Step 3. 
```bash
#!/bin/bash -l
#SBATCH --job-name=MiSTI_Convert
#SBATCH --time=54:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=40G
#SBATCH --partition=batch
#SBATCH --output=logs/misti_convert_%j.log

# ----------------------------
# Paths
# ----------------------------
BASE_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD"
SFS_DIR="${BASE_DIR}/Pairwise_SFS"
CONVERTER="${BASE_DIR}/MiSTI/utils/ANGSDSFS.py"

cd $SFS_DIR

# ----------------------------
# Conversion Loop
# ----------------------------
# This finds every .real.sfs file we just made and converts it
for f in *.real.sfs; do
    # Extract IDs (e.g., SRR17129394 and SRR17134087)
    # This assumes filenames are ID1_ID2.real.sfs
    ID1=$(echo $f | cut -d'_' -f1)
    ID2=$(echo $f | cut -d'_' -f2 | cut -d'.' -f1)
    
    echo "Converting Pair: $ID1 and $ID2..."
    
    # Run the python converter
    # Argument 1: The raw SFS file
    # Argument 2: Name of Individual 1 (must match PSMC filename)
    # Argument 3: Name of Individual 2 (must match PSMC filename)
    python $CONVERTER $f $ID1 $ID2 > ${ID1}_${ID2}.mi.sfs
done

echo "Done! You now have 6 .mi.sfs files ready for the final split-time search."
```

#### I will generate the 6 pairwise SFS files.


python3 ~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/MiSTI/utils/calc_time.py \
--funits /home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/dama.units.txt \
/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/SRR17129394.psmc \
/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/SRR17134087.psmc \
> /home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/Index_files/Final_94_87RESULT
