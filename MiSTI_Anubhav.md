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
#SBATCH --job-name=Dama_6Pairs_Direct
#SBATCH --time=58:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=13
#SBATCH --mem=32G
#SBATCH --partition=batch
#SBATCH --array=0-5
#SBATCH --output=logs/dama_pair_%A_%a.log

# ----------------------------
# Paths & References
# ----------------------------
module load angsd

# Direct path to the specific realSFS version you requested
REALSFS_PATH="/software/ngsTools/1.0.2/ngsTools/ANGSD/angsd/misc/realSFS"

OUT_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD"
REF="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Dama_gazelle_hifiasm-ULONT_primary.fasta"
BASE="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Haplotype_Caller/Five_samples/rmdup_Dama"

cd $OUT_DIR
mkdir -p Individual_SAFs Pairwise_SFS logs

# ----------------------------
# Hardcoded Full Paths for the 6 Pairs
# ----------------------------
PAIRS=(
    "${BASE}/SRR17129394_mapped_sorted_RG_rmdup.bam|${BASE}/SRR17134085_mapped_sorted_RG_rmdup.bam"
    "${BASE}/SRR17129394_mapped_sorted_RG_rmdup.bam|${BASE}/SRR17134086_mapped_sorted_RG_rmdup.bam"
    "${BASE}/SRR17129394_mapped_sorted_RG_rmdup.bam|${BASE}/SRR17134087_mapped_sorted_RG_rmdup.bam"
    "${BASE}/SRR17134085_mapped_sorted_RG_rmdup.bam|${BASE}/SRR17134088_mapped_sorted_RG_rmdup.bam"
    "${BASE}/SRR17134086_mapped_sorted_RG_rmdup.bam|${BASE}/SRR17134087_mapped_sorted_RG_rmdup.bam"
    "${BASE}/SRR17134087_mapped_sorted_RG_rmdup.bam|${BASE}/SRR17134088_mapped_sorted_RG_rmdup.bam"
)

# Get the pair for the current SLURM task ID (0 through 5)
CURRENT_PAIR=${PAIRS[$SLURM_ARRAY_TASK_ID]}

# Split the string into individual BAM paths
BAM_A=$(echo "$CURRENT_PAIR" | cut -d'|' -f1)
BAM_B=$(echo "$CURRENT_PAIR" | cut -d'|' -f2)

# Extract Sample IDs for file naming
ID_A=$(basename "$BAM_A" | cut -d'_' -f1)
ID_B=$(basename "$BAM_B" | cut -d'_' -f1)

echo "------------------------------------------------------"
echo "Job Array ID: ${SLURM_ARRAY_TASK_ID}"
echo "Pair: $ID_A vs $ID_B"
echo "------------------------------------------------------"

# ----------------------------
# STEP 1: Individual SAF Generation
# ----------------------------
for BAM in "$BAM_A" "$BAM_B"; do
    ID=$(basename "$BAM" | cut -d'_' -f1)
    if [ ! -f "Individual_SAFs/${ID}.saf.idx" ]; then
        echo "Generating SAF for $ID..."
        angsd -P 13 -i "$BAM" -anc "$REF" -ref "$REF" -out "Individual_SAFs/${ID}" \
              -dosaf 1 -gl 1 -C 50 -minQ 20 -minmapq 30
    else
        echo "SAF for $ID already exists. Skipping."
    fi
done

# ----------------------------
# STEP 2: Pairwise 2DSFS Generation
# ----------------------------
echo "Generating Pairwise SFS using: $REALSFS_PATH"

$REALSFS_PATH "Individual_SAFs/${ID_A}.saf.idx" "Individual_SAFs/${ID_B}.saf.idx" -P 13 \
> "Pairwise_SFS/${ID_A}_${ID_B}.real.sfs"

echo "Done. Results in Pairwise_SFS/${ID_A}_${ID_B}.real.sfs"
```
#### Step 3. Converting the .sfs to .mi.sfs files.
```bash
#!/bin/bash -l
#SBATCH --job-name=MiSTI_Convert_Gazelle
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --partition=batch
#SBATCH --array=0-5
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/SFS/sfs_to_mi_sfs/logs/misti_convert_%A_%a.log
#SBATCH --error=/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/SFS/sfs_to_mi_sfs/logs/misti_convert_%A_%a.err

# ----------------------------
# 1. PRE-FLIGHT CHECK
# ----------------------------
# Ensure the log directory exists before the script does heavy lifting
mkdir -p /home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/SFS/sfs_to_mi_sfs/logs

# ----------------------------
# 2. THE 6 PAIRS
# ----------------------------
PAIRS=(
"SRR17129394 SRR17134087"
"SRR17129394 SRR17134088"
"SRR17134085 SRR17134087"
"SRR17134085 SRR17134088"
"SRR17134086 SRR17134087"
"SRR17134086 SRR17134088"
)

# ----------------------------
# 3. PATHS & SETUP
# ----------------------------
SFS_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/SFS"
OUT_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/SFS/sfs_to_mi_sfs"
CONVERTER="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/MiSTI/utils/ANGSDSFS.py"

mkdir -p $OUT_DIR

CURRENT_PAIR=${PAIRS[$SLURM_ARRAY_TASK_ID]}
ID1=$(echo $CURRENT_PAIR | awk '{print $1}')
ID2=$(echo $CURRENT_PAIR | awk '{print $2}')

# ----------------------------
# 4. CONVERSION
# ----------------------------
INPUT_SFS="${SFS_DIR}/${ID1}_${ID2}.sfs"

echo "Processing Job Array ID: ${SLURM_ARRAY_TASK_ID}"
echo "Input: ${INPUT_SFS}"

if [ -f "$INPUT_SFS" ]; then
    # Using python3 is usually safer on clusters
    python3 $CONVERTER $INPUT_SFS $ID1 $ID2 > ${OUT_DIR}/${ID1}_${ID2}.mi.sfs
    echo "SUCCESS: Created ${OUT_DIR}/${ID1}_${ID2}.mi.sfs"
else
    echo "ERROR: File $INPUT_SFS not found!"
    exit 1
fi
```

#### Step 4. Calculate the Timescale (calc_time.py)
```bash
#!/bin/bash -l
#SBATCH --job-name=MiSTI_Timescale_Gazelle
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=batch
#SBATCH --array=0-5
#SBATCH --output=/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/SFS/sfs_to_mi_sfs/Timescale_calculation/logs/timescale_%A_%a.log
#SBATCH --error=/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/SFS/sfs_to_mi_sfs/Timescale_calculation/logs/timescale_%A_%a.err

# 1. PATHS
BASE="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD"
PSMC_DIR="${BASE}/Real_SFS/PSMC_files"
WORKDIR="${BASE}/SFS/sfs_to_mi_sfs/Timescale_calculation"
MISTI_UTILS="${BASE}/MiSTI/utils"

# 2. DEFINE PAIRS
PAIRS=(
"SRR17129394 SRR17134087"
"SRR17129394 SRR17134088"
"SRR17134085 SRR17134087"
"SRR17134085 SRR17134088"
"SRR17134086 SRR17134087"
"SRR17134086 SRR17134088"
)

CURRENT_PAIR=${PAIRS[$SLURM_ARRAY_TASK_ID]}
ID1=$(echo $CURRENT_PAIR | awk '{print $1}')
ID2=$(echo $CURRENT_PAIR | awk '{print $2}')

# 3. EXECUTION
# We MUST be in utils for the script to find the modified migrationIO.py correctly
cd $MISTI_UTILS

echo "Processing $ID1 vs $ID2 with Gazelle Defaults"

# Running without --funits to ensure it relies on the hardcoded defaults
python3 calc_time.py "${PSMC_DIR}/${ID1}.psmc" "${PSMC_DIR}/${ID2}.psmc" > "${WORKDIR}/timescale.${ID1}_${ID2}.txt"

echo "Done. Results in ${WORKDIR}/timescale.${ID1}_${ID2}.txt"
```
#### Step 5. Generating the Migration Bands
- Before running MiSTI, I need to tell when migration could have happened.
- My collaborator used specific time windows(0-1000 years ago).
- This script reads the "timescale.txt files" and creates a small helper file for each pair.
```bash
#!/bin/bash
# ----------------------------
# Create Migration Band Definitions
# ----------------------------
BASE="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD"
SFS_DIR="${BASE}/Pairwise_SFS"

PAIRS=("SRR17129394_SRR17134087" "SRR17129394_SRR17134088" "SRR17134085_SRR17134087" "SRR17134085_SRR17134088" "SRR17134086_SRR17134087" "SRR17134086_SRR17134088")

for pair in "${PAIRS[@]}"; do
    # This looks at the timescale file and finds which PSMC steps correspond to 0-1000 years
    # It outputs a string like "-mi 1 0 3 0.00 1 -mi 2 0 3 0.00 1"
    awk '$2>=0 && $2<1000{print $1}' ${SFS_DIR}/timescale.${pair//_/.}.txt | sed -e 1b -e '$!d' | awk -v ORS=" " '{print $0}' | awk '{print "-mi 1",$1,$2,"0.00 1 -mi 2",$1,$2,"0.00 1"}' > ${SFS_DIR}/${pair}.migBand
done
```
#### Step 6. Final MiSTI optimization
- It takes your .mi.sfs, the psmc files and the migration bands and then tests the 55 split times for each of my species six pairs. 
```bash
#!/bin/bash -l
#SBATCH --job-name=Dama_MiSTI_Final
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=13
#SBATCH --mem=40G
#SBATCH --partition=batch
#SBATCH --array=0-5
#SBATCH --output=/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/logs/misti_final_%A_%a.log

# 1. Setup Pair
PAIRS=("SRR17129394 SRR17134087" "SRR17129394 SRR17134088" "SRR17134085 SRR17134087" "SRR17134085 SRR17134088" "SRR17134086 SRR17134087" "SRR17134086 SRR17134088")
CURRENT_PAIR=${PAIRS[$SLURM_ARRAY_TASK_ID]}
ID1=$(echo $CURRENT_PAIR | awk '{print $1}')
ID2=$(echo $CURRENT_PAIR | awk '{print $2}')

# 2. Paths
BASE="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD"
PSMC_DIR="${BASE}/Real_SFS/PSMC_files"
SFS_DIR="${BASE}/Pairwise_SFS"
MISTI_PY="${BASE}/MiSTI/misti/misti.py"
UNITS="${BASE}/MiSTI/misti/setunits.txt"
RESULT_DIR="${BASE}/Results"

mkdir -p $RESULT_DIR

# 3. RUN 55 SPLIT POINTS IN PARALLEL
# We use the migration bands defined by your collaborator's leopard/wolf logic
# -j 13 means it runs 13 split points at a time on your 13 allocated CPUs
parallel --header : -j 13 python ${MISTI_PY} \
    -uf --bsSize 10 --hetloss 0.0 0.0 --funits ${UNITS} \
    ${PSMC_DIR}/${ID1}.psmc ${PSMC_DIR}/${ID2}.psmc \
    ${SFS_DIR}/${ID1}_${ID2}.mi.sfs {per} \
    -o ${RESULT_DIR}/${ID1}_${ID2}_{per}.mi \
    -mi 1 0 1 0.00 1 -mi 2 0 1 0.00 1 \
    -mi 1 2 14 0.00 1 -mi 2 2 14 0.00 1 \
    -mi 1 15 18 0.00 1 -mi 2 15 18 0.00 1 \
    -mi 1 19 23 0.00 1 -mi 2 19 23 0.00 1 \
    -mi 1 24 31 0.00 1 -mi 2 24 31 0.00 1 \
    -mi 1 32 33 0.00 1 -mi 2 32 33 0.00 1 \
    -mi 1 34 41 0.00 1 -mi 2 34 41 0.00 1 \
    -mi 1 42 51 0.00 1 -mi 2 42 51 0.00 1 \
    -mi 1 52 55 0.00 1 -mi 2 52 55 0.00 1 \
    ">>" ${RESULT_DIR}/${ID1}_${ID2}.opt.out ::: per $(seq 1 55)
```
#### Step 7.  Final Step: Extracting the Data for Excel 
- Once the job in the Step 6 finishes, I will have teh file named Results/SRR17129394_SRR17134087.opt.out.
- You need to extract the Likelihoods and Years to find the peak split time.
```bash
# Example for one pair
echo "Step Year Likelihood"
paste <(awk '{print $1, $2}' /home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Pairwise_SFS/timescale.SRR17129394.SRR17134087.txt | grep -v "Units") \
      <(grep "llh =" /home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Results/SRR17129394_SRR17134087.opt.out | awk '{print $15}') \
      | column -t
```
---
- Copy the output from the Step 7 into excel and create a scatterplot and add a 5th-degree polynomial trendline.
- The peak of that line is the estimted divergence time
---

