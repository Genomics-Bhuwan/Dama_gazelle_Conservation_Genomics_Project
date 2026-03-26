#### Using MiSTI for divergence between two species.
- https://github.com/Genomics-HSE/MiSTI

#### Step 1. Run PSMC (Li, Durbin 2011) on two genomes.
#### Step 2 a. Run ANGSD and SFS to get it and then combine. Generate joint SFS
```bash
#!/bin/bash -l
#SBATCH --job-name=Dama_SAF_Only
#SBATCH --time=54:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --partition=batch
#SBATCH --array=0-4
#SBATCH --output=logs/saf_gen_%A_%a.log

module load angsd

OUT_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD"
REF="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Dama_gazelle_hifiasm-ULONT_primary.fasta"
BASE="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Haplotype_Caller/Five_samples/rmdup_Dama"

mkdir -p ${OUT_DIR}/Individual_SAFs

# The 5 unique BAM files from your list
BAM_FILES=(
    "${BASE}/SRR17129394_mapped_sorted_RG_rmdup.bam"
    "${BASE}/SRR17134085_mapped_sorted_RG_rmdup.bam"
    "${BASE}/SRR17134086_mapped_sorted_RG_rmdup.bam"
    "${BASE}/SRR17134087_mapped_sorted_RG_rmdup.bam"
    "${BASE}/SRR17134088_mapped_sorted_RG_rmdup.bam"
)

CURRENT_BAM=${BAM_FILES[$SLURM_ARRAY_TASK_ID]}
ID=$(basename "$CURRENT_BAM" | cut -d'_' -f1)

echo "Processing Sample: $ID"

# Generate SAF
angsd -P 14 \
    -i "$CURRENT_BAM" \
    -anc "$REF" \
    -ref "$REF" \
    -out "${OUT_DIR}/Individual_SAFs/${ID}" \
    -dosaf 1 \
    -gl 1 \
    -C 50 \
    -minQ 20 \
    -minmapq 30

```
#### Step 2 b. SFS
```bash
#!/bin/bash -l
#SBATCH --job-name=Dama_6Pairs_SFS
#SBATCH --time=54:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=84G
#SBATCH --partition=batch
#SBATCH --array=0-5
#SBATCH --output=logs/sfs_pair_%A_%a.log
#SBATCH --error=logs/sfs_pair_%A_%a.err

# --- Environment ---
module load angsd

# PATHS
REALSFS="/software/ngsTools/1.0.2/ngsTools/ANGSD/angsd/misc/realSFS"
SAF_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Individual_SAFs"
OUT_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/SFS"

mkdir -p $OUT_DIR/logs

# --- Define the 6 Cross-Population Pairs ---
PAIRS=(
  "SRR17129394|SRR17134087"
  "SRR17129394|SRR17134088"
  "SRR17134085|SRR17134087"
  "SRR17134085|SRR17134088"
  "SRR17134086|SRR17134087"
  "SRR17134086|SRR17134088"
)

# Get the pair for THIS array task
CURRENT_PAIR=${PAIRS[$SLURM_ARRAY_TASK_ID]}
ID_A=$(echo "$CURRENT_PAIR" | cut -d'|' -f1)
ID_B=$(echo "$CURRENT_PAIR" | cut -d'|' -f2)

echo "Calculating SFS for: Addra ($ID_A) vs Mhorr ($ID_B)"

# --- Run realSFS with -nSites to prevent memory crash ---
# -P 16 matches cpus-per-task
# -nSites 50000000 limits memory usage by processing in blocks
$REALSFS ${SAF_DIR}/${ID_A}.saf.idx ${SAF_DIR}/${ID_B}.saf.idx \
-P 16 -nSites 50000000 > ${OUT_DIR}/${ID_A}_${ID_B}.sfs

echo "Finished Pair: $ID_A vs $ID_B"

```
#### Step 3. Convert joint realSFS to MiSTI file format.
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

#### Step 4. Align the Time Scales
- Since, we have 3 Addra vs 2 Mhorr. We assesed the index time interval for each pair.
- MiSTI merges the discrete time windows from your two separate PSMC runs into a single timeline.
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

#### Step 5. Generating the migration band for six pairs 
- It creates a .migBand file for every pair, setting the migration rate to 0.00 (testing for isolation).
-mi 1: Migration from Addra to Mhorr.
-mi 2: Migration from Mhorr to Addra.
- 0.00: We are testing if they were totally isolated during these blocks.
- 1: This tells MiSTI that this parameter is "active" for the model.
```bash
#!/bin/bash -l
#SBATCH --job-name=MiSTI_Bands_Gazelle
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=batch
#SBATCH --output=/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/SFS/sfs_to_mi_sfs/Timescale_calculation/Migration_band/logs/bands_%j.log

# 1. Define Paths
TIMESCALE_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/SFS/sfs_to_mi_sfs/Timescale_calculation"
MIGBAND_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/SFS/sfs_to_mi_sfs/Timescale_calculation/Migration_band"

mkdir -p $MIGBAND_DIR
cd $TIMESCALE_DIR

# 2. Define the pairs
PAIRS=(
  "SRR17129394 SRR17134087"
  "SRR17129394 SRR17134088"
  "SRR17134085 SRR17134087"
  "SRR17134085 SRR17134088"
  "SRR17134086 SRR17134087"
  "SRR17134086 SRR17134088"
)

# 3. Loop to create bands based on verified Gazelle history indices
for p in "${PAIRS[@]}"; do
    ID_A=$(echo "$p" | awk '{print $1}')
    ID_B=$(echo "$p" | awk '{print $2}')
    OUTFILE="${MIGBAND_DIR}/${ID_A}_${ID_B}.migBand"

    echo "Creating migration bands for: $ID_A vs $ID_B"
    
    # Reset/Create the file
    > $OUTFILE

    # --- BLOCK 1: Modern Era (0 - 5,000 years ago) ---
    # WHY INDICES 0-1: Index 1 ends at ~4,000 years in your data. 
    # This captures the period after the Sahara dried out, forcing 
    # gazelles into fragmented populations.
    echo "# Block 1: Modern Isolation (0-5kya)" >> $OUTFILE
    echo "-mi 1 0 1 0.00 1" >> $OUTFILE
    echo "-mi 2 0 1 0.00 1" >> $OUTFILE

    # --- BLOCK 2: African Humid Period / Green Sahara (5,000 - 11,000 years ago) ---
    # WHY INDICES 2-3: Index 3 ends at ~8,000-8,400 years. 
    # This captures the "Green Sahara" where high connectivity was 
    # ecologically possible due to increased rainfall.
    echo "# Block 2: Green Sahara (5k-11kya)" >> $OUTFILE
    echo "-mi 1 2 3 0.00 1" >> $OUTFILE
    echo "-mi 2 2 3 0.00 1" >> $OUTFILE

    # --- BLOCK 3: Last Glacial Maximum (11,000 - 25,000 years ago) ---
    # WHY INDICES 4-7: Index 7 ends at ~23,000-24,000 years. 
    # This is the hyper-arid period and the most likely time for the 
    # original genetic split between Addra and Mhorr.
    echo "# Block 3: LGM Split (11k-25kya)" >> $OUTFILE
    echo "-mi 1 4 7 0.00 1" >> $OUTFILE
    echo "-mi 2 4 7 0.00 1" >> $OUTFILE

    # --- BLOCK 4: Deep Ancestry (25,000 - 100,000 years ago) ---
    # WHY INDICES 8-21: Index 21 ends at ~90,000-96,000 years. 
    # This allows MiSTI to account for ancestral population size 
    # fluctuations before the primary subspecies divergence.
    echo "# Block 4: Deep Ancestry (25k-100kya)" >> $OUTFILE
    echo "-mi 1 8 21 0.00 1" >> $OUTFILE
    echo "-mi 2 8 21 0.00 1" >> $OUTFILE

done

echo "✅ All 6 migration band files (.migBand) are ready with calibrated indices!"
```
#### Step 6. Verify the MAIN program.
- Before running the main program, let's verify the files are correct.
- Every .migBand file should have 8 lines (2 lines per time block).
```bash
wc -l *.migBand
```
#### Step 7. 
- We are running a GRID search.
- We will test 60 different "Split times from 1 to 60" for each of your 6 pairs.
- The goal is to find which specific year in history gives the highest likelihood(llh) score- meaning it is the most statistically probable time that Addra and Mhorr gazelles diverged.
- Since you have high-depth data (25-30x), I have set the --hetloss (False Negative Rate) to 0.02 (a very safe 2% correction since we have high quality resequenced samples).
```bash
#!/bin/bash -l
#SBATCH --job-name=MiSTI_Final_Fix
#SBATCH --time=58:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH --partition=batch
#SBATCH --array=0-5
#SBATCH --output=misti_fix_%A_%a.log

cd /home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/MiSTI

PAIRS=(
  "SRR17129394 SRR17134087"
  "SRR17129394 SRR17134088"
  "SRR17134085 SRR17134087"
  "SRR17134085 SRR17134088"
  "SRR17134086 SRR17134087"
  "SRR17134086 SRR17134088"
)

ID_A=$(echo ${PAIRS[$SLURM_ARRAY_TASK_ID]} | awk '{print $1}')
ID_B=$(echo ${PAIRS[$SLURM_ARRAY_TASK_ID]} | awk '{print $2}')

PSMC_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files"
SFS_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/SFS/sfs_to_mi_sfs"
MIG_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/SFS/sfs_to_mi_sfs/Timescale_calculation/Migration_band"
OUT_DIR="${MIG_DIR}/MiSTI_Final_Results"

mkdir -p "$OUT_DIR"
FINAL_LOG="${OUT_DIR}/final_results_${ID_A}_${ID_B}.log"

for splitTime in {1..80}
do
    # 1. Clean migration parameters
    MIG_PARAMS=$(grep -v '^#' "${MIG_DIR}/${ID_A}_${ID_B}.migBand")

    # 2. Run Python command
    # Removed --bsSize because it is not in your MiSTI.py script
    python3 MiSTI.py \
        -o "${OUT_DIR}/Result_${ID_A}_${ID_B}_${splitTime}.mi" \
        -uf \
        --hetloss 0.02 0.02 \
        $MIG_PARAMS \
        "${PSMC_DIR}/${ID_A}.psmc" \
        "${PSMC_DIR}/${ID_B}.psmc" \
        "${SFS_DIR}/${ID_A}_${ID_B}.mi.sfs" \
        "$splitTime" >> "$FINAL_LOG" 2>&1

    echo "Finished Index ${splitTime}" >> "$FINAL_LOG"
done
```
#### Step 8: Extract the likelihoods:
```bash
grep "llh =" final_results_SRR17129394_SRR17134087.log
```
