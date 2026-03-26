#### Using MiSTI for divergence between two species.
- https://github.com/Genomics-HSE/MiSTI

#### Step 1. Run PSMC (Li, Durbin 2011) on two genomes.
#### Step 2. Run ANGSD and SFS to get it and then combine. Generate joint SFS
```bash
#!/bin/bash -l
#SBATCH --job-name=ANGSD_MiSTI
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --partition=batch
#SBATCH --output=angsd_misti_%A.log

# ----------------------------
# Paths & References
# ----------------------------
OUT_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD"
REF="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Dama_gazelle_hifiasm-ULONT_primary.fasta"

cd $OUT_DIR

# ----------------------------
# STEP 1: Generate SAF files for each population
# ----------------------------
echo "[$(date)] Generating SAF for Addra..."
angsd -P 24 -bam addra_bams.txt -anc $REF -ref $REF -out addra_pop \
      -dosaf 1 -gl 1 -C 50 -minQ 20 -minmapq 30

echo "[$(date)] Generating SAF for Mhorr..."
angsd -P 24 -bam mhorr_bams.txt -anc $REF -ref $REF -out mhorr_pop \
      -dosaf 1 -gl 1 -C 50 -minQ 20 -minmapq 30

# ----------------------------
# STEP 2: Generate Joint SFS (realSFS)
# ----------------------------
#!/bin/bash -l
#SBATCH --job-name=Dama_6Pairs_SFS
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
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
# Format: Addra|Mhorr
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

# --- Run realSFS ---
# -P 16 uses 16 cores to speed up the calculation
$REALSFS ${SAF_DIR}/${ID_A}.saf.idx ${SAF_DIR}/${ID_B}.saf.idx -P 16 > ${OUT_DIR}/${ID_A}_${ID_B}.sfs

echo "Finished Pair: $ID_A vs $ID_B"
```
#### Step 3. Convert joint realSFS to MiSTI file format.
```bash
# 1. Define the pairs (Addra vs Mhorr)
PAIRS=(
  "SRR17129394|SRR17134087"
  "SRR17129394|SRR17134088"
  "SRR17134085|SRR17134087"
  "SRR17134085|SRR17134088"
  "SRR17134086|SRR17134087"
  "SRR17134086|SRR17134088"
)

# 2. Run the conversion loop
for p in "${PAIRS[@]}"; do
    # Extract the individual IDs from the pair string
    ID_A=$(echo "$p" | cut -d'|' -f1)
    ID_B=$(echo "$p" | cut -d'|' -f2)

    echo "Converting Pair: $ID_A vs $ID_B..."

    # Use the MiSTI utility to format the SFS
    # Syntax: python ANGSDSFS.py [input_sfs] [pop1_name] [pop2_name]
    python ANGSDSFS.py ${ID_A}_${ID_B}.sfs $ID_A $ID_B > ${ID_A}_${ID_B}.mi.sfs
done

echo "✅ All 6 pairs converted to .mi.sfs format."
```

#### Step 4. Align the Time Scales
- Since, we have 3 Addra vs 2 Mhorr. We assesed the index time interval for each pair.
- MiSTI merges the discrete time windows from your two separate PSMC runs into a single timeline.
```bash
# 1. Define the pairs (Addra vs Mhorr)
PAIRS=(
  "SRR17129394|SRR17134087"
  "SRR17129394|SRR17134088"
  "SRR17134085|SRR17134087"
  "SRR17134085|SRR17134088"
  "SRR17134086|SRR17134087"
  "SRR17134086|SRR17134088"
)

# 2. Path to your PSMC files
PSMC_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files"

# 3. Create the units file (Gazelle mutation rate and generation time)
echo "2.96e-09 5.85" > gazelle_units.txt

# 4. Loop to calculate the timescale for every pair
for p in "${PAIRS[@]}"; do
    ID_A=$(echo "$p" | cut -d'|' -f1)
    ID_B=$(echo "$p" | cut -d'|' -f2)

    echo "Calculating timescale for: $ID_A vs $ID_B"
    
    # Run calc_time.py using your specific PSMC files
    python calc_time.py --funits gazelle_units.txt \
    ${PSMC_DIR}/${ID_A}.psmc ${PSMC_DIR}/${ID_B}.psmc > timescale.${ID_A}_${ID_B}.txt
done

echo "✅ Timescale generation complete for all 6 pairs!"
```

#### Step 5. Generating the migration band for six pairs 
- It creates a .migBand file for every pair, setting the migration rate to 0.00 (testing for isolation).
-mi 1: Migration from Addra to Mhorr.
-mi 2: Migration from Mhorr to Addra.
- 0.00: We are testing if they were totally isolated during these blocks.
- 1: This tells MiSTI that this parameter is "active" for the model.
```bash
# 1. Define the pairs
PAIRS=(
  "SRR17129394|SRR17134087"
  "SRR17129394|SRR17134088"
  "SRR17134085|SRR17134087"
  "SRR17134085|SRR17134088"
  "SRR17134086|SRR17134087"
  "SRR17134086|SRR17134088"
)

# 2. Loop to create bands based on Gazelle history
for p in "${PAIRS[@]}"; do
    ID_A=$(echo "$p" | cut -d'|' -f1)
    ID_B=$(echo "$p" | cut -d'|' -f2)
    TIMESCALE="timescale.${ID_A}_${ID_B}.txt"
    OUTFILE="${ID_A}_${ID_B}.migBand"

    echo "Creating migration bands for: $ID_A vs $ID_B"
    
    # Reset/Create the file
    > $OUTFILE

    # --- BLOCK 1: Modern Era (0 - 5,000 years ago) ---
    # Finds the start and end index in the timescale for these years
    # BLOCK 1: 0 - 5,000 years ago (The Late Holocene Aridification) ---
    # Context: The "Desertification of the Sahara." After the Green Sahara ended, 
    # populations were forced into "refugia" (pockets of vegetation), 
    # likely leading to the modern isolation of Addra and Mhorr.
    awk '$2>=0 && $2<5000{print $1}' $TIMESCALE | sed -e 1b -e '$!d' | awk -v ORS=" " '{print $0}' | awk '{print "-mi 1",$0,"0.00 1\n-mi 2",$0,"0.00 1"}' >> $OUTFILE

    # --- BLOCK 2: African Humid Period / Green Sahara (5,000 - 11,000 years ago) ---
     Context: Increased monsoon rainfall created a "Savannah-like" Sahara with lakes and rivers. 
    # High connectivity was possible. We test 0.00 migration here to see if 
    # they were already separate species even when the habitat was open.
    awk '$2>=5000 && $2<11000{print $1}' $TIMESCALE | sed -e 1b -e '$!d' | awk -v ORS=" " '{print $0}' | awk '{print "-mi 1",$0,"0.00 1\n-mi 2",$0,"0.00 1"}' >> $OUTFILE

    # --- BLOCK 3: Last Glacial Maximum / Hyper-Arid (11,000 - 25,000 years ago) ---
    Context: A hyper-arid period. The Sahara was even larger and drier than today. 
    # This is a strong candidate for the original "Split Time" where 
    # northern (Mhorr) and eastern (Addra) populations were cut off by sand seas.
    awk '$2>=11000 && $2<25000{print $1}' $TIMESCALE | sed -e 1b -e '$!d' | awk -v ORS=" " '{print $0}' | awk '{print "-mi 1",$0,"0.00 1\n-mi 2",$0,"0.00 1"}' >> $OUTFILE

    # --- BLOCK 4: Early Pleistocene (25,000 - 100,000 years ago) ---
     # Context: Earlier cycles of wet and dry periods. This helps MiSTI understand 
    # the "Deep Ancestry" and the size of the original ancestral population 
    # before the subspecies split.
    awk '$2>=25000 && $2<100000{print $1}' $TIMESCALE | sed -e 1b -e '$!d' | awk -v ORS=" " '{print $0}' | awk '{print "-mi 1",$0,"0.00 1\n-mi 2",$0,"0.00 1"}' >> $OUTFILE
done

echo "✅ All 6 migration band files (.migBand) are ready!"
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
- Since you have high-depth data (25-30x), I have set the --hetloss (False Negative Rate) to 0.05 (a very safe 5% correction).
```bash
#!/bin/bash -l
#SBATCH --job-name=Gazelle_MiSTI
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --partition=batch
#SBATCH --array=0-5
#SBATCH --output=logs/misti_%A_%a.log
#SBATCH --error=logs/misti_%A_%a.err

# --- Define the Pairs ---
PAIRS=(
  "SRR17129394|SRR17134087"
  "SRR17129394|SRR17134088"
  "SRR17134085|SRR17134087"
  "SRR17134085|SRR17134088"
  "SRR17134086|SRR17134087"
  "SRR17134086|SRR17134088"
)

# Get current pair IDs
CURRENT_PAIR=${PAIRS[$SLURM_ARRAY_TASK_ID]}
ID_A=$(echo "$CURRENT_PAIR" | cut -d'|' -f1)
ID_B=$(echo "$CURRENT_PAIR" | cut -d'|' -f2)

# PATHS to your files
PSMC_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files"
MIG_BAND="${ID_A}_${ID_B}.migBand"
SFS_FILE="${ID_A}_${ID_B}.mi.sfs"

echo "Running MiSTI Grid Search for $ID_A vs $ID_B"

# --- The Grid Search Loop ---
# We test split times 1 through 60. 
# Each result is appended to a master log file for that pair.

for splitTime in {1..60}
do
  echo "Testing Split Time Index: ${splitTime}"
  
  python MiSTI.py -uf --bsSize 10 --hetloss 0.05 0.05 \
  --funits gazelle_units.txt \
  ${PSMC_DIR}/${ID_A}.psmc ${PSMC_DIR}/${ID_B}.psmc \
  ${SFS_FILE} ${splitTime} \
  $(cat ${MIG_BAND}) \
  -o Result_${ID_A}_${ID_B}_${splitTime}.mi >> final_results_${ID_A}_${ID_B}.log
done

echo "✅ MiSTI Analysis Complete for $ID_A vs $ID_B"
```
#### Step 8: Extract the likelihoods:
```bash
grep "llh =" final_results_SRR17129394_SRR17134087.log
```
