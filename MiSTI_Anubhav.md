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
#SBATCH --job-name=Dama_6Pairs_Hardcoded
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=13
#SBATCH --mem=32G
#SBATCH --partition=batch
#SBATCH --array=0-5
#SBATCH --output=logs/dama_pair_%A_%a.log

# ----------------------------
# 1. DEFINE YOUR 6 PAIRS MANUALLY HERE
# ----------------------------
# Format: "BAM_A BAM_B"
PAIRS=(
"/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Haplotype_Caller/Five_samples/rmdup_Dama/SRR17129394_mapped_sorted_RG_rmdup.bam /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Haplotype_Caller/Five_samples/rmdup_Dama/SRR17134087_mapped_sorted_RG_rmdup.bam"
"/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Haplotype_Caller/Five_samples/rmdup_Dama/SRR17129394_mapped_sorted_RG_rmdup.bam /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Haplotype_Caller/Five_samples/rmdup_Dama/SRR17134088_mapped_sorted_RG_rmdup.bam"
"/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Haplotype_Caller/Five_samples/rmdup_Dama/SRR17134085_mapped_sorted_RG_rmdup.bam /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Haplotype_Caller/Five_samples/rmdup_Dama/SRR17134087_mapped_sorted_RG_rmdup.bam"
"/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Haplotype_Caller/Five_samples/rmdup_Dama/SRR17134085_mapped_sorted_RG_rmdup.bam /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Haplotype_Caller/Five_samples/rmdup_Dama/SRR17134088_mapped_sorted_RG_rmdup.bam"
"/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Haplotype_Caller/Five_samples/rmdup_Dama/SRR17134086_mapped_sorted_RG_rmdup.bam /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Haplotype_Caller/Five_samples/rmdup_Dama/SRR17134087_mapped_sorted_RG_rmdup.bam"
"/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Haplotype_Caller/Five_samples/rmdup_Dama/SRR17134086_mapped_sorted_RG_rmdup.bam /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Haplotype_Caller/Five_samples/rmdup_Dama/SRR17134088_mapped_sorted_RG_rmdup.bam"
)

# ----------------------------
# 2. SELECT THE CURRENT PAIR
# ----------------------------
# SLURM_ARRAY_TASK_ID will be 0, 1, 2, 3, 4, or 5
CURRENT_PAIR=${PAIRS[$SLURM_ARRAY_TASK_ID]}
BAM_A=$(echo $CURRENT_PAIR | awk '{print $1}')
BAM_B=$(echo $CURRENT_PAIR | awk '{print $2}')

# Set IDs (Strip the long suffix)
ID_A=$(basename $BAM_A _mapped_sorted_RG_rmdup.bam)
ID_B=$(basename $BAM_B _mapped_sorted_RG_rmdup.bam)

# ----------------------------
# 3. SET PATHS
# ----------------------------
OUT_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD"
REF="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Dama_gazelle_hifiasm-ULONT_primary.fasta"

mkdir -p $OUT_DIR/Individual_SAFs $OUT_DIR/Pairwise_SFS $OUT_DIR/logs
cd $OUT_DIR

echo "Running Pair Index $SLURM_ARRAY_TASK_ID: $ID_A vs $ID_B"

# ----------------------------
# 4. STEP 1: SAF Generation
# ----------------------------
for BAM in $BAM_A $BAM_B; do
    ID=$(basename $BAM _mapped_sorted_RG_rmdup.bam)
    if [ ! -f Individual_SAFs/${ID}.saf.idx ]; then
        angsd -P 13 -i $BAM -anc $REF -ref $REF -out Individual_SAFs/${ID} \
              -dosaf 1 -gl 1 -C 50 -minQ 20 -minmapq 30
    fi
done

# ----------------------------
# 5. STEP 2: 2DSFS Generation
# ----------------------------
realSFS Individual_SAFs/${ID_A}.saf.idx Individual_SAFs/${ID_B}.saf.idx -P 13 \
> Pairwise_SFS/${ID_A}_${ID_B}.real.sfs

echo "Completed SFS for $ID_A and $ID_B"
```
#### Step 3. Converting the .sfs to .mi.sfs files.
```bash
#!/bin/bash -l
#SBATCH --job-name=MiSTI_Convert_Hardcoded
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=40G
#SBATCH --partition=batch
#SBATCH --array=0-5
#SBATCH --output=/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/logs/misti_convert_%A_%a.log

# ----------------------------
# 1. HARDCODED PAIRS (Same as Step 1)
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
# 2. PATHS & SETUP
# ----------------------------
BASE_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD"
SFS_DIR="${BASE_DIR}/Pairwise_SFS"
CONVERTER="${BASE_DIR}/MiSTI/utils/ANGSDSFS.py"

# Select the pair based on the Task ID (0-5)
CURRENT_PAIR=${PAIRS[$SLURM_ARRAY_TASK_ID]}
ID1=$(echo $CURRENT_PAIR | awk '{print $1}')
ID2=$(echo $CURRENT_PAIR | awk '{print $2}')

cd $SFS_DIR

# ----------------------------
# 3. CONVERSION
# ----------------------------
# Filename from Step 1 was: ID1_ID2.real.sfs
INPUT_SFS="${ID1}_${ID2}.real.sfs"

echo "------------------------------------------------------"
echo "Job Array ID: ${SLURM_ARRAY_TASK_ID}"
echo "Converting SFS for Pair: $ID1 and $ID2"
echo "Input File: $INPUT_SFS"
echo "------------------------------------------------------"

if [ -f "$INPUT_SFS" ]; then
    # Run the python converter
    # Argument 1: The raw SFS file
    # Argument 2: Name of Individual 1 (Matches PSMC filename)
    # Argument 3: Name of Individual 2 (Matches PSMC filename)
    python $CONVERTER $INPUT_SFS $ID1 $ID2 > ${ID1}_${ID2}.mi.sfs
    echo "Successfully created ${ID1}_${ID2}.mi.sfs"
else
    echo "ERROR: Input file $INPUT_SFS not found in $SFS_DIR"
    exit 1
fi

echo "Done!"
```

#### Step 4. Calculate the Timescale (calc_time.py)
```bash
#!/bin/bash -l
#SBATCH --job-name=Dama_Step2_Prep
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=batch
#SBATCH --array=1-6
#SBATCH --output=/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/logs/prep_%A_%a.log

# ----------------------------
# 4. Calculate Timescale (The Step You Asked For)
# ----------------------------
#!/bin/bash -l
#SBATCH --job-name=MiSTI_Timescale
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=batch
#SBATCH --array=0-5
#SBATCH --output=/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/logs/timescale_%A_%a.log

# ----------------------------
# 1. HARDCODED PAIRS (Matches your 6 specific comparisons)
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
# 2. FULL PATHS
# ----------------------------
BASE="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD"
PSMC_DIR="${BASE}/Real_SFS/PSMC_files"
SFS_DIR="${BASE}/Pairwise_SFS"
CALC_PY="${BASE}/MiSTI/utils/calc_time.py"
UNITS="${BASE}/MiSTI/misti/setunits.txt"

# Select the pair
CURRENT_PAIR=${PAIRS[$SLURM_ARRAY_TASK_ID]}
ID1=$(echo $CURRENT_PAIR | awk '{print $1}')
ID2=$(echo $CURRENT_PAIR | awk '{print $2}')

echo "------------------------------------------------------"
echo "Calculating Timescale for: $ID1 vs $ID2"
echo "PSMC 1: ${PSMC_DIR}/${ID1}.psmc"
echo "PSMC 2: ${PSMC_DIR}/${ID2}.psmc"
echo "------------------------------------------------------"

# ----------------------------
# 3. RUN CALC_TIME.PY
# ----------------------------
# This generates the text file mapping steps to years
python $CALC_PY --funits $UNITS \
    ${PSMC_DIR}/${ID1}.psmc \
    ${PSMC_DIR}/${ID2}.psmc \
    > ${SFS_DIR}/timescale.${ID1}.${ID2}.txt

if [ $? -eq 0 ]; then
    echo "Successfully generated timescale.${ID1}.${ID2}.txt"
else
    echo "ERROR: Calculation failed for $ID1 and $ID2"
    exit 1
fi
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

