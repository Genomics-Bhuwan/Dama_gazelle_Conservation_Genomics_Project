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

