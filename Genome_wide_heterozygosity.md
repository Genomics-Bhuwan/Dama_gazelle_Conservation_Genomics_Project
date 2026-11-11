```bash
#!/bin/bash -l
#SBATCH --job-name=Dama_ANGSD
#SBATCH --array=0-4
#SBATCH --time=80:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G
#SBATCH --partition=batch
#SBATCH --output=Dama_ANGSD_%A_%a.log
#SBATCH --error=Dama_ANGSD_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

set -euo pipefail

# ============================================================
# 1. LOAD ANGSD & SAMTOOLS ON LOCAL HPC
# ============================================================

module load gcc-6.3.0
module load angsd
module load samtools-1.22.1

echo "Software modules loaded:"
module list

echo "Testing ANGSD availability..."
angsd -GL 1 > /dev/null 2>&1 || true
realSFS > /dev/null 2>&1 || true
thetaStat > /dev/null 2>&1 || true

# ============================================================
# 2. DEFINE DIRECTORIES
# ============================================================

BASE_DIR="/home/bistbs/downsampled_10X"
BAM_DIR="${BASE_DIR}"
REFERENCE="${BASE_DIR}/Dama_gazelle_hifiasm-ULONT_primary.fasta"
OUTPUT_DIR="${BASE_DIR}/ANGSD"

mkdir -p "${OUTPUT_DIR}"

# ============================================================
# 3. DEFINE SAMPLES AND SELECT CURRENT ARRAY TASK
# ============================================================

SAMPLES=(
    "SRR17129394_10X"
    "SRR17134085_10X"
    "SRR17134086_10X"
    "SRR17134087_10X"
    "SRR17134088_10X"
)

# Pick current sample based on SLURM array index (0-4)
SAMPLE="${SAMPLES[${SLURM_ARRAY_TASK_ID}]}"
THREADS="${SLURM_CPUS_PER_TASK:-15}"

echo "================================================"
echo "RUNNING ARRAY TASK ID: ${SLURM_ARRAY_TASK_ID}"
echo "PROCESSING SAMPLE: ${SAMPLE}"
echo "ALLOCATED THREADS: ${THREADS}"
echo "Date: $(date)"
echo "================================================"

# ============================================================
# 4. CHECK INPUT FILES FOR THIS SAMPLE
# ============================================================

echo "Checking reference genome..."

if [[ ! -s "$REFERENCE" ]]; then
    echo "ERROR: Reference genome missing at $REFERENCE"
    exit 1
fi

if [[ ! -s "${REFERENCE}.fai" ]]; then
    echo "ERROR: Reference index missing!"
    echo "Run 'samtools faidx ${REFERENCE}' first."
    exit 1
fi

BAM="${BAM_DIR}/${SAMPLE}.bam"

if [[ ! -s "$BAM" ]]; then
    echo "ERROR: BAM missing: $BAM"
    exit 1
fi

if [[ ! -s "${BAM}.bai" ]]; then
    echo "ERROR: BAM index missing: ${BAM}.bai"
    exit 1
fi

echo "BAM and index found: $BAM"

# ============================================================
# 5. CHECK BAM INTEGRITY AND SCAFFOLD COMPATIBILITY
# ============================================================

if command -v samtools >/dev/null 2>&1; then

    cut -f1 "${REFERENCE}.fai" | sort -u > "${OUTPUT_DIR}/reference_scaffolds.txt"

    echo "Checking BAM integrity: $SAMPLE"

    samtools quickcheck -v "$BAM"

    samtools view -H "$BAM" | \
        awk -F'\t' '
        $1=="@SQ" {
            for(i=2;i<=NF;i++) {
                if($i ~ /^SN:/) {
                    sub(/^SN:/,"",$i)
                    print $i
                }
            }
        }' | sort -u > "${OUTPUT_DIR}/${SAMPLE}_bam_scaffolds.txt"

    comm -12 \
        "${OUTPUT_DIR}/reference_scaffolds.txt" \
        "${OUTPUT_DIR}/${SAMPLE}_bam_scaffolds.txt" \
        > "${OUTPUT_DIR}/${SAMPLE}_common_scaffolds.txt"

    if [[ ! -s "${OUTPUT_DIR}/${SAMPLE}_common_scaffolds.txt" ]]; then
        echo "ERROR: No matching scaffolds for $SAMPLE"
        exit 1
    fi

    echo "Scaffold compatibility passed: $SAMPLE"

else
    echo "WARNING: samtools not available. Skipping optional BAM checks."
fi

# ============================================================
# 6. RUN ANGSD FOR ASSIGNED SAMPLE
# ============================================================

SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE}"
mkdir -p "${SAMPLE_DIR}"
PREFIX="${SAMPLE_DIR}/${SAMPLE}"

# --------------------------------------------------------
# 6A. GENOTYPE LIKELIHOODS AND SAF ESTIMATION
# --------------------------------------------------------
echo "Step 1: Estimating genotype likelihoods and SAF for ${SAMPLE}..."

angsd \
    -P "${THREADS}" \
    -i "$BAM" \
    -ref "$REFERENCE" \
    -anc "$REFERENCE" \
    -GL 1 \
    -doSaf 1 \
    -doCounts 1 \
    -baq 1 \
    -C 50 \
    -minQ 30 \
    -minMapQ 30 \
    -minInd 1 \
    -out "$PREFIX"

if [[ ! -s "${PREFIX}.saf.idx" ]]; then
    echo "ERROR: SAF index missing for $SAMPLE"
    exit 1
fi

echo "SAF estimation completed."

#!/bin/bash -l
#SBATCH --job-name=Dama_SFS_Theta
#SBATCH --array=0-4
#SBATCH --time=80:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G
#SBATCH --partition=batch
#SBATCH --output=Dama_SFS_Theta_%A_%a.log
#SBATCH --error=Dama_SFS_Theta_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

set -euo pipefail

# ============================================================
# 1. LOAD SOFTWARE MODULES
# ============================================================

module purge
module load gcc-6.3.0
module load angsd

echo "Loaded modules:"
module list

# ANGSD installation paths
ANGSD_DIR="/software/ngsTools/1.0.2/ngsTools/ANGSD/angsd"
ANGSD_MISC="${ANGSD_DIR}/misc"

REAL_SFS="${ANGSD_MISC}/realSFS"
THETA_STAT="${ANGSD_MISC}/thetaStat"

# Add misc directory to PATH
export PATH="${ANGSD_MISC}:${PATH}"

# Verify executables
if [[ ! -x "$REAL_SFS" ]]; then
    echo "ERROR: realSFS not found"
    exit 1
fi

if [[ ! -x "$THETA_STAT" ]]; then
    echo "ERROR: thetaStat not found"
    exit 1
fi

echo "realSFS: $REAL_SFS"
echo "thetaStat: $THETA_STAT"

# ============================================================
# 2. DEFINE DIRECTORIES
# ============================================================

BASE_DIR="/home/bistbs/downsampled_10X"
OUTPUT_DIR="${BASE_DIR}/ANGSD"

cd "$OUTPUT_DIR"

THREADS="${SLURM_CPUS_PER_TASK:-15}"

# ============================================================
# 3. DEFINE FIVE SAMPLES
# ============================================================

SAMPLES=(
    "SRR17129394_10X"
    "SRR17134085_10X"
    "SRR17134086_10X"
    "SRR17134087_10X"
    "SRR17134088_10X"
)

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

if (( TASK_ID < 0 || TASK_ID >= ${#SAMPLES[@]} )); then
    echo "ERROR: Invalid array task ID: $TASK_ID"
    exit 1
fi

SAMPLE="${SAMPLES[$TASK_ID]}"

SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE}"
PREFIX="${SAMPLE_DIR}/${SAMPLE}"

echo "================================================"
echo "SAMPLE: ${SAMPLE}"
echo "ARRAY TASK: ${TASK_ID}"
echo "THREADS: ${THREADS}"
echo "STARTED: $(date)"
echo "================================================"

mkdir -p "$SAMPLE_DIR"

# ============================================================
# 4. CHECK EXISTING SAF FILES
# ============================================================

SAF_IDX="${PREFIX}.saf.idx"
SAF_GZ="${PREFIX}.saf.gz"
SAF_POS="${PREFIX}.saf.pos.gz"

for FILE in "$SAF_IDX" "$SAF_GZ" "$SAF_POS"; do

    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing or empty SAF file: $FILE"
        exit 1
    fi

done

echo "All SAF input files are present."

# ============================================================
# 5. ESTIMATE FOLDED SITE FREQUENCY SPECTRUM
# ============================================================

SFS="${PREFIX}.folded.sfs"
SFS_TMP="${SFS}.tmp"

echo "================================================"
echo "STEP 1: FOLDED SFS"
echo "================================================"

if [[ -s "$SFS" ]]; then

    echo "Existing folded SFS found. Skipping estimation."

else

    rm -f "$SFS_TMP"

    "$REAL_SFS" \
        "$SAF_IDX" \
        -cores "$THREADS" \
        -fold 1 \
        > "$SFS_TMP"

    if [[ ! -s "$SFS_TMP" ]]; then
        echo "ERROR: Folded SFS estimation failed."
        exit 1
    fi

    mv "$SFS_TMP" "$SFS"

    echo "Folded SFS completed."

fi

# ============================================================
# 6. ESTIMATE THETA VALUES
# ============================================================

THETAS_GZ="${PREFIX}.thetas.gz"
THETAS_IDX="${PREFIX}.thetas.idx"

echo "================================================"
echo "STEP 2: SAF TO THETA"
echo "================================================"

if [[ -s "$THETAS_GZ" && -s "$THETAS_IDX" ]]; then

    echo "Existing theta files found. Skipping calculation."

else

    # Remove incomplete outputs before restarting
    rm -f "$THETAS_GZ" "$THETAS_IDX"

    "$REAL_SFS" saf2theta \
        "$SAF_IDX" \
        -sfs "$SFS" \
        -fold 1 \
        -outname "$PREFIX"

    if [[ ! -s "$THETAS_GZ" || ! -s "$THETAS_IDX" ]]; then
        echo "ERROR: Theta estimation failed."
        exit 1
    fi

    echo "Theta estimation completed."

fi

# ============================================================
# 7. GENOME-WIDE THETA STATISTICS
# ============================================================

THETA_STATS="${THETAS_IDX}.pestPG"

echo "================================================"
echo "STEP 3: GENOME-WIDE THETA STATISTICS"
echo "================================================"

if [[ -s "$THETA_STATS" ]]; then

    echo "Existing genome-wide statistics found. Skipping."

else

    "$THETA_STAT" do_stat "$THETAS_IDX"

    if [[ ! -s "$THETA_STATS" ]]; then
        echo "ERROR: Genome-wide theta statistics missing."
        exit 1
    fi

    echo "Genome-wide theta statistics completed."

fi

# ============================================================
# 8. SLIDING WINDOW THETA ANALYSIS
# ============================================================

WINDOW_PREFIX="${PREFIX}.50kb_10kb"
WINDOW_STATS="${WINDOW_PREFIX}.pestPG"

echo "================================================"
echo "STEP 4: SLIDING WINDOW ANALYSIS"
echo "Window size: 50000 bp"
echo "Step size: 10000 bp"
echo "================================================"

if [[ -s "$WINDOW_STATS" ]]; then

    echo "Existing sliding-window statistics found. Skipping."

else

    "$THETA_STAT" do_stat \
        "$THETAS_IDX" \
        -win 50000 \
        -step 10000 \
        -outnames "$WINDOW_PREFIX"

    if [[ ! -s "$WINDOW_STATS" ]]; then
        echo "ERROR: Sliding-window output missing."
        exit 1
    fi

    echo "Sliding-window analysis completed."

fi

# ============================================================
# 9. FINAL REPORT
# ============================================================

echo "================================================"
echo "ANALYSIS COMPLETED"
echo "SAMPLE: ${SAMPLE}"
echo "FINISHED: $(date)"
echo "================================================"

echo "Generated outputs:"

ls -lh \
    "$SFS" \
    "$THETAS_GZ" \
    "$THETAS_IDX" \
    "$THETA_STATS" \
    "$WINDOW_STATS"

echo "ALL DONE: ${SAMPLE}"
```
