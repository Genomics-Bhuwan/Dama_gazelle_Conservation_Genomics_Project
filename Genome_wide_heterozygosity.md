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

# --------------------------------------------------------
# 6B. ESTIMATE FOLDED SITE FREQUENCY SPECTRUM
# --------------------------------------------------------
echo "Step 2: Estimating folded SFS for ${SAMPLE}..."

realSFS \
    "${PREFIX}.saf.idx" \
    -P "${THREADS}" \
    -fold 1 \
    > "${PREFIX}.folded.sfs"

if [[ ! -s "${PREFIX}.folded.sfs" ]]; then
    echo "ERROR: SFS missing for $SAMPLE"
    exit 1
fi

echo "Folded SFS completed."

# --------------------------------------------------------
# 6C. ESTIMATE THETA STATISTICS
# --------------------------------------------------------
echo "Step 3: Estimating theta statistics for ${SAMPLE}..."

realSFS saf2theta \
    "${PREFIX}.saf.idx" \
    -sfs "${PREFIX}.folded.sfs" \
    -fold 1 \
    -outname "${PREFIX}"

if [[ ! -s "${PREFIX}.thetas.idx" ]]; then
    echo "ERROR: Theta index missing for $SAMPLE"
    exit 1
fi

echo "Theta estimation completed."

# --------------------------------------------------------
# 6D. GENOME-WIDE THETA STATISTICS
# --------------------------------------------------------
echo "Step 4: Genome-wide theta statistics for ${SAMPLE}..."

thetaStat do_stat \
    "${PREFIX}.thetas.idx"

# --------------------------------------------------------
# 6E. SLIDING WINDOW ANALYSIS
# --------------------------------------------------------
echo "Step 5: Sliding window analysis for ${SAMPLE}..."

thetaStat do_stat \
    "${PREFIX}.thetas.idx" \
    -win 50000 \
    -step 10000 \
    -outnames "${PREFIX}.50kb_10kb.theta"

# --------------------------------------------------------
# 6F. EXTRACT PER-SITE THETA VALUES
# --------------------------------------------------------
echo "Step 6: Extracting per-site theta values for ${SAMPLE}..."

thetaStat print \
    "${PREFIX}.thetas.idx" \
    2>/dev/null \
    > "${PREFIX}.per_site_theta.txt"

# ============================================================
# 7. FINAL TASK REPORT
# ============================================================

echo "================================================"
echo "ANGSD ANALYSIS COMPLETED FOR SAMPLE: ${SAMPLE}"
echo "Finished: $(date)"
echo "================================================"
```
