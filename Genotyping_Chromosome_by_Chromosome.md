####Reference genome assembly indexing and creating the .dict file.
```bash
#!/bin/bash
#SBATCH --job-name=Dama_Reference_Index
#SBATCH --output=/anvil/projects/x-bio260092/Dama_Gazelle/Dama_Gazelle_Variant_Analysis/logs/reference_index_%j.out
#SBATCH --error=/anvil/projects/x-bio260092/Dama_Gazelle/Dama_Gazelle_Variant_Analysis/logs/reference_index_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=32G
#SBATCH --time=05:00:00
#SBATCH --partition=shared
#SBATCH --account=bio260092
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

module purge
module load gcc/11.2.0
module load samtools/1.12
module load gatk/4.1.8.1

BASE="/anvil/projects/x-bio260092/Dama_Gazelle"

REF_GZ="${BASE}/Dama_gazelle_hifiasm-ULONT_primary.fasta.gz"

REF_DIR="${BASE}/Dama_Gazelle_Variant_Analysis/reference"
REF="${REF_DIR}/Dama_gazelle_hifiasm-ULONT_primary.fasta"

LOGDIR="${BASE}/Dama_Gazelle_Variant_Analysis/logs"

mkdir -p "${REF_DIR}"
mkdir -p "${LOGDIR}"

echo "============================================================"
echo "DAMA GAZELLE REFERENCE INDEXING"
echo "============================================================"
echo "Start: $(date)"
echo "Node: $(hostname)"
echo ""
echo "Compressed reference:"
echo "${REF_GZ}"
echo ""
echo "Uncompressed reference:"
echo "${REF}"
echo "============================================================"

# ------------------------------------------------------------
# 1. Check compressed reference
# ------------------------------------------------------------

if [[ ! -f "${REF_GZ}" ]]; then
    echo "ERROR: Compressed reference not found:"
    echo "${REF_GZ}"
    exit 1
fi

echo "Compressed reference found:"
ls -lh "${REF_GZ}"

# ------------------------------------------------------------
# 2. Check/create uncompressed FASTA
# ------------------------------------------------------------

if [[ ! -f "${REF}" ]]; then

    echo ""
    echo "Uncompressed FASTA does not exist."
    echo "Decompressing reference..."

    TMP_REF="${REF}.tmp"

    rm -f "${TMP_REF}"

    gunzip -c "${REF_GZ}" > "${TMP_REF}"

    mv "${TMP_REF}" "${REF}"

    echo "Reference decompression completed."

else

    echo ""
    echo "Uncompressed FASTA already exists:"
    ls -lh "${REF}"

fi

# ------------------------------------------------------------
# 3. Make sure FASTA is not empty/corrupt
# ------------------------------------------------------------

if [[ ! -s "${REF}" ]]; then
    echo "ERROR: Reference FASTA is empty."
    exit 1
fi

echo ""
echo "Reference FASTA:"
ls -lh "${REF}"

# ------------------------------------------------------------
# 4. Create FASTA index
# ------------------------------------------------------------

if [[ ! -s "${REF}.fai" ]]; then

    echo ""
    echo "Creating FASTA index (.fai)..."

    rm -f "${REF}.fai"

    samtools faidx "${REF}"

    echo "FASTA index created."

else

    echo ""
    echo "FASTA index already exists:"
    ls -lh "${REF}.fai"

fi

# ------------------------------------------------------------
# 5. Create GATK sequence dictionary
# ------------------------------------------------------------

DICT="${REF_DIR}/Dama_gazelle_hifiasm-ULONT_primary.dict"

# Remove previously failed zero-byte dictionary
if [[ -f "${DICT}" && ! -s "${DICT}" ]]; then
    echo ""
    echo "Removing zero-byte dictionary:"
    echo "${DICT}"
    rm -f "${DICT}"
fi

if [[ ! -s "${DICT}" ]]; then

    echo ""
    echo "Creating GATK sequence dictionary..."

    gatk --java-options "-Xmx16G" CreateSequenceDictionary \
        -R "${REF}" \
        -O "${DICT}"

    echo "GATK sequence dictionary created."

else

    echo ""
    echo "GATK sequence dictionary already exists:"
    ls -lh "${DICT}"

fi

# ------------------------------------------------------------
# 6. Verify reference
# ------------------------------------------------------------

echo ""
echo "============================================================"
echo "REFERENCE FILES"
echo "============================================================"

ls -lh \
    "${REF}" \
    "${REF}.fai" \
    "${DICT}"

echo ""
echo "============================================================"
echo "FIRST 20 CONTIGS"
echo "============================================================"

head -20 "${REF}.fai"

echo ""
echo "============================================================"
echo "REFERENCE INDEXING COMPLETED"
echo "============================================================"
echo "End: $(date)"
echo "============================================================"
```
################################
####Joint Genotyping using Chromosomes.
```bash
#!/bin/bash
#SBATCH --job-name=Dama_ChrJoint
#SBATCH --output=/anvil/projects/x-bio260092/Dama_Gazelle/Dama_Gazelle_Variant_Analysis/logs/chr_joint_%A_%a.out
#SBATCH --error=/anvil/projects/x-bio260092/Dama_Gazelle/Dama_Gazelle_Variant_Analysis/logs/chr_joint_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=70G
#SBATCH --time=68:00:00
#SBATCH --partition=shared
#SBATCH --account=bio260092
#SBATCH --array=1-17
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

module purge
module load gcc/11.2.0
module load samtools/1.12
module load gatk/4.1.8.1

# ============================================================
# PATHS
# ============================================================

BASE="/anvil/projects/x-bio260092"

ANALYSIS="${BASE}/Dama_Gazelle/Dama_Gazelle_Variant_Analysis"

REF="${ANALYSIS}/reference/Dama_gazelle_hifiasm-ULONT_primary.fasta"

GVCF_DIR="${BASE}"

GENOMICSDB_DIR="${ANALYSIS}/01_GenomicsDB"

JOINT_DIR="${ANALYSIS}/02_Joint_Genotyping"

LOG_DIR="${ANALYSIS}/logs"

mkdir -p "${GENOMICSDB_DIR}"
mkdir -p "${JOINT_DIR}"
mkdir -p "${LOG_DIR}"

# ============================================================
# CHROMOSOME
# ============================================================

CHR="${SLURM_ARRAY_TASK_ID}"

echo "============================================================"
echo "DAMA GAZELLE CHROMOSOME JOINT GENOTYPING"
echo "============================================================"
echo "Chromosome:      ${CHR}"
echo "SLURM Job ID:    ${SLURM_JOB_ID}"
echo "Array Task ID:   ${SLURM_ARRAY_TASK_ID}"
echo "Node:            $(hostname)"
echo "Start:           $(date)"
echo "============================================================"

# ============================================================
# CHECK REFERENCE
# ============================================================

echo "Checking reference..."

if [[ ! -s "${REF}" ]]; then
    echo "ERROR: Reference FASTA not found:"
    echo "${REF}"
    exit 1
fi

if [[ ! -s "${REF}.fai" ]]; then
    echo "ERROR: FASTA index not found:"
    echo "${REF}.fai"
    exit 1
fi

DICT="${REF%.fasta}.dict"

if [[ ! -s "${DICT}" ]]; then
    echo "ERROR: GATK sequence dictionary not found:"
    echo "${DICT}"
    exit 1
fi

echo "Reference:"
ls -lh "${REF}"
echo "FASTA index:"
ls -lh "${REF}.fai"
echo "Sequence dictionary:"
ls -lh "${DICT}"

# ============================================================
# CHECK CHROMOSOME
# ============================================================

if ! grep -q "^${CHR}" "${REF}.fai"; then
    echo "ERROR: Chromosome ${CHR} was not found in reference."
    exit 1
fi

echo "Chromosome ${CHR} confirmed in reference."

# ============================================================
# GVCFs
# ============================================================

GVCF1="${GVCF_DIR}/SRR17129394.g.vcf.gz"
GVCF2="${GVCF_DIR}/SRR17134085.g.vcf.gz"
GVCF3="${GVCF_DIR}/SRR17134086.g.vcf.gz"
GVCF4="${GVCF_DIR}/SRR17134087.g.vcf.gz"
GVCF5="${GVCF_DIR}/SRR17134088.g.vcf.gz"

echo ""
echo "Checking GVCFs..."

for GVCF in \
    "${GVCF1}" \
    "${GVCF2}" \
    "${GVCF3}" \
    "${GVCF4}" \
    "${GVCF5}"
do
    if [[ ! -s "${GVCF}" ]]; then
        echo "ERROR: GVCF not found:"
        echo "${GVCF}"
        exit 1
    fi

    if [[ ! -s "${GVCF}.tbi" ]]; then
        echo "ERROR: GVCF index not found:"
        echo "${GVCF}.tbi"
        exit 1
    fi

    echo "OK: $(basename "${GVCF}")"
done

# ============================================================
# SAMPLE MAP
# ============================================================

SAMPLE_MAP="${ANALYSIS}/dama_gazelle.sample_map"

cat > "${SAMPLE_MAP}" << EOF
SRR17129394	${GVCF1}
SRR17134085	${GVCF2}
SRR17134086	${GVCF3}
SRR17134087	${GVCF4}
SRR17134088	${GVCF5}
EOF

echo ""
echo "Sample map:"
cat "${SAMPLE_MAP}"

# ============================================================
# CHROMOSOME-SPECIFIC OUTPUTS
# ============================================================

DB="${GENOMICSDB_DIR}/chr${CHR}"

OUT="${JOINT_DIR}/chr${CHR}.joint_genotyped.vcf.gz"

# ============================================================
# GENOMICSDBIMPORT
# ============================================================

echo ""
echo "============================================================"
echo "Running GenomicsDBImport"
echo "Chromosome: ${CHR}"
echo "Database:   ${DB}"
echo "============================================================"

if [[ -d "${DB}" ]]; then
    echo "ERROR: GenomicsDB workspace already exists:"
    echo "${DB}"
    echo ""
    echo "Delete it manually if you want to rerun this chromosome."
    exit 1
fi

gatk --java-options "-Xmx60G" GenomicsDBImport \
    --genomicsdb-workspace-path "${DB}" \
    --sample-name-map "${SAMPLE_MAP}" \
    --reader-threads 8 \
    --batch-size 50 \
    -L "${CHR}"

echo ""
echo "GenomicsDBImport completed for chromosome ${CHR}."

# ============================================================
# GENOTYPE GVCFs
# ============================================================

echo ""
echo "============================================================"
echo "Running GenotypeGVCFs"
echo "Chromosome: ${CHR}"
echo "Output:     ${OUT}"
echo "============================================================"

gatk --java-options "-Xmx60G" GenotypeGVCFs \
    -R "${REF}" \
    -V "gendb://${DB}" \
    -L "${CHR}" \
    -O "${OUT}"

# ============================================================
# INDEX VCF
# ============================================================

echo ""
echo "Indexing chromosome ${CHR} VCF..."

gatk IndexFeatureFile \
    -I "${OUT}"

# ============================================================
# FINAL CHECK
# ============================================================

echo ""
echo "============================================================"
echo "CHROMOSOME ${CHR} COMPLETED"
echo "============================================================"

ls -lh "${OUT}"
ls -lh "${OUT}.tbi"

echo ""
echo "End: $(date)"
echo "============================================================"
```
