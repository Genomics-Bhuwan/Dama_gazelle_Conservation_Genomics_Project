# Variant Calling using DeepVariant

##### DeepVariant is a deep learning-based variant caller that takes aligned reads (in BAM or CRAM format), produces pileup image tensors from them, classifies each tensor using a convolutional neural network, and finally reports the results in a standard VCF or gVCF file.

```bash
#!/bin/sh
# ----------------Parameters---------------------- #
#$ -S /bin/sh
#$ -q sThC.q
#$ -l mres=10G,h_data=5G,h_vmem=5G
#$ -cwd
#$ -j y
#$ -N DeepVariant
#$ -o deepvariant_$JOB_ID.log

# ----------------Modules------------------------- #
module load docker

# ----------------Input/Output Paths---------------- #
INPUT_DIR=/localscratch/bistbs/DeepVariant/input
OUTPUT_DIR=/localscratch/bistbs/DeepVariant/output

BAM=${INPUT_DIR}/merged_rmdup.bam
REF=${INPUT_DIR}/Dama_gazelle_hifiasm-ULONT_primary.fasta
VCF_OUTPUT=${OUTPUT_DIR}/deepvariant_output.vcf.gz
GVCF_OUTPUT=${OUTPUT_DIR}/deepvariant_output.g.vcf.gz
LOG_DIR=${OUTPUT_DIR}/logs

mkdir -p ${OUTPUT_DIR} ${LOG_DIR}

# ----------------Run DeepVariant------------------ #
BIN_VERSION="1.9.0"

echo "🚀 Starting DeepVariant job on $HOSTNAME at $(date)"

docker run -v "${INPUT_DIR}":"/input" \
           -v "${OUTPUT_DIR}":"/output" \
           google/deepvariant:${BIN_VERSION} \
           /opt/deepvariant/bin/run_deepvariant \
           --model_type=WGS \
           --ref=/input/$(basename ${REF}) \
           --reads=/input/$(basename ${BAM}) \
           --output_vcf=/output/$(basename ${VCF_OUTPUT}) \
           --output_gvcf=/output/$(basename ${GVCF_OUTPUT}) \
           --num_shards=24 \
           --vcf_stats_report=true \
           --disable_small_model=true \
           --logging_dir=/output/logs \
           --dry_run=false

echo "✅ DeepVariant finished at $(date)"
```

