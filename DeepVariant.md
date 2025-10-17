# Variant Calling using DeepVariant

##### DeepVariant is a deep learning-based variant caller that takes aligned reads (in BAM or CRAM format), produces pileup image tensors from them, classifies each tensor using a convolutional neural network, and finally reports the results in a standard VCF or gVCF file.

```bash
#!/bin/bash
# ----------------Parameters---------------------- #
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=20G,h_data=10G,h_vmem=10G
#$ -cwd
#$ -j y
#$ -N DeepVariant
#$ -o deepvariant.log

# ----------------Modules / Paths----------------- #
APPTAINER=/usr/bin/apptainer
DEEPVARIANT_IMAGE=docker://google/deepvariant:1.9.0

# ----------------User Variables------------------ #
INPUT_DIR=/localscratch/bistbs
OUTPUT_DIR=/scratch/bistbs/DeepVariant_output
REF_GENOME=YOUR_REFERENCE.fasta        # <-- Replace with your reference genome
BAM_FILE=YOUR_SAMPLE.bam               # <-- Replace with your BAM file

# Make output directories if they don't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}/logs

# ----------------Run DeepVariant---------------- #
echo "Starting DeepVariant at $(date)"

$APPTAINER run \
  -B ${INPUT_DIR}:/input \
  -B ${OUTPUT_DIR}:/output \
  ${DEEPVARIANT_IMAGE} \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/input/${REF_GENOME} \
    --reads=/input/${BAM_FILE} \
    --output_vcf=/output/${BAM_FILE%.bam}.dv.vcf \
    --output_gvcf=/output/${BAM_FILE%.bam}.dv.g.vcf \
    --num_shards=$(nproc) \
    --vcf_stats_report=true \
    --disable_small_model=true \
    --logging_dir=/output/logs

echo "DeepVariant finished at $(date)"

```

