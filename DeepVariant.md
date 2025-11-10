# Variant Calling using DeepVariant

##### DeepVariant is a deep learning-based variant caller that takes aligned reads (in BAM or CRAM format), produces pileup image tensors from them, classifies each tensor using a convolutional neural network, and finally reports the results in a standard VCF or gVCF file.

#### My rmdup_.bam file holds five samples in it but DeepVariant runs for one bam for each sample at a time. Therefore, I am firstly seperating all the five samples and making 5 .bam files. Then will be running DeepVariant.

```bash
#!/bin/bash

# ----------------Modules / Paths----------------- #
APPTAINER=/usr/bin/apptainer
DEEPVARIANT_IMAGE=docker://google/deepvariant:1.9.0

# ----------------User Variables------------------ #
INPUT_DIR=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates
OUTPUT_DIR=/localscratch/bistbs/DeepVariant
REF_GENOME=Dama_gazelle_hifiasm-ULONT_primary.fasta
BAM_FILE=all_samples_merged_rmdup.bam

# ----------------Create output directories------------------ #
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

