```bash
PROJECT=/home/bistbs/Dama_gazell_SNPeff
SNPEFF=${PROJECT}/snpEff
GENOME_ID="Dama_gazelle"
GENOME_NAME="Dama_gazelle"

GENOME_FASTA="${PROJECT}/Dama_gazelle_hifiasm-ULONT_primary.fasta"
GENOME_GFF="${PROJECT}/Addra_complete.genomic.gff"

# 1. Create directory structure
mkdir -p ${SNPEFF}/data/${GENOME_ID}

# 2. Stage genome and annotation files
cp ${GENOME_GFF} ${SNPEFF}/data/${GENOME_ID}/genes.gff
cp ${GENOME_FASTA} ${SNPEFF}/data/${GENOME_ID}/sequences.fa

# 3. Add genome configuration
if ! grep -q "^${GENOME_ID}\.genome" ${SNPEFF}/snpEff.config; then
    echo "${GENOME_ID}.genome : ${GENOME_NAME}" >> ${SNPEFF}/snpEff.config
fi

# 4. Build database
cd ${SNPEFF}
java -Xmx70g -jar snpEff.jar build \
    -gff3 \
    -v \
    -noCheckProtein \
    -noCheckCds \
    ${GENOME_ID}
```



```bash
#!/bin/bash
#SBATCH --job-name=DamaGazelle_SnpEff
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=24:00:00

set -euo pipefail

# Ensure local Java 21 binary is loaded
export PATH=$HOME/Dama_gazell_SNPeff/jdk-21.0.2+13/bin:$PATH

PROJECT=/home/bistbs/Dama_gazell_SNPeff
SNPEFF=${PROJECT}/snpEff

GENOME_ID="Dama_gazelle"
INPUT_VCF="${PROJECT}/dama_gazelle.SNPs.final.biallelic.minQ30.vcf.gz"
OUTPUT_VCF="${PROJECT}/dama_gazelle.annotated.vcf.gz"
REPORT_HTML="${PROJECT}/snpeff_report.html"

# Run Variant Annotation
echo "Annotating variants..."
cd ${SNPEFF}

java -Xmx70g -jar snpEff.jar ann \
    -v \
    -stats ${REPORT_HTML} \
    ${GENOME_ID} \
    ${INPUT_VCF} \
    | bgzip > ${OUTPUT_VCF}

# Index output VCF
echo "Indexing output VCF..."
tabix -p vcf ${OUTPUT_VCF}

echo "Done! Output generated at: ${OUTPUT_VCF}"
```
\
#####Find out the number of variants
set -euo pipefail

# Define paths for Dama gazelle execution
WORKING_DIR="/home/bistbs/Dama_gazell_SNPeff"
SNPSIFT_JAR="${WORKING_DIR}/snpEff/SnpSift.jar"
JAVA_EXEC="${WORKING_DIR}/jdk-21.0.2+13/bin/java"

INPUT_VCF="${WORKING_DIR}/dama_gazelle.annotated.vcf"
OUTPUT_DIR="${WORKING_DIR}/Variants_by_Consequences"

# 1. Create the output directory if it doesn't exist
echo "Creating output directory..."
mkdir -p ${OUTPUT_DIR}

# 2. Filter for Missense Variants
echo "Extracting Missense variants..."
${JAVA_EXEC} -jar ${SNPSIFT_JAR} filter "ANN[*].EFFECT has 'missense_variant'" \
  ${INPUT_VCF} \
  > ${OUTPUT_DIR}/Dama_gazelle_missense_sites.vcf

# 3. Filter for Synonymous Variants
echo "Extracting Synonymous variants..."
${JAVA_EXEC} -jar ${SNPSIFT_JAR} filter "ANN[*].EFFECT has 'synonymous_variant'" \
  ${INPUT_VCF} \
  > ${OUTPUT_DIR}/Dama_gazelle_synonymous_sites.vcf

# 4. Filter for Loss of Function (LoF) & Splicing
echo "Extracting Loss of Function (LoF) & Splicing variants..."
${JAVA_EXEC} -jar ${SNPSIFT_JAR} filter "(ANN[*].EFFECT has 'transcript_ablation') | (ANN[*].EFFECT has 'splice_donor_variant') | (ANN[*].EFFECT has 'splice_acceptor_variant') | (ANN[*].EFFECT has 'stop_gained') | (ANN[*].EFFECT has 'frameshift_variant') | (ANN[*].EFFECT has 'inframe_insertion') | (ANN[*].EFFECT has 'inframe_deletion') | (ANN[*].EFFECT has 'splice_region_variant')" \
  ${INPUT_VCF} \
  > ${OUTPUT_DIR}/Dama_gazelle_lof_sites.vcf

# 5. Filter for Intergenic Variants
echo "Extracting Intergenic variants..."
${JAVA_EXEC} -jar ${SNPSIFT_JAR} filter "ANN[*].EFFECT has 'intergenic_region'" \
  ${INPUT_VCF} \
  > ${OUTPUT_DIR}/Dama_gazelle_intergenic_sites.vcf

# -----------------------------------------------------------------
# Verification Block
# -----------------------------------------------------------------
echo "---------------------------------------------------"
echo "Filtering Complete! Total variant counts extracted:"
echo "---------------------------------------------------"
echo -n "Missense Sites:   " && grep -v "^#" ${OUTPUT_DIR}/Dama_gazelle_missense_sites.vcf | wc -l
echo -n "Synonymous Sites: " && grep -v "^#" ${OUTPUT_DIR}/Dama_gazelle_synonymous_sites.vcf | wc -l
echo -n "LoF/Splicing:     " && grep -v "^#" ${OUTPUT_DIR}/Dama_gazelle_lof_sites.vcf | wc -l
echo -n "Intergenic Sites: " && grep -v "^#" ${OUTPUT_DIR}/Dama_gazelle_intergenic_sites.vcf | wc -l
echo "---------------------------------------------------"
```
