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
