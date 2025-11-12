# 🧬 Dama gazelle — ADMIXTURE Workflow

**Author:** Bhuwan Singh Bist  
**Date:** 2025-11-12  
**Purpose:** Complete workflow for preparing and running ADMIXTURE analysis on *Dama gazelle* using a high-quality biallelic SNP VCF.

---

##  Step 1: Directory Setup

```bash
mkdir -p /scratch/bistbs/Population_Genomic_Analysis/Admixture
cd /scratch/bistbs/Population_Genomic_Analysis/Admixture
```

## Step 2: Input Data
- Biallelic SNPs only + Indels removed + High quality (QUAL ≥ 30) + Missingness ≤ 7%
```bash
/scratch/bistbs/Population_Genomic_Analysis/Admixture/Dama_gazelle_biallelic_snps.vcf
```
## Step 3: Filter VCF with VCF-tools to work on maximum missingness.
```bash
module load vcf-tools

VCF=/scratch/bistbs/Population_Genomic_Analysis/Admixture/Dama_gazelle_biallelic_snps.vcf
OUT=/scratch/bistbs/Population_Genomic_Analysis/Admixture/Dama_gazelle_filtered

vcftools --vcf ${VCF} \
  --minQ 30 \
  --max-missing 0.8 \
  --remove-indels \
  --recode --recode-INFO-all \
  --out ${OUT}
```

#### Step 4: Convert VCF to PLINK Format
```bash
module load plink

plink --vcf /scratch/bistbs/Population_Genomic_Analysis/Admixture/Dama_gazelle_filtered.recode.vcf \
  --double-id \
  --allow-extra-chr \
  --make-bed \
  --out /scratch/bistbs/Population_Genomic_Analysis/Admixture/Dama_gazelle_admixture
```

Output files:

Dama_gazelle_admixture.bed
Dama_gazelle_admixture.bim
Dama_gazelle_admixture.fam

⚙️ Step 5: Run ADMIXTURE
module load admixture

cd /scratch/bistbs/Population_Genomic_Analysis/Admixture

for K in {2..6}; do
   admixture --cv Dama_gazelle_admixture.bed $K | tee log${K}.out
done


Outputs:

Dama_gazelle_admixture.${K}.Q – Individual ancestry proportions

Dama_gazelle_admixture.${K}.P – Allele frequencies per cluster

log${K}.out – Cross-validation logs

📊 Step 6: Determine Best K
grep -h "CV error" log*.out


Example output:

CV error (K=2): 0.541
CV error (K=3): 0.498
CV error (K=4): 0.503
CV error (K=5): 0.510
CV error (K=6): 0.517


👉 Best K = 3 (lowest CV error).

🎨 Step 7: Visualize ADMIXTURE Results in R
# Load packages
library(ggplot2)
library(tidyverse)

# Set working directory
setwd("/scratch/bistbs/Population_Genomic_Analysis/Admixture")

# Read ADMIXTURE output
qdata <- read.table("Dama_gazelle_admixture.3.Q")
fam <- read.table("Dama_gazelle_admixture.fam")

# Combine data
qdata$ID <- fam$V2
colnames(qdata) <- c("Cluster1", "Cluster2", "Cluster3", "ID")

# Plot ancestry proportions
ggplot(qdata, aes(x = ID, y = 1, fill = Cluster1)) +
  geom_bar(aes(y = Cluster1), stat = "identity") +
  geom_bar(aes(y = Cluster2), stat = "identity", position = "stack") +
  geom_bar(aes(y = Cluster3), stat = "identity", position = "stack") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 7)) +
  labs(title = "Dama gazelle ADMIXTURE (K=3)", x = "Individuals", y = "Ancestry Proportion")
