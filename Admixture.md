# 🧬 Dama gazelle — ADMIXTURE and PCA Workflow

**Author:** Bhuwan Singh Bist  
**Date:** 2025-11-12  
**Purpose:** Complete workflow for preparing and running ADMIXTURE as well as PCA analysis on *Dama gazelle* using a high-quality biallelic SNP VCF.

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
- minQ 30 --> keep sites that have at least a minimum quality score of 30. Meaning that the base quality accuracy at this postiton is 99.9% correct.
- remove-indels --> this removes indels and we are interested in only SNPs. 
- max-missing 0.8 --> keep positions where there is data for at least 80% of individuals .
  
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

#### Step 4: Convert VCF to PLINK Format: Plink (https://zzz.bwh.harvard.edu/plink/)

# ---------------------------------------------------------------
# PLINK pipeline for Dama gazelle population genomics
# Converts VCF to PLINK format, performs LD pruning, and prepares 
# files for downstream analyses like PCA or Admixture
# ---------------------------------------------------------------

# -------------------------------
# Step A: Load PLINK module
# -------------------------------
# download or install the plink. I had already downloaded and am not downloading now.
# -------------------------------
# Step B: LD pruning preparation
# Convert VCF to PLINK format and prune SNPs in LD
# -------------------------------
```bash
VCF=/scratch/bistbs/Population_Genomic_Analysis/PCA/Dama_gazelle_biallelic_snps_autosomes.vcf
OUT=Dama_gazelle_LDprune

plink --vcf $VCF \
      --indep-pairwise 50 5 0.5 \
      --out $OUT \
      --const-fid 0 \
      --allow-extra-chr

```
# -------------------------------
# Step C: Create LD-pruned binary PLINK dataset
# Only extract SNPs kept after LD pruning
# -------------------------------
```bash
plink --vcf /scratch/bistbs/Population_Genomic_Analysis/Admixture/Dama_gazelle_biallelic_snps_filtered.recode.vcf \
      --double-id \
      --allow-extra-chr \
      --set-missing-var-ids @:# \
      --extract Dama_gazelle_plink.prune.in \
      --make-bed \
      --out Dama_gazelle_LDpruned
```
# ---------------------------------------------------------------
# Notes:
# - 'Dama_gazelle_plink.prune.in' = SNPs kept after LD pruning
# - 'Dama_gazelle_LDpruned' = ready for PCA/Admixture analyses
# - Adjust window size / step size / r² threshold depending on species' LD decay
# ---------------------------------------------------------------

---
Output files:

Dama_gazelle_admixture.bed
Dama_gazelle_admixture.bim
Dama_gazelle_admixture.fam
---
# ---------------------------------------------------------------
# Step D. Replace the first column of the .bim file with "0" for ADMIXTURE
awk '{$1="0"; print $0}' Dama_gazelle_LDpruned.bim > Dama_gazelle_LDpruned.bim.tmp
mv Dama_gazelle_LDpruned.bim.tmp Dama_gazelle_LDpruned.bim
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Step E. ADMIXTURE: Individual Admixture Proportions with Admixture and visualization using MapMixture.
# Now we can infer individual admixture proportions with Admixture. This is a genetic clustering program to define populations and assign individuals to them.  
# We will use the plink files we just make to run Admixture. 


# ====== Step E.1: ADMIXTURE Analysis with Cross-Validation ======
```bash
# Base name of your PLINK binary files (no extension)
FILE="Dama_gazelle_LDpruned"

ADMIXTURE_PATH="/scratch/bistbs/Population_Genomic_Analysis/Admixture/Admixture_mapmixture/admixture_linux-1.3.0/admixture"

# Run ADMIXTURE for K = 1 to 10 with cross-validation
for K in {1..10}; do
    echo "Running ADMIXTURE for K=${K}..."
    ${ADMIXTURE_PATH} --cv ${FILE}.bed $K > log${K}.out
done

# Summarize CV errors to find the best K
grep "CV error" log*.out | awk '{print $3, $4, $5}' | sed -e 's/(//;s/)//;s/://;s/K=//' > ${FILE}_cv_errors.txt

# Follow the Steps from MapMixture Package for rest.

```



##### PCA
- run PCA with PLINK (var-wts gives weighted PCA like SmartPCA)
- add --allow-extra-chr because you have nonstandard chromosome names
```bash
VCF=/scratch/bistbs/Population_Genomic_Analysis/PCA/Dama_gazelle_biallelic_snps_autosomes.vcf
OUT=Dama_gazelle_LDprune

plink --vcf $VCF \
      --indep-pairwise 50 5 0.5 \
      --out $OUT \
      --const-fid 0 \
      --allow-extra-chr
```

VCF=/scratch/bistbs/Population_Genomic_Analysis/PCA/Dama_gazelle_biallelic_snps_autosomes.vcf
OUT=Dama_gazelle_LDprune
FINAL=Dama_gazelle_LDpruned

##### Extract LD-pruned SNPs and create PLINK binary files
```bash
plink --vcf $VCF \
      --extract ${OUT}.prune.in \
      --make-bed \
      --out $FINAL \
      --const-fid 0 \
      --allow-extra-chr
```
##### Get the final data for visualization

```bash
cd /scratch/bistbs/Population_Genomic_Analysis/PCA
plink --bfile Dama_gazelle_LDpruned \
      --pca var-wts \
      --const-fid 0 \
      --allow-extra-chr \
      --out Dama_gazelle_PCA

```

##### Plot the eigenvector to get the samples in Principal Components.
##### Plot the eigenvalues to get the screeplot to see the variance explained by each principal component.
```bash
# Install necessary packages (run once)
install.packages("ggplot2")
install.packages("ggrepel")

# Load libraries
library(ggplot2)
library(ggrepel)

# Set working directory
setwd("F:/Collaborative_Projects/Dama_Gazelle_Project/PCA")

# ---- READ DATA ----
dat <- read.csv("Dama_gazelle_PCA.csv", header = TRUE)

# ---- PCA SCATTER PLOT ----
p <- ggplot(dat, aes(x = PC1, y = PC2, color = Sample, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(show.legend = FALSE) +
  theme_classic() +
  labs(
    title = "PCA Plot of Dama Gazelle",
    x = "Principal Component 1",
    y = "Principal Component 2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# Display PCA scatter plot
print(p)

# ---- SAVE PCA PLOT ----
ggsave("Dama_gazelle_PCA_plot_300dpi.jpeg", p, width = 7, height = 5, dpi = 300)
ggsave("Dama_gazelle_PCA_plot_600dpi.jpeg", p, width = 7, height = 5, dpi = 600)
ggsave("Dama_gazelle_PCA_plot_300dpi.tiff", p, width = 7, height = 5, dpi = 300, compression = "lzw")
ggsave("Dama_gazelle_PCA_plot_600dpi.tiff", p, width = 7, height = 5, dpi = 600, compression = "lzw")
ggsave("Dama_gazelle_PCA_plot.pdf", p, width = 7, height = 5)

# ---- EIGENVALUE BARPLOT ----
Eigenvalues <- c(3.60779, 1.05618, 0.755254, 0.392603, -0.00613124)
bar_data <- data.frame(
  PC = paste0("PC", 1:length(Eigenvalues)),
  Eigenvalue = Eigenvalues
)

bp <- ggplot(bar_data, aes(x = PC, y = Eigenvalue)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_classic() +
  labs(
    title = "Scree Plot of Eigenvalues",
    x = "Principal Components",
    y = "Eigenvalues"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Display Scree (Eigenvalue) plot
print(bp)

# ---- SAVE EIGENVALUE PLOT ----
ggsave("Dama_gazelle_Eigenvalue_plot_300dpi.jpeg", bp, width = 7, height = 5, dpi = 300)
ggsave("Dama_gazelle_Eigenvalue_plot_600dpi.jpeg", bp, width = 7, height = 5, dpi = 600)
ggsave("Dama_gazelle_Eigenvalue_plot_300dpi.tiff", bp, width = 7, height = 5, dpi = 300, compression = "lzw")
ggsave("Dama_gazelle_Eigenvalue_plot_600dpi.tiff", bp, width = 7, height = 5, dpi = 600, compression = "lzw")
ggsave("Dama_gazelle_Eigenvalue_plot.pdf", bp, width = 7, height = 5)

cat("✅ All plots saved in JPEG, TIFF (300 & 600 dpi), and PDF formats.\n")
```


