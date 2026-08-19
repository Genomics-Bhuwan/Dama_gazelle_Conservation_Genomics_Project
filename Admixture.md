#### Run PCA for the Dama gazelle dataset

```bash

# 1. Load the modules
module purge
module load biocontainers/default
module load plink/1.90b6.21

# 2. Step 1: Run LD Pruning
plink --vcf dama_gazelle.SNPs.final.biallelic.vcf.gz --indep-pairwise 50 5 0.2 --out Dama_gazelle_prune --const-fid 0 --allow-extra-chr --threads 4

# 3. Step 2: Extract Pruned SNPs & Create PLINK Binary Bed Files
plink --vcf dama_gazelle.SNPs.final.biallelic.vcf.gz --extract Dama_gazelle_prune.prune.in --make-bed --out Dama_gazelle_LDpruned --const-fid 0 --allow-extra-chr --threads 4

# 4. Step 3: Run PCA
plink --bfile Dama_gazelle_LDpruned --pca var-wts --out Dama_gazelle_PCA --const-fid 0 --allow-extra-chr --threads 4
```
