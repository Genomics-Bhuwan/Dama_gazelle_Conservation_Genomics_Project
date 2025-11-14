#### Fst calculation

#### Filtration before final Fst calculation.
```bash
vcftools --vcf /scratch/bistbs/Population_Genomic_Analysis/Fst/Dama_gazelle_biallelic_snps_filtered.recode.vcf \
  --max-missing 0.9 \
  --minQ 30 \
  --minDP 5 \
  --recode --stdout \
  | bcftools view -Oz -o /scratch/bistbs/Population_Genomic_Analysis/Fst/filtered_for_fst.vcf.gz
```
---
-------------------------------------------------------------
Windowed FST (Weir & Cockerham) – Explanation
-------------------------------------------------------------
- FST can be calculated per SNP (per-site FST), but individual SNP
- estimates are extremely noisy, especially with small sample sizes.

##### To obtain smoother and more biologically interpretable estimates,
##### we calculate FST in sliding genomic windows.

######   --fst-window-size 50000
- This defines the size of each genomic window (here 50 kb).
- All SNPs that fall inside the 50 kb interval are used to
##### compute one averaged FST value for that window.

#####   --fst-window-step 10000
- After computing FST in the first 50 kb window, we move the
##### window forward by 10 kb (the “step size”) and compute the
##### next window. This creates overlapping sliding windows.

##### Why use windowed FST?
 • Reduces noise caused by single SNPs with unreliable frequencies.
 • Provides a smoother profile of differentiation across the genome.
 • More robust when sample sizes are small (like n=3 vs n=2).
 • Helps identify genomic regions of high differentiation.

##### Interpretation:
-  The output ".windowed.weir.fst" file contains one FST estimate
- per window, including:
#####       CHROM   BIN_START   BIN_END   N_VARIANTS   MEAN_FST

#####   MEAN_FST is the average genetic differentiation in that window.

##### Recommended:
  • Increase window size (e.g., 100 kb–500 kb) if sample sizes are very small.
  • Drop windows with "nan" (no informative SNPs).
  • Compute mean FST across all windows for a genome-wide estimate.
##### -------------------------------------------------------------
---

#### Fst calculation using Weir and Cockerham calculation.
```bash
vcftools --gzvcf /scratch/bistbs/Population_Genomic_Analysis/Fst/filtered_for_fst.vcf.gz \
  --weir-fst-pop /scratch/bistbs/Population_Genomic_Analysis/Fst/popA.txt \
  --weir-fst-pop /scratch/bistbs/Population_Genomic_Analysis/Fst/popB.txt \
  --fst-window-size 50000 \
  --fst-window-step 10000 \
  --out /scratch/bistbs/Population_Genomic_Analysis/Fst/Dama_Addra_vs_Mohrr_windowed_50kb_step10kb
```


```
