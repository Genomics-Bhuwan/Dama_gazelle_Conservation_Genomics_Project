##### SNP density tutorial
- SNP density is a tool used for visualizing the genetic diversity across the genome.
- It represents the number of SNPs in a heatmap.
- It is used for comparing genetic diversity between different populations or closely related species.
- To calculate the SNP density, count the number of SNPs over a certain length of the genome, typically ranging from 100kb to 1 megabase.
- Use VCFtools with the snpden function(Danecek et al., 2011).
- Vcf as input file and a table with the snp count for reach region interval.


##### Code
- Get only autosomes. Remove unaligned scaffolds and cytotypes.
- Input .vcf

```bash
/scratch/bistbs/Population_Genomic_Analysis/Heterozygosity/SNPden/Dama_gazelle_biallelic_snps_autosomes.vcf

# Samples in your VCF
SAMPLES=("SRR17129394" "SRR17134085" "SRR17134086" "SRR17134087" "SRR17134088")

# Parameters
MAF_THRESHOLD=0.1          # Minimum allele frequency
WINDOW_SIZE=1000000        # SNP density window size in bp (1 Mb)

# Create working directories
mkdir -p snpden/plots

# Step 1: Extract heterozygous sites for all samples
echo "Extracting heterozygous sites (MAF >= $MAF_THRESHOLD)..."
vcftools --vcf $VCF --recode --out snpden/AllSamples_hetsites --maf $MAF_THRESHOLD

# Step 2: Calculate SNP density for each individual
echo "Calculating SNP density per sample..."
for sample in "${SAMPLES[@]}"; do
    echo "Processing sample: $sample"

    # Extract sample-specific SNPs from heterozygous VCF
    echo $sample > snpden/${sample}_keep.txt
    vcftools --vcf snpden/AllSamples_hetsites.recode.vcf \
             --keep snpden/${sample}_keep.txt \
             --recode --out snpden/${sample}_hetsites

    # Calculate SNP density in windows
    vcftools --vcf snpden/${sample}_hetsites.recode.vcf \
             --SNPdensity $WINDOW_SIZE \
             --out snpden/${sample}_hetsites

    # Add sample name to the output file
    awk -v sample="$sample" 'NR==1{print $0"\tIndiv"} NR>1{print $0"\t"sample}' \
        snpden/${sample}_hetsites.snpden > snpden/${sample}_hetsites_id.snpden
done

# Step 3: Combine all sample SNP density files into one
echo "Combining all samples into a single SNP density file..."
head -n 1 snpden/${SAMPLES[0]}_hetsites_id.snpden > snpden/Dama_gazelle_hetsites.snpden
tail -q -n +2 snpden/*_id.snpden >> snpden/Dama_gazelle_hetsites.snpden

echo "✅ SNP density pipeline complete. Combined file:"
echo "   snpden/Dama_gazelle_hetsites.snpden"
```

##### Visualization of SNP density plot using painted chromosomes.

```bash

# -----------------------------
# 1️⃣ Set working directory
# -----------------------------
setwd("F:/Collaborative_Projects/Dama_Gazelle_Project/Heterozygosity/SNPden")  

# -----------------------------
# 2️⃣ Load packages
# -----------------------------
library(tidyverse)
library(gdata)

# -----------------------------
# 3️⃣ Read input files
# -----------------------------
snpden <- read.table("Dama_gazelle_hetsites.snpden", header = FALSE)
names_vec <- as.character(read.table("Dama_gazelle_IDs.txt")$V1)

# -----------------------------
# 4️⃣ Prepare data
# -----------------------------
colnames(snpden) <- c("CHROM", "BIN_START", "SNP_COUNT", "VARIANTS.KB", "Indiv")
snpden.master <- snpden

# Keep only chromosomes 1-17
valid_chroms <- as.character(1:17)
snpden.master <- snpden.master[snpden.master$CHROM %in% valid_chroms, ]
snpden.master$CHROM <- factor(snpden.master$CHROM, levels = valid_chroms)

# Convert to numeric
snpden.master$BIN_START <- as.numeric(as.character(snpden.master$BIN_START))
snpden.master$VARIANTS.KB <- as.numeric(as.character(snpden.master$VARIANTS.KB))

# -----------------------------
# 5️⃣ Assign species
# -----------------------------
snpden.master$Species <- case_when(
  snpden.master$Indiv %in% c("SRR17129394", "SRR17134085", "SRR17134086") ~ "Addra",
  snpden.master$Indiv %in% c("SRR17134087", "SRR17134088") ~ "Mhorr",
  TRUE ~ "Unknown"
)

# -----------------------------
# 6️⃣ Bin SNP densities (7 discrete bins)
# -----------------------------
snpden.master$groups <- cut(
  snpden.master$VARIANTS.KB,
  breaks = c(0, 0.1, 0.25, 0.5, 1, 2, 3, Inf),
  include.lowest = TRUE,
  labels = c("0-0.1","0.1-0.25","0.25-0.5","0.5-1","1-2","2-3","3+")
)

# Factorize individuals
snpden.master$Indiv <- factor(snpden.master$Indiv, levels = names_vec)

# -----------------------------
# 7️⃣ Create plots directory
# -----------------------------
if(!dir.exists("plots")) dir.create("plots")

# -----------------------------
# 8️⃣ Okabe-Ito color palette (7 colors)
# -----------------------------
okabe_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                  "#0072B2", "#D55E00", "#CC79A7")

# -----------------------------
# 9️⃣ Plot per species with bounding box around the whole plot
# -----------------------------
for(spec in unique(snpden.master$Species)) {
  
  snpden.spec <- subset(snpden.master, Species == spec)
  
  plot_title <- paste0("SNP Density - ", spec, " Gazelle")
  
  snpden_plot <- ggplot(snpden.spec, aes(x = BIN_START, y = Indiv)) +
    geom_tile(aes(fill = groups)) +
    facet_grid(CHROM ~ ., switch = 'y') +
    labs(x = 'Chromosome Position (Mb)', y = '', title = plot_title) +
    scale_fill_manual(values = okabe_colors, name = "Variants/kb") +
    scale_x_continuous(
      labels = function(x) paste0(x/1e6, " Mb"),
      breaks = seq(0, 250000000, 50000000),
      expand = c(0, 0)
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 12, face = "bold", color = "black"),
      axis.title.x = element_text(size = 14, face = "bold", color = "black"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # bounding box
      strip.text.y.left = element_text(angle = 0, size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    )
  
  # -----------------------------
  # Save plots
  # -----------------------------
  ggsave(filename = paste0("plots/", spec, "_combined_okabe_boundingBox.jpeg"),
         plot = snpden_plot,
         dpi = 600,
         units = 'cm',
         width = 28,
         height = 18,
         bg = "white")
  
  ggsave(filename = paste0("plots/", spec, "_combined_okabe_boundingBox.pdf"),
         plot = snpden_plot,
         device = "pdf",
         units = 'cm',
         width = 28,
         height = 18)
}
```
