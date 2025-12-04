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
cd /scratch/bistbs/Population_Genomic_Analysis/Heterozygosity/SNPden/Sergei_SNPden

VCF=Dama_gazelle_biallelic_snps_autosomes.vcf

for sample in SRR17129394 SRR17134085 SRR17134086 SRR17134087 SRR17134088
do
    echo "Processing $sample ..."

    bcftools view \
        -s $sample \
        -g het \
        -v snps \
        $VCF \
        -Oz -o ${sample}_HET.vcf.gz

    bcftools index ${sample}_HET.vcf.gz
done

###Checking the coutns of the heterozygous SNPs for each sample.
for s in SRR17129394 SRR17134085 SRR17134086 SRR17134087 SRR17134088
do
    echo -n "$s: "
    bcftools view -H ${s}_HET.vcf.gz | wc -l
done

###### I want to have the bin size of 1 Mb.
SNP density window = 1 Mb
WINDOW=1000000

for sample in SRR17129394 SRR17134085 SRR17134086 SRR17134087 SRR17134088
do
    echo "Calculating SNP density for $sample ..."

    vcftools --gzvcf ${sample}_HET.vcf.gz \
             --SNPdensity $WINDOW \
             --out ${sample}_density

done

####Verify per-sample SNP density differs. Previously, all the SNP density plot were same. I want to make it different.
#### This checks how many windows each sample has:

for s in SRR17129394 SRR17134085 SRR17134086 SRR17134087 SRR17134088
do
    echo -n "$s: "
    awk 'NR>1 {print $4}' ${s}_density.snpden | sort -u | wc -l
done
### Combine all samples into one .snpden file.
echo -e "CHROM\tBIN_START\tBIN_END\tVARIANTS.KB\tIndiv" > combined_snpden.txt

for s in SRR17129394 SRR17134085 SRR17134086 SRR17134087 SRR17134088
do
    awk -v name=$s 'NR>1 {print $0"\t"name}' ${s}_density.snpden >> combined_snpden.txt
done

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
#####################################################################################
######################################################################################
#######################################################################################
# SNP density plot using MACE

#### Step 1. Get only heterozygous VCF only using bcftools
```bash
/scratch/bistbs/Population_Genomic_Analysis/Heterozygosity/SNPden/Sergei_SNPden/Dama_gazelle_biallelic_snps_autosomes.vcf

bcftools view -g het -v snps \
  Dama_gazelle_biallelic_snps_autosomes.vcf \
  -Oz -o Dama_gazelle_het_snps.vcf.gz

bcftools index Dama_gazelle_het_snps.vcf.gz
```
##### Step 2. Extract the Unique Scaffold IDs or the Chromosome number from the vcf file.
```bash
bcftools query -f '%CHROM\n' Dama_gazelle_het_snps.vcf.gz | sort | uniq > whitelist.txt
```

##### Step 3. Create the syn file. 
```bash
cd /scratch/bistbs/Population_Genomic_Analysis/Heterozygosity/SNPden/Sergei_SNPden
awk '{print $0 "\tChr"$0}' whitelist.txt > syn.txt
```
- This will produce syn.txt like this:
---
1	Chr1
2	Chr2
3	Chr3
4	Chr4
5	Chr5
6	Chr6
7	Chr7
8	Chr8
9	Chr9
10	Chr10
11	Chr11
12	Chr12
13	Chr13
14	Chr14
15	Chr15
16	Chr16
17	Chr17
---
##### Step 4. Ordering teh list
- This is a text file listing the scaffold names in the order you want them to be plotted.
- would be one per line.
```bash
cd /scratch/bistbs/Population_Genomic_Analysis/Heterozygosity/SNPden/Sergei_SNPden

# Create orderlist in numeric order using the syn file
sort -k2V syn.txt | awk '{print $2}' > orderlist.txt
```
##### Step 5. Create a scaffold length file. 
- It provides the length of each scaffold in the reference genome essential for plotting heterozygosity densities along the genome.
- Since, MACE plots the heterozygosity along the scaffolds.
- Without the actual lengths, MACE cannot scale the x-axis correctly for each scaffold.
- It wouldnot know where one scaffold ends and the next begins.
```bash
cd /scratch/bistbs/Population_Genomic_Analysis/Heterozygosity/SNPden/Sergei_SNPden

# Extract scaffold lengths
cut -f1,2 /scratch/bistbs/Population_Genomic_Analysis/Heterozygosity/SNPden/Dama_gazelle_hifiasm-ULONT_primary.fasta.fai > len.txt
```

##### Step 6. Running the MACE from its script.https://github.com/mahajrod/mace
- Donwload MACE into the respective folder where you want to plot the heterozygosity density plot.

```bash
cd /scratch/bistbs/Population_Genomic_Analysis/Heterozygosity/SNPden/Sergei_SNPden/MACE/scripts
```
-  Run the Python script
```bash
python draw_variant_window_densities.py \
  -i /scratch/bistbs/Population_Genomic_Analysis/Heterozygosity/SNPden/Sergei_SNPden/Dama_gazelle_het_snps.vcf.gz \
  -o /scratch/bistbs/Population_Genomic_Analysis/Heterozygosity/SNPden/Sergei_SNPden/Dama_gazelle_density \
  -n /scratch/bistbs/Population_Genomic_Analysis/Heterozygosity/SNPden/Sergei_SNPden/len.txt \
  --scaffold_white_list /scratch/bistbs/Population_Genomic_Analysis/Heterozygosity/SNPden/Sergei_SNPden/whitelist.txt \
  --scaffold_syn_file /scratch/bistbs/Population_Genomic_Analysis/Heterozygosity/SNPden/Sergei_SNPden/syn.txt \
  --scaffold_ordered_list /scratch/bistbs/Population_Genomic_Analysis/Heterozygosity/SNPden/Sergei_SNPden/orderlist.txt


```
- This will get you png or svg files as your output.
