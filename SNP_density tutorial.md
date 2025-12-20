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
- The reason I am keeping the window of 50kB is cause the variant size for my samples is too small in large windows.
- Therefore, I tried with the bare minimum and was able to get it with the code given below.

```bash
# --- 1. CONFIGURATION ---
VCF="Dama_gazelle_biallelic_snps_autosomes.vcf"
WINDOW=50000     
SAMPLES=("SRR17129394" "SRR17134085" "SRR17134086" "SRR17134087" "SRR17134088")

echo "Starting Heterozygosity Density Analysis..."
echo "Input VCF: $VCF"
echo "Window Size: $WINDOW bp"

# --- 2. INDIVIDUAL PROCESSING LOOP ---
for sample in "${SAMPLES[@]}"
do
    echo "------------------------------------------"
    echo "Processing $sample ..."
    
    # Extract only Heterozygous sites for the specific sample
    bcftools view -s "$sample" -g het -c1 "$VCF" -Oz -o "${sample}_HET.vcf.gz"
    
    # Index the new VCF (required for many downstream tools)
    bcftools index "${sample}_HET.vcf.gz"

    # Calculate SNP density in windows
    vcftools --gzvcf "${sample}_HET.vcf.gz" --SNPdensity "$WINDOW" --out "${sample}_density"
done

# --- 3. FINAL VERIFICATION & SUMMARY ---
echo ""
echo "--- FINAL VERIFICATION: HETEROZYGOUS SITE COUNTS ---"
printf "%-15s | %-15s\n" "Sample" "Total Het Sites"
echo "------------------------------------------"

for s in "${SAMPLES[@]}"
do
    # Sum the 3rd column of the .snpden file (SNP count per window)
    COUNT=$(awk 'NR>1 {sum+=$3} END {print sum}' "${s}_density.snpden")
    printf "%-15s | %-15s\n" "$s" "$COUNT"
done

# --- 4. DATA AGGREGATION (Optional but recommended) ---
# Creates a single file with an added column for the Sample ID
echo "Collating results into 'all_samples_density.txt'..."
echo -e "CHROM\tBIN_START\tSNP_COUNT\tVAR_PRI_MI\tSAMPLE" > all_samples_density.txt
for s in "${SAMPLES[@]}"
do
    awk -v sam="$s" 'NR>1 {print $0 "\t" sam}' "${s}_density.snpden" >> all_samples_density.txt
done

echo "Done."

```

##### Visualization of SNP density plot using painted chromosomes.

```bash
library(ggplot2)
library(dplyr)
library(readr)

# Input file (CSV)
INPUT_FILE <- "F:/Collaborative_Projects/Dama_Gazelle_Project/Heterozygosity/Heterozygosity.csv"

# Read CSV
data <- read_csv(INPUT_FILE, show_col_types = FALSE)

# Explicit species mapping
data <- data %>%
  mutate(
    Species = case_when(
      Sample %in% c("SRR17129394", "SRR17134085", "SRR17134086") ~ "Addra gazelle",
      Sample %in% c("SRR17134087", "SRR17134088") ~ "Mohrr gazelle",
      TRUE ~ "Unknown"
    )
  )

# Sort by heterozygosity
data <- data %>% arrange(Heterozygosity)

# Preserve sorted x-axis order
data$Sample <- factor(data$Sample, levels = data$Sample)

# Plot with species color labels
plot <- ggplot(data, aes(x = Sample, y = Heterozygosity, color = Species)) +
  geom_point(size = 3) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Sample",
    y = "Genome-wide heterozygosity",
    title = "Individual Heterozygosity Across Samples"
  ) +
  scale_color_manual(
    values = c(
      "Addra gazelle" = "skyblue",
      "Mohrr gazelle" = "orange",
      "Unknown" = "grey"
    )
  )

# Save plot
ggsave("heterozygosity_plot_labeled.jpeg", plot, width = 12, height = 6, dpi = 300)
# Save plot as PDF
ggsave("heterozygosity_plot_labeled.pdf", plot, width = 12, height = 6)


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
