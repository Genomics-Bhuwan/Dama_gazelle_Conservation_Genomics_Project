#### Heterozygosity tutorial
---
- We will use two strategies to calcualte the average heterozygosity.
- The first step is using bcftools and vcftools.
- The second step will be exploring the heterozygosity variation across the genome using ANGSD(Analysis of Next Generation Sequencing Data).
- This calculates the heterozygosity of each sample in a vcf file.
- Firstly, it extracts the list of the sample names and calculates the number of heterozygous variants for each sample.
- The results are saved as tab-separated file with the name "heterozygosity.tsv" and two columns aka "Sample" and "Heterozygosity".
 --- 



##### Method A. Average heterozygosity

```bash
#Replace these with the actual file names
VCF_FILE="/scratch/bistbs/Population_Genomic_Analysis/Heterozygosity/Dama_gazelle_biallelic_snps_autosomes.vcf"
OUTPUT_FILE="Dama_heterozygosity_v5.tsv"
GENOME_LENGTH=3108406478

SAMPLES=$(bcftools query -l $VCF_FILE)

echo -e "Sample\tHeterozygous_sites\tHeterozygosity" > $OUTPUT_FILE

for SAMPLE in $SAMPLES; do
  HETEROZYGOUS=$(bcftools query -s $SAMPLE -f '[%GT\n]' $VCF_FILE | grep -F "0/1" | wc -l)
  HETEROZYGOSITY=$(echo "scale=7; $HETEROZYGOUS / $GENOME_LENGTH" | bc)
  echo -e "$SAMPLE\t$HETEROZYGOUS\t$HETEROZYGOSITY" >> $OUTPUT_FILE
done


```
##### Method A. Visualization for the Heterozygosity.

```r
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



### B. Genome-wide heterozygosity using Window-based approach.
- This approach allows us to visualize the variation throughtout the genome with the ability to focus on regions of particularly higher or lower diversity.
- It is also a comparision with other methods such as ROH.
- ANGSD peforms it. It is a software designed to analyze low-depth NGS data. It handles large-scale sequencing data particularly for non-model organisms.
---
Key capabilities of ANGSD include:
---
1. Genotype calling: ANGSD can estimate genotype likelihoods from sequencing data without actually calling genotypes, which helps reduce biases and errors introduced by hard genotype calls.
2. SNP discovery: The software can identify Single Nucleotide Polymorphisms (SNPs) and estimate their allele frequencies while accounting for sequencing errors and varying levels of sequencing depth.
3. Genetic association studies: ANGSD can perform Genome-Wide Association Studies (GWAS) and estimate genotype-phenotype associations using mixed linear models or generalized linear models.
4. Population structure: ANGSD can estimate population structure and admixture proportions using principal component analysis (PCA) or model-based approaches.
5. Selection scans: The software can detect signals of positive or balancing selection, using various statistics like FST, Tajima's D, and nucleotide diversity.
6. Handling of various file formats: ANGSD can work with various input file formats, including BAM, CRAM, and VCF files, and output results in multiple formats suitable for downstream analyses.

Overall, ANGSD is a versatile and powerful tool for analyzing NGS data, particularly for non-model organisms or low-coverage sequencing projects. It provides a comprehensive suite of analysis tools that cater to various research objectives in population genetics and genomics.
---
#### How to run heterozygosity using ANGSD

```bash
module load bioinformatics/angsd/0.921

angsd -P <threads> -i <input_bam_file> -anc <ancestral_fasta_file> -dosaf <dosaf_value> -gl <genotype_likelihood_method> -C <base_quality_adjustment> -minQ <min_base_quality> -minmapq <min_mapping_quality> -fold <fold_value> -out <output_file> -ref <reference_fasta_file> -r <region_of_interest>
```

- `P <threads>` - Sets the number of threads to be used in parallel.
- `i <input_bam_file>` - Specifies the input file as a BAM file.
- `anc <ancestral_fasta_file>` - Specifies the ancestral fasta reference file.
- `dosaf <dosaf_value>` - Computes the Site Frequency Spectrum (SFS) based on the genotype likelihoods.
- `gl <genotype_likelihood_method>` - Specifies the method used for calculating genotype likelihoods.
- `C <base_quality_adjustment>` - Adjusts the base quality score by a specified value before using it.
- `minQ <min_base_quality>` - Sets the minimum base quality score required.
- `minmapq <min_mapping_quality>` - Sets the minimum mapping quality score required.
- `fold <fold_value>` - Indicates whether you are analyzing folded SFS or unfolded SFS.
- `out <output_file>` - Specifies the output file path and name.
- `ref <reference_fasta_file>` - Specifies the reference fasta file.
- `r <region_of_interest>` - Specifies the region of interest.

```bash
realSFS -nsites <number_of_sites> <input_saf_idx_file> > <output_est_ml_file>
```

- `nsites <number_of_sites>` - Specifies the number of sites to be considered for the estimation. Replace `<number_of_sites>` with the desired number.
- `<input_saf_idx_file>` - Specifies the input .saf.idx file, which is the output from the ANGSD program that contains information about the SFS.
- `> <output_est_ml_file>` - Specifies the output file path and name for the maximum likelihood estimate of the SFS.

ANGSD does not differentiate the chromosomes if you run the whole genome at once, that is why we need to use the ‘region’ variable when running ANGSD to specify the chromosome/scaffold. This will be useful to look the heterozygosity throughout the genome. 
##### This is the code I used for running the ANGSD and SFS for Dama gazelle for the five samples as given in the code below.
```bash
#!/bin/bash -l
#SBATCH --job-name=ANGSD_Het
#SBATCH --time=300:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --partition=batch
#SBATCH --output=angsd_het_%A.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

#Set variables
# ANGSD executable (full path)
ANGSD_EXE="/scratch/bistbs/Population_Genomic_Analysis/PSMC/angsd/angsd"
REAL_SFS_EXE="/scratch/bistbs/Population_Genomic_Analysis/PSMC/angsd/misc/realSFS"  # adjust if different

# Directories and reference files
bam_dir="/scratch/bistbs/Population_Genomic_Analysis/PSMC"
output_dir="${bam_dir}/heterozygosity"
ancestral_fasta="${bam_dir}/Dama_gazelle_hifiasm-ULONT_primary.fasta"
reference_fasta="${bam_dir}/Dama_gazelle_hifiasm-ULONT_primary.fasta"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# List of BAM files
samples=("SRR17129394.bam" "SRR17134085.bam" "SRR17134086.bam" "SRR17134087.bam" "SRR17134088.bam")

# Loop over each sample
for bam_file in "${samples[@]}"; do
    SAMPLE=$(basename "$bam_file" .bam)
    
    # Loop over autosomes 1-17 (matching FASTA headers)
    for i in {1..17}; do
        CHR="$i"
        echo "[$(date)] Processing $SAMPLE chromosome $CHR..."
        
        # Run ANGSD to calculate SAF
        $ANGSD_EXE -P 24 \
            -i "${bam_dir}/${bam_file}" \
            -anc "$ancestral_fasta" \
            -ref "$reference_fasta" \
            -dosaf 1 \
            -gl 1 \
            -C 50 \
            -minQ 20 \
            -minmapq 30 \
            -out "${output_dir}/${SAMPLE}.${CHR}" \
            -r "$CHR"
        
        # Estimate folded SFS
        $REAL_SFS_EXE -fold 1 "${output_dir}/${SAMPLE}.${CHR}.saf.idx" > "${output_dir}/${SAMPLE}.${CHR}.est.ml"
    done
done

echo "[$(date)] All ANGSD runs completed."
```

- Replace the `<placeholders>` with the desired values for your specific analysis, and update the paths for input and output files as needed. The `$SAMPLE` variable should also be set to the appropriate sample name.

- This script will loop through scaffolds 1 to 17, running the ANGSD command and then the realSFS command for each scaffold. The results will be saved in the specified output directory.

- Now we can add the sample name and scaffold number for each line in our output. This will make our work easier when we want to plot our results. 

```bash
#!/bin/bash
#!/bin/bash

# ----------------------------
# Annotate realSFS output files with number of lines, sample name, and scaffold
# ----------------------------

output_dir="/scratch/bistbs/Population_Genomic_Analysis/PSMC/heterozygosity"

# List of BAM/sample names (without .bam)
samples=("SRR17129394" "SRR17134085" "SRR17134086" "SRR17134087" "SRR17134088")

# Loop over each sample
for SAMPLE in "${samples[@]}"; do
    echo "[$(date)] Annotating $SAMPLE files..."

    # Loop over scaffolds 1 to 17
    for i in {1..17}; do
        input_file="${output_dir}/${SAMPLE}.${i}.est.ml"
        output_file="${output_dir}/${SAMPLE}.${i}.est.ml.annotated"

        if [ -f "$input_file" ]; then
            # Count number of lines
            num_lines=$(wc -l < "$input_file")

            # Annotate each line with line count, sample, and scaffold
            awk -v lines="$num_lines" -v sample="$SAMPLE" -v scaffold="$i" '{print lines, sample, scaffold, $0}' "$input_file" > "$output_file"

            # Replace original file with annotated version
            mv "$output_file" "$input_file"

            echo "[$(date)] Annotated $input_file"
        else
            echo "[$(date)] WARNING: File $input_file not found, skipping."
        fi
    done
done

echo "[$(date)] All annotation completed!"

```

Concatenate files

```bash
#!/bin/bash

# ----------------------------
# Concatenate all realSFS output files for a given sample
# ----------------------------

input_directory="/scratch/bistbs/Population_Genomic_Analysis/PSMC/heterozygosity"
output_directory="/scratch/bistbs/Population_Genomic_Analysis/PSMC/heterozygosity/concatenated"
mkdir -p "$output_directory"

# Output file
output_file="${output_directory}/all_samples_est_ml_concatenated.txt"

# Remove output file if it already exists
[ -f "$output_file" ] && rm "$output_file"

# List of samples
samples=("SRR17129394" "SRR17134085" "SRR17134086" "SRR17134087" "SRR17134088")

# Loop over each sample and scaffold 1-17
for SAMPLE in "${samples[@]}"; do
    for i in {1..17}; do
        input_file="${input_directory}/${SAMPLE}.${i}.est.ml"
        if [ -f "$input_file" ]; then
            cat "$input_file" >> "$output_file"
        else
            echo "WARNING: $input_file not found, skipping."
        fi
    done
done

echo "All files concatenated into $output_file"

```

##### Plot the results for one chromosome using R

```bash
library(tidyverse)
library(viridis)
library(scales)

# Read the concatenated file
het_master <- read.table("/path/to/file/all_est_ml_concatenated.txt")

# Process and plot
het_master %>%
  rename(sample = V2,
         chromosome = V3) %>%
  mutate(
    heterozygosity = V5 / (V4 + V5),
    position = ((V1 * 200000) - 200000)  # rough genomic position
  ) %>%
  ggplot(aes(x = position, y = heterozygosity)) +
  geom_line(color = "grey", alpha = 0.5) +
  geom_point(aes(colour = factor(sample)), size = 1) +  # color by sample
  scale_color_viridis(discrete = TRUE) +
  facet_grid(sample ~ chromosome, scales = "free_x") +  # facets by sample and chromosome
  labs(x = NULL, y = "Heterozygosity\n") +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(labels = comma) +
  theme_minimal() +
  theme(
    legend.position = "top",
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing.x = unit(0, "line"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.25)
  )

```

