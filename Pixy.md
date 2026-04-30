#### Dxy and Pxy within and between the population of Addra and Mhorr gazelle.
- Link: GITHUB: https://github.com/ksamuk/pixy
- https://pixy.readthedocs.io/en/latest/; 
- Paper: https://doi.org/10.1111/1755-0998.13326

#### Step 1. Get the monomorphic sites or invariants site as well as variants site.
```bash
#!/bin/bash -l
#SBATCH --time=75:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --partition=batch
#SBATCH --job-name=Monomorphic_site
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bistbs@miamioh.edu

# Move to working directory
cd /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Pixy
module load gatk-4.1.2.0 

# Define your paths for clarity
REFERENCE="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Dama_gazelle_hifiasm-ULONT_primary.fasta"
INPUT_GVCF="all_samples_combined.g.vcf.gz"
OUTPUT_VCF="all_samples_allsites_genotyped.vcf.gz"

# Run GATK GenotypeGVCFs
# We increase memory to 90G given the T2T assembly size
gatk --java-options "-Xmx90g" GenotypeGVCFs \
  -R $REFERENCE \
  -V $INPUT_GVCF \
  -O $OUTPUT_VCF \
  --include-non-variant-sites true
```


#### Step 2. Go for variant filtration using vcf-tools
```bash
vcftools --gzvcf /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Pixy/all_samples_allsites_genotyped.vcf.gz \
  --remove-indels \
  --max-missing 0.8 \
  --min-meanDP 10 \
  --max-meanDP 100 \
  --minQ 30 \
  --recode --stdout | gzip -c > filtered_allsites.vcf.gz
```

#### Step 3. Get the names of the samples of Addra and Mhorr
- Make a file and seperate them.
- Calcualte the pi dxy fst among these six pair of population using below code or revise it if needed.
  ---
Specifically, the dXY estimator can be used to estimate the absolute divergence between populations, along with Fst. I suggest we test the pairwise differences between each of the 3 addra and 2 mhorr gazelles (6x pairwise comparisons) as well as within each subspecies (between the 3 Addra and between the two Mhorr gazelles).
  ---


#### Step 4. Run Pixy and Dxy
- Run it based on chromosomes.
 #### Step 4.a Make the populations.txt
```bash
SRR17129394	Addra
SRR17134085	Addra
SRR17134086	Addra
SRR17134087	Mhorr
SRR17134088	Mhorr
```
#### Step 4.b. Run the pixy

```bash
pixy \
--vcf /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Pixy/VCF_filtration/filtered_all_sites.vcf.gz \
--populations populations.txt \
--stats pi dxy fst \
--chromosomes '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17' \
--window_size 10000 \
--n_cores 4 \
--output_prefix dama_pixy
```



#################################################################################################################################

#### Klaus recommended me to find out the outlier loci in the Dxy and correspond those loci or SNPs with the genes associated with them
```bash
- Identifying the genes in Dxy outlier regions to find out which genes are in those spikes on chromosomes 6 7 14 and 15.
- We need the GFF/GTF annotation file from our dama gazelle reference genome.
- Workflow:
- Filter outliers: We take the pxy_dxy.txt file and pull the top 1% or specific threshold of windows.
- Overlap with Annotation using a tool called bedtools intersect to see which gene names in the GFF file overlap with your high Dxy window coordinates. 
```
```bash
# ============================================================
# STEP 1: EXTRACT DXY OUTLIERS FOR GENE ANNOTATION
# ============================================================

# 1. Load required libraries
library(dplyr)

# 2. Define the file path (using the path you provided)
# Note: In R, use forward slashes / even on Windows
input_file <- "F:/Collaborative_Projects/Dama_Gazelle_Project/Pixy/dama_pixy_dxy.txt"
output_bed <- "F:/Collaborative_Projects/Dama_Gazelle_Project/Pixy/dxy_outliers.bed"

# 3. Load the data
# We use check.names=F in case there are odd characters in headers
dxy_data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 4. Filter out windows with NA values
# Pixy produces NAs where there are no sites to compare
dxy_clean <- dxy_data %>%
  filter(!is.na(avg_dxy))

# 5. Calculate the Outlier Threshold
# We take the top 1% of the data. 
# Based on your plot, these will capture those high spikes on Chr 6, 7, 14, 15.
threshold_value <- quantile(dxy_clean$avg_dxy, 0.99)

# 6. Extract the outliers
dxy_outliers <- dxy_clean %>%
  filter(avg_dxy >= threshold_value) %>%
  arrange(desc(avg_dxy)) # Sort by highest Dxy first

# 7. Format for BED file (Chrom, Start, End)
# BED format is 0-based for the start position, but Pixy is 1-based.
# To be safe for bedtools, we subtract 1 from window_pos_1.
bed_output <- dxy_outliers %>%
  mutate(start_0 = window_pos_1 - 1) %>%
  select(chromosome, start_0, window_pos_2)

# 8. Save the BED file
write.table(bed_output, 
            file = output_bed, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

# 9. Summary Report
cat("--- Dxy Outlier Extraction Summary ---\n")
cat("Total windows analyzed: ", nrow(dxy_clean), "\n")
cat("Outlier Threshold (99th percentile): ", threshold_value, "\n")
cat("Number of outlier windows found: ", nrow(dxy_outliers), "\n")
cat("Outliers saved to: ", output_bed, "\n\n")

# Print count per chromosome to verify those spikes
cat("Outlier count per chromosome:\n")
print(table(dxy_outliers$chromosome))
```

#### Intersection of the outliers vs. Annotation file of .gff/gtf files
```bash
# 1. Move to the directory
cd /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Pixy/Gene_Ontology/

# 2. Load bedtools (Check your cluster's specific name, often 'bedtools' or 'bedtools2')
module load bedtools

# 3. Run the intersection
# -a: Your outlier windows
# -b: Your genome annotation
# -wa: Keep the info from your BED file (Chr, Start, End)
# -wb: Append the info from the GFF file (Gene names, types, etc.)
bedtools intersect \
-a dxy_outliers.bed \
-b Addra_complete.genomic.gff \
-wa -wb > dxy_outlier_genes_full_info.txt


# 3.List the name of the genes
# This pulls the 'gene=' field from the GFF portion of the output
grep "gene=" dxy_outlier_genes_full_info.txt | \
awk -F'gene=' '{print $2}' | \
cut -d';' -f1 | \
sort | \
uniq > final_candidate_genes.txt

# 4. Annotate the genes using gene Ontology using ShinyGo

```







######################################################################################################################################################
######################################################################################################################################################
##### Running the pixy to estimte the Dxy and Fst between the Dama gazelles and Grant's gazelle. 
##### It will provide a baseline between two well-differentiated species to compare the Dxy and Fst values between these species versus between Addra and Mhorr gazelle.

#### Step 1. a Call for genotypes for variants sites as well as non-variants sites.
```bash
#!/bin/bash -l
#SBATCH --time=85:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=90G
#SBATCH --partition=batch
#SBATCH --job-name=AllSites
#SBATCH --output=allsites_%j.out
#SBATCH --error=allsites_%j.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bistbs@miamioh.edu

# Note: %j will automatically insert the unique Job ID into the filename

# Move to your current working directory
cd /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Pixy/Klaus_Pxy_recommendation

module load gatk-4.1.2.0 

# Define paths
REFERENCE="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Dama_gazelle_hifiasm-ULONT_primary.fasta"
INPUT_GVCF="Joint_Genotyping_combined.g.vcf.gz"
OUTPUT_VCF="dama_baseline_allsites.vcf.gz"

# Run GATK GenotypeGVCFs
# Including non-variant sites as requested by your collaborator
gatk --java-options "-Xmx90g" GenotypeGVCFs \
  -R $REFERENCE \
  -V $INPUT_GVCF \
  -O $OUTPUT_VCF \
  --include-non-variant-sites true
```
#### Step 2. Use vcf-tools for removing the ambigious calls.
```bash
#!/bin/bash
#SBATCH --job-name=vcf_filter
#SBATCH --output=vcf_filter_%j.out
#SBATCH --error=vcf_filter_%j.err
#SBATCH --time=70:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --partition=batch

# Load required module
module load vcf-tools

# Define input and output
INPUT_VCF="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Pixy/Klaus_Pxy_recommendation/Remove_one_sample/dama_baseline_filtered.vcf.gz"
OUTPUT_VCF="filtered_allsites.vcf.gz"

# Run filtering
vcftools --gzvcf $INPUT_VCF \
  --remove-indels \
  --max-missing 0.8 \
  --min-meanDP 10 \
  --max-meanDP 100 \
  --minQ 30 \
  --recode --stdout | gzip -c > $OUTPUT_VCF

echo "Filtering complete: $OUTPUT_VCF"
```
#### Step 3. Get the names of the samples of Addra, Mhorr and Grant's gazelle
- Make a file and seperate them.
- Calcualte the pi dxy fst among these six pair of population using below code or revise it if needed.
- Specifically, the dXY estimator can be used to estimate the absolute divergence between populations, along with Fst.
- I suggest we test the pairwise differences between each of the 3 addra and 2 mhorr gazelles (6x pairwise comparisons) as well as within each subspecies (between the 3 Addra and between the two Mhorr gazelles along with Grant's gazelle).
  
#### Step 4. Run Pixy and Dxy
- Run it based on chromosomes(I am running with 1-17 autosomes only).
 #### Step 4.a Make the  individual populations.txt for individual pariwise comparions.
```bash
SRR17129394	Addra1
SRR17134085	Addra2
SRR17134086	Addra3
SRR17134087	Mhorr1
SRR17134088	Mhorr2
SRR6878810  Grant's
SRR6894844	Thompson's

```
#### Step 4.b. Run the pixy

```bash
pixy \
mkdir pixy_results

--vcf dama_pixy_filtered.vcf.gz \
--populations indiv_populations.txt \
--stats pi dxy fst \
--chromosomes '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17' \
--window_size 10000 \
--n_cores 20 \
--output_folder ./pixy_results \
--output_prefix dama
```

