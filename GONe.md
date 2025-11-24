##### GONE2 for estimating the past Ne using LD method.https://github.com/esrud/GONE2
##### Go the the path where there is vcf file.
- Since, I have a vcf file with two subspecies: Addra and Mohrr.
- I am splitting the vcf file based on sub-species and will run the GONE2 seperately.

##### Path to the GONE2
```bash
/scratch/bistbs/Population_Genomic_Analysis/GONE2
```
##### For addra
- Addra
```bash
bcftools view -S <(echo -e "SRR17129394\nSRR17134085\nSRR17134086") \
    -o Dama_gazelle_Addra.vcf \
    /scratch/bistbs/Population_Genomic_Analysis/GONE2/Dama_gazelle_biallelic_snps_autosomes.vcf
```

##### Mhorr
```bash
bcftools view -S <(echo -e "SRR17134087\nSRR17134088") \
    -o Dama_gazelle_Mhorr.vcf \
    /scratch/bistbs/Population_Genomic_Analysis/GONE2/Dama_gazelle_biallelic_snps_autosomes.vcf
```

##### Run GONE2
- Run for Addra and save outputs in Output_Addra
```bash
#!/bin/bash -l
#SBATCH --job-name=Mohrr_GONe
#SBATCH --time=300:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=batch
#SBATCH --output=logs/MOhrr_%A.out
#SBATCH --error=logs/MOhrr_%A.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu


#!/bin/bash
#======================================================
# GONE2 Analysis for Addra Gazelle (All SNPs)
#======================================================

# File paths
VCF_FILE="/scratch/bistbs/Population_Genomic_Analysis/GONE2/Dama_gazelle_Addra.vcf"
GONE2_DIR="/scratch/bistbs/Population_Genomic_Analysis/GONE2/GONE2"
OUT_FILE="/scratch/bistbs/Population_Genomic_Analysis/GONE2/Output_Addra"

# Parameters
RECOMB_RATE=1.1   # cM/Mb
THREADS=16
GENOTYPE_TYPE=0   # 0 = unphased diploid

#---------------------------------------------
# Step 1: Compile GONE2 if not already compiled
#---------------------------------------------
cd $GONE2_DIR
if [ ! -f gone2 ]; then
    echo "Compiling GONE2..."
    g++ -fopenmp gone2.cpp lib/*.cpp -o gone2 -DMAXLOCI=9416350 -DMAXIND=3 
else
    echo "GONE2 already compiled."
fi

#---------------------------------------------
# Step 2: Prepare VCF file (all SNPs)
#---------------------------------------------
PREP_VCF="/scratch/bistbs/Population_Genomic_Analysis/GONE2/Addra_allSNPs.vcf"

echo "Preparing VCF file with all SNPs..."
grep "^#" $VCF_FILE > header.txt
grep -v "^#" $VCF_FILE >> header.txt
mv header.txt $PREP_VCF

echo "VCF prepared: $PREP_VCF"

#---------------------------------------------
# Step 3: Run GONE2
#---------------------------------------------
echo "Running GONE2..."
$GONE2_DIR/gone2 -g $GENOTYPE_TYPE -r $RECOMB_RATE -t $THREADS $PREP_VCF -o $OUT_FILE &

echo "GONE2 job started in background."
echo "Monitor progress with:"
echo "  cat ${OUT_FILE}_GONE_progress.tmp"

```

# Recompile GONE2 for your number of individuals
- The reason why I am using MAXIND 10 is coz the maximum capacity should be greater than the number of individuals we have in the vcf file.
make clean
make MAXLOCI=9416350 MAXIND=10 gone


### Run the GONe program for the MOhrr gazelle
```bash
GONE2_DIR=/scratch/bistbs/Population_Genomic_Analysis/GONE2/Mohrr_gazelle/GONE2
PREP_VCF=/scratch/bistbs/Population_Genomic_Analysis/GONE2/Dama_gazelle_Mhorr.vcf
GENOTYPE_TYPE=0      # 0 = unphased diploid
RECOMB_RATE=1.1      # constant recombination rate in cM/Mb
THREADS=16
OUT_FILE=/scratch/bistbs/Population_Genomic_Analysis/GONE2/Output_Mhorr

# Run in background
$GONE2_DIR/gone2 -g $GENOTYPE_TYPE -r $RECOMB_RATE -t $THREADS $PREP_VCF -o $OUT_FILE &
echo "GONE2 job started in background."
echo "Monitor progress with:"
echo "  cat ${OUT_FILE}_GONE_progress.tmp"

```

###Visualization of the "Estimated Effective population size" against the generation in time using Linkage Disequilibrium method.

```bash
# Load required packages
library(tidyverse)

# Set file paths for GONe outputs
addra_file <- "/scratch/bistbs/Population_Genomic_Analysis/GONE2/Output_Addra/Dama_gazelle_Addra_GONE_Ne"
mhorr_file <- "/scratch/bistbs/Population_Genomic_Analysis/GONE2/Output_Mhorr/Dama_gazelle_Mhorr_GONE_Ne"

# Read the data
addra <- read.table(addra_file, header = TRUE)
mhorr <- read.table(mhorr_file, header = TRUE)

# Inspect first few rows to check column names
head(addra)
head(mhorr)

# Add species column
addra$Species <- "Addra"
mhorr$Species <- "Mhorr"

# Combine into one dataframe
ne_data <- rbind(addra, mhorr)

# Example: Assuming columns are named 'Generation' and 'Ne'
# If column names are different, replace them accordingly

# Plot Ne vs Generation
ggplot(ne_data, aes(x = Generation, y = Ne, color = Species)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  scale_y_log10() +  # optional: log scale if Ne varies widely
  scale_x_reverse() +  # generations often plotted backward in time
  labs(
    x = "Generation (backward in time)",
    y = "Estimated Effective Population Size (Ne)",
    title = "Historical Effective Population Size for Dama Gazelle Sub-species"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")
``
