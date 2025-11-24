##### GONE2 for estimating the past Ne using LD method.
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
#SBATCH --job-name=Addra_GONe
#SBATCH --time=300:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=batch
#SBATCH --output=logs/Addra_%A.out
#SBATCH --error=logs/Addra_%A.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu


###RUn
module purge
module load gcc-14.2.0

### Run again
cd /scratch/bistbs/Population_Genomic_Analysis/GONE2/GONE2
make clean

###Since, the SNPs were greater than 20000, therefore, I had to increase the limit to this. Coz, SNPs in my file were aroudn 10-11 millions.
make MAXLOCI=10000000 MAXIND=3 gone

### Run the GONe program for the Addra gazelle
/scratch/bistbs/Population_Genomic_Analysis/GONE2/GONE2/gone2 \
    -g 0 -r 1.1 -t 16 \
    /scratch/bistbs/Population_Genomic_Analysis/GONE2/Dama_gazelle_Addra.vcf \
    -o /scratch/bistbs/Population_Genomic_Analysis/GONE2/Output_Addra/Dama_gazelle_Addra
```

- Run for Mhorr and save outputs in Output_Mhorr
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


###RUn
module purge
module load gcc-14.2.0

### Run again
cd /scratch/bistbs/Population_Genomic_Analysis/GONE2/GONE2
make clean

###Since, the SNPs were greater than 20000, therefore, I had to increase the limit to this. Coz, SNPs in my file were aroudn 10-11 millions.
make MAXLOCI=10000000 MAXIND=2 gone


### Run the GONe program for the MOhrr gazelle
/scratch/bistbs/Population_Genomic_Analysis/GONE2/GONE2/gone2 \
    -g 0 -r 1.1 -t 8 \
    /scratch/bistbs/Population_Genomic_Analysis/GONE2/Dama_gazelle_Mhorr.vcf \
    -o /scratch/bistbs/Population_Genomic_Analysis/GONE2/Output_Mhorr/Dama_gazelle_Mhorr
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
