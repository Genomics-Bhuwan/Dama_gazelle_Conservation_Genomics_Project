##### hPSMC
---
-!/bin/bash -l
- SBATCH --job-name=hPSMC_all
- SBATCH --time=120:00:00
- SBATCH --cpus-per-task=10
- SBATCH --mem=128G
- SBATCH --partition=batch
- SBATCH --output=hpsmc_run_%A_%a.log
- SBATCH --mail-type=END,FAIL
- SBATCH --mail-user=bistbs@miamioh.edu
---

###### Step 1 & 2: Haploidize BAMs and concatenate chromosome FASTAs
---
- Step A. I am haplodizing for sample Addra: SRR17134085.bam only. It is to be done indivdiual sample wise in first few steps as also recommended by the author: https://github.com/jacahill/hPSMC?tab=readme-ov-file
-  After this, I will do hte haplodizing for another sample Mohrr:SRR17134085.bam. Then I will run hPSMC in these two samples.
-  I will repeat the same for other combination with Addra vs Mohrr sub-species.
---

```bash
REF=/scratch/bistbs/Population_Genomic_Analysis/PSMC/Dama_gazelle_hifiasm-ULONT_primary.fasta
BAM=/scratch/bistbs/Population_Genomic_Analysis/PSMC/SRR17134085.bam
CHROM_FILE=/scratch/bistbs/Population_Genomic_Analysis/PSMC/chromosomes.txt
PU2FA=/scratch/bistbs/Population_Genomic_Analysis/PSMC/Chrom-Compare/pu2fa
OUTDIR=/scratch/bistbs/Population_Genomic_Analysis/hPSMC/85
MAX_COVERAGE=100
mkdir -p $OUTDIR
```
##### Haploidize each chromosome/scaffold
```bash
    module load bcftools-1.15
module load samtools-1.22.1
module load parallel
cd /scratch/bistbs/Population_Genomic_Analysis/PSMC

SAMPLE="SRR17129394"
REF="Dama_gazelle_hifiasm-ULONT_primary.fasta"
PU2FA="/scratch/bistbs/Population_Genomic_Analysis/PSMC/Chrom-Compare/pu2fa"

export SAMPLE REF PU2FA  # Export for parallel

# Run each chromosome in parallel using 24 threads
cat chromosomes.txt | parallel -j 24 --env SAMPLE,REF,PU2FA '
    echo "Processing $SAMPLE chromosome {}..."
    samtools mpileup -s -f $REF -q30 -Q30 -r {} $SAMPLE.bam | \
    $PU2FA -c {} -C 100 > haploidized/${SAMPLE}_chr{}.fa
'
```

###### For sample 87
```bash
#!/bin/bash -l
#SBATCH --job-name=hPSMC_sample87
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --partition=batch
#SBATCH --output=psmc_sample87.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

module load bcftools-1.15
module load samtools-1.22.1
module load parallel

cd /scratch/bistbs/Population_Genomic_Analysis/PSMC
mkdir -p haploidized

SAMPLE="SRR17134087"
REF="Dama_gazelle_hifiasm-ULONT_primary.fasta"
PU2FA="/scratch/bistbs/Population_Genomic_Analysis/PSMC/Chrom-Compare/pu2fa"

export SAMPLE REF PU2FA  # Export for parallel

# Run each chromosome in parallel using 24 threads
cat chromosomes.txt | parallel -j 24 '
    echo "Processing '"$SAMPLE"' chromosome {}..."
    samtools mpileup -s -f $REF -q30 -Q30 -r {} $SAMPLE.bam | \
    $PU2FA -c {} -C 100 > haploidized/${SAMPLE}_chr{}.fa
'
```

##### Step 3.  Combine the chromosomes for each sapmles independently.
```bash
cd /scratch/bistbs/Population_Genomic_Analysis/PSMC/haploidized

cat SRR17129394_chr*.fa > SRR17129394_all.fa
cat SRR17134087_chr*.fa > SRR17134087_all.fa
```
##### Step 4.  Run Python for two samples to combine.
- The python script for hPSMC got corrupted that is why I had to fix it and it started running. 
```bash
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/psmcfa_from_2_fastas_try.py \
    -b 10 \
    -m 5 \
    SRR17129394_all.fa \
    SRR17134087_all.fa \
    > hPSMC.psmcfa

```

##### Step 5. PSMC for hPSMC
```bash
#!/bin/bash -l
#SBATCH --job-name=hPSMC_all
#SBATCH --time=150:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=90G
#SBATCH --partition=batch
#SBATCH --output=hpsmc_run_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

## Demographic modelling using HPSMC
```bash
#!/bin/bash -l
#SBATCH --job-name=hPSMC_all
#SBATCH --time=150:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=90G
#SBATCH --partition=batch
#SBATCH --output=hpsmc_run_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

## Demographic modelling using HPSMC

/scratch/bistbs/Population_Genomic_Analysis/PSMC/psmc/psmc \
  -N25 -t15 -r5 -p "4+25*2+4+6" \
  -o /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output/hPSMC.psmc \
  /scratch/bistbs/Population_Genomic_Analysis/PSMC/haploidized/hPSMC.psmcfa
```
##### Step 6. Plotting.
```bash
# Required libraries
library(ggplot2)
library(dplyr)

# Load your PSMC parser
source("F:/Collaborative_Projects/Dama_Gazelle_Project/hPSMC/plotPsmc.r") ### This is the source code Bhuwan has. Or you could google it.

# Working directory
setwd("F:/Collaborative_Projects/Dama_Gazelle_Project/hPSMC/")

# Parameters
mu <- 2.96e-9
g  <- 2.85

# Input file — ONLY ONE
psmc_file <- "hPSMC.psmc"

# ----------------------------
# READ PSMC RESULT USING AVAILABLE FUNCTION
# ----------------------------
res <- psmc.result(
  file = psmc_file,
  mu = mu,
  g = g,
  i.iteration = 25
)

df  <- bind_rows(res, .id = "iter")

# Add labels (only one)
df$SampleID   <- "hPSMC"
df$Subspecies <- "NA"
df$Label      <- "hPSMC"

# Keep only the main estimate (iteration = 1)
df_main <- df %>% filter(iter == "1")

# Remove first 8 points (PSMC burn-in)
df_main <- df_main %>%
  slice(9:n())

# Line color
line_color <- "#0072B2"

# ----------------------------
# PLOT WITH CUSTOM AXES
# ----------------------------
p_main <- ggplot(df_main, aes(x=YearsAgo, y=Ne)) +
  geom_step(linewidth = 1.6, direction = "hv", color = line_color) +
  
  scale_x_log10(
    limits = c(2e4, 2e6),  # X-axis from 20 Kya to 2 Mya
    breaks = c(2e4,5e4,1e5,5e5,1e6,2e6),
    labels = c('20 Kya','50 Kya','100 Kya','500 Kya','1 Mya','2 Mya')
  ) +
  scale_y_log10(
    limits = c(7e3, 1.1e5),  # Y-axis from 7k to 110k
    breaks = c(7e3,1e4,2e4,5e4,1e5),
    labels = c('7k','10k','20k','50k','100k')
  ) +
  
  annotation_logticks(sides = "bl") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour="black", fill=NA, linewidth=1),
    axis.title = element_text(size=14, face="bold"),
    axis.text = element_text(size=12)
  ) +
  labs(
    x = "Time",
    y = "Effective population size (Ne)",
    title = "PSMC Plot – hPSMC"
  ) +
  annotate(
    "text",
    x = 2e6,
    y = 8e3,  # slightly above bottom
    label = "(μ = 2.96e-9 & g = 5.85)",
    hjust = 1,
    vjust = 0,
    size = 5,
    fontface = "bold"
  )

print(p_main)

# Save plots
ggsave("hPSMC_PSMC_Main_20K_2M_7k_110k.pdf", p_main, width=14, height=8)
ggsave("hPSMC_PSMC_Main_20K_2M_7k_110k.jpeg", p_main, width=14, height=8, dpi=300)
```

##### Step 7. Pre-divergence analysis
- Use this paper as reference for looking the script for hPSMC.
- https://academic.oup.com/zoolinnean/article/204/3/zlaf059/8194500#525155301
- Ancestral $N_e$ (Observed): $55,000$Mutation Rate ($\mu$): $2.96 \times 10^{-9}$ per site per generationGeneration Time ($g$): $5.85$ yearsTotal Sequence Length ($L$): $255,813,793$ basesTarget Divergence Range: 1 kya to 4 Mya



- Run simulations of divergence without post-divergence migration to compare to the hPSMC plot.
- To compare and interpret hPSMC results, we need to compare our data to simulations.
- Since, the visual interpretation of hPSMC plots is suspectible to user bias and not replicable.
- To conduct simulations, the user should estimate the ancestral population size and estimate the a recent and ancient bound for when the sample might have diverged.
- Run it uing Python2 not Python3. I used Python 2.7.5


```bash
#!/bin/bash -l
#SBATCH --job-name=hPSMC_all
#SBATCH --time=150:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=90G
#SBATCH --partition=batch
#SBATCH --output=hpsmc_run_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# PARAMETER SUMMARY:
# N_e = 55,000, mu = 2.96e-9, g = 5.85 years
# L = 255,813,793 bases (I counted teh sequence length for the .psmcfa and used this length)
# Scaled Mutation Rate: -t 166583.8
# Divergence Range: 1 kya to 4 Mya

---

### Begin MSPRIME/ms simulations (10 simulations) ###
# CORRECTED PATH FOR mspms 

/home/users/bistbs/miniconda3/bin/mspms 4 40 -p 8 -t 166583.8 -r 60000.0 5000000 -I 2 2 2 -ej 0.000777 2 1 > ./hPSMC_sim_1kya.ms_sim &
/home/users/bistbs/miniconda3/bin/mspms 4 40 -p 8 -t 166583.8 -r 60000.0 5000000 -I 2 2 2 -ej 0.003885 2 1 > ./hPSMC_sim_5kya.ms_sim &
/home/users/bistbs/miniconda3/bin/mspms 4 40 -p 8 -t 166583.8 -r 60000.0 5000000 -I 2 2 2 -ej 0.011655 2 1 > ./hPSMC_sim_15kya.ms_sim &
/home/users/bistbs/miniconda3/bin/mspms 4 40 -p 8 -t 166583.8 -r 60000.0 5000000 -I 2 2 2 -ej 0.03108 2 1 > ./hPSMC_sim_40kya.ms_sim &
/home/users/bistbs/miniconda3/bin/mspms 4 40 -p 8 -t 166583.8 -r 60000.0 5000000 -I 2 2 2 -ej 0.0777 2 1 > ./hPSMC_sim_100kya.ms_sim &
/home/users/bistbs/miniconda3/bin/mspms 4 40 -p 8 -t 166583.8 -r 60000.0 5000000 -I 2 2 2 -ej 0.19425 2 1 > ./hPSMC_sim_250kya.ms_sim &
/home/users/bistbs/miniconda3/bin/mspms 4 40 -p 8 -t 166583.8 -r 60000.0 5000000 -I 2 2 2 -ej 0.3885 2 1 > ./hPSMC_sim_500kya.ms_sim &
/home/users/bistbs/miniconda3/bin/mspms 4 40 -p 8 -t 166583.8 -r 60000.0 5000000 -I 2 2 2 -ej 0.777 2 1 > ./hPSMC_sim_1mya.ms_sim &
/home/users/bistbs/miniconda3/bin/mspms 4 40 -p 8 -t 166583.8 -r 60000.0 5000000 -I 2 2 2 -ej 1.554 2 1 > ./hPSMC_sim_2mya.ms_sim &
/home/users/bistbs/miniconda3/bin/mspms 4 40 -p 8 -t 166583.8 -r 60000.0 5000000 -I 2 2 2 -ej 3.108 2 1 > ./hPSMC_sim_4mya.ms_sim &
wait

---
---
### Convert ms to psmcfa format (10 conversions) ###
python2 /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_1kya.ms_sim > ./hPSMC_sim_1kya.ms_sim.psmcfa &
python2 /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_5kya.ms_sim > ./hPSMC_sim_5kya.ms_sim.psmcfa &
python2 /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_15kya.ms_sim > ./hPSMC_sim_15kya.ms_sim.psmcfa &
python2 /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_40kya.ms_sim > ./hPSMC_sim_40kya.ms_sim.psmcfa &
python2 /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_100kya.ms_sim > ./hPSMC_sim_100kya.ms_sim.psmcfa &
python2 /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_250kya.ms_sim > ./hPSMC_sim_250kya.ms_sim.psmcfa &
python2 /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_500kya.ms_sim > ./hPSMC_sim_500kya.ms_sim.psmcfa &
python2 /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_1mya.ms_sim > ./hPSMC_sim_1mya.ms_sim.psmcfa &
python2 /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_2mya.ms_sim > ./hPSMC_sim_2mya.ms_sim.psmcfa &
python2 /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_4mya.ms_sim > ./hPSMC_sim_4mya.ms_sim.psmcfa &
wait

---
---
### Run PSMC (10 PSMC runs) ###
/scratch/bistbs/Population_Genomic_Analysis/PSMC/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ./hPSMC_sim_1kya.ms_sim.psmc ./hPSMC_sim_1kya.ms_sim.psmcfa &
/scratch/bistbs/Population_Genomic_Analysis/PSMC/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ./hPSMC_sim_5kya.ms_sim.psmc ./hPSMC_sim_5kya.ms_sim.psmcfa &
/scratch/bistbs/Population_Genomic_Analysis/PSMC/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ./hPSMC_sim_15kya.ms_sim.psmc ./hPSMC_sim_15kya.ms_sim.psmcfa &
/scratch/bistbs/Population_Genomic_Analysis/PSMC/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ./hPSMC_sim_40kya.ms_sim.psmc ./hPSMC_sim_40kya.ms_sim.psmcfa &
/scratch/bistbs/Population_Genomic_Analysis/PSMC/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ./hPSMC_sim_100kya.ms_sim.psmc ./hPSMC_sim_100kya.ms_sim.psmcfa &
/scratch/bistbs/Population_Genomic_Analysis/PSMC/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ./hPSMC_sim_250kya.ms_sim.psmc ./hPSMC_sim_250kya.ms_sim.psmcfa &
/scratch/bistbs/Population_Genomic_Analysis/PSMC/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ./hPSMC_sim_500kya.ms_sim.psmc ./hPSMC_sim_500kya.ms_sim.psmcfa &
/scratch/bistbs/Population_Genomic_Analysis/PSMC/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ./hPSMC_sim_1mya.ms_sim.psmc ./hPSMC_sim_1mya.ms_sim.psmcfa &
/scratch/bistbs/Population_Genomic_Analysis/PSMC/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ./hPSMC_sim_2mya.ms_sim.psmc ./hPSMC_sim_2mya.ms_sim.psmcfa &
/scratch/bistbs/Population_Genomic_Analysis/PSMC/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ./hPSMC_sim_4mya.ms_sim.psmc ./hPSMC_sim_4mya.ms_sim.psmcfa &
wait
---
---

### Estimate Divergence time with hPSMC ###
```bash
ls hPSMC_sim_*psmc | \
python2 /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/compare_sims_to_data.py \
    -i /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output/hPSMC.psmc \
    -N 15000 \
    -g 5 \
    -o hPSMC_gazelle_assessment_FINAL.txt

```





