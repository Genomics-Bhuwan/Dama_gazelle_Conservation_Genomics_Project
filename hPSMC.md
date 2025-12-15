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
- Run it uing Python2 not Python3. I used Python 2.7.5
- Just change Ne, mutation rate and generation time only. Rest keep same. They are intentionally kept there.

#### This is the hPSMC_quantify_split_time_Dama_gazelle.py that needs to be run in python2.
##### This is to be done for 94-87 sample.
```bash
### General Structure ###
## 1) Run hPSMC
## 2) Look at hPSMC output, what is the approximate Ne before it begins to rise to infinite, what is the approximate time range
## 3) Simulate population of presplit Ne with range of split times covering the observed split
## 4) Determine which simulations bracket the hPSMC plot
###
## We will start here with phase 3

### IMPORTS ###
import sys,os
from optparse import OptionParser
import datetime

### DEFAULTS ###
## Set by Flag ##
out="./hPSMC_sim_"
Ne=55000
lower=0
upper=10000000
sims=11
par=1
PSMC="psmc"
ms="ms"
hPSMC="./"

## Hard Coded ##
## ms command ##
N_CHRMS=4
N_REPS=40
GENERATION_TIME=5.85
N_SITES=5000000
MU=0.00000000296  # 10^-9  Mutation rate per year
RECOMB_S=0.00000001  # 10^-8 1CM/MB per site per generation 

## ms2psmcfa command ##
BIN_SIZE=10

### OPTIONS ###
usage = "usage: python %prog [options] input.psmc > outfile.sh"
parser = OptionParser(usage=usage)

parser.add_option("-o", "--out", type="str", dest="out", help="output directory for simulations and prefix all files for the run, default=\"./hPSMC_sim_\"")
parser.add_option("-N", "--Ne", type="int", dest="Ne", help="The ancestral population size to simulate, default=10,000")
parser.add_option("-l", "--lower", type="int", dest="lower", help="lower bound for simulations, the most recent divergence time to be simulated")
parser.add_option("-u", "--upper", type="int", dest="upper", help="upper bound for simulations, the most ancient divergence time to be simulated")
parser.add_option("-s", "--sims", type="int", dest="sims", help="the number of simulations to conduct, simulations will evenly split between high and low, minimum value=2, minimum meaningful value=3")
parser.add_option("-p", "--parallel", type="int", dest="par", help="Number of simulations to run simultaneously")
parser.add_option("-P", "--PSMC", type="str", dest="PSMC", help="If the psmc executable is not in your path give it's location, default = \"psmc\"")
parser.add_option("-m", "--ms", type="str", dest="ms", help="If the ms executable is not in your path give it's location, default = \"ms\"")
parser.add_option("-H", "--hPSMC", type="str", dest="hPSMC", help="If the hPSMC directory is not in your path give it's location, NOTE:  Just the directory not the script.  default = \"./\"") # make a def to check for "/"

(options, args) = parser.parse_args()

### SET OPTIONS ###
if options.out != None:
	out=options.out
if options.Ne != None:
	Ne=options.Ne
if options.lower != None:
	lower=options.lower
if options.upper != None:
	upper=options.upper
if options.sims != None:
	sims=options.sims
	if sims<2:
		sys.stderr.write("Error:  at least 2 simulations are required to estimate divergence time, please reset -s flag to an int >= 2, exiting")
		sys.exit()
	elif sims<3:
		sys.stderr.write("Warning:  This test will only determine whether the input hPSMC run falls between the lower and upper bound.  For more detailed estimation increase the number of simulations with -s")
if options.par != None:
	par=options.par
	if par<1:
		sys.stderr.write("Hint:  For systems with multiple processors increasing the number of simulations run in parallel will improve runtime.  use the -p flag to reset")
if options.PSMC != None:
	PSMC=options.PSMC
if options.ms != None:
	ms=options.ms
if options.hPSMC != None:
	hPSMC=options.hPSMC


### FUNCTIONS ###	

def ms_command(YEARS):
	THETA = 4 * Ne * GENERATION_TIME * MU * N_SITES	
	RECOMB = 4 * Ne * RECOMB_S * N_SITES
	SPLIT = (float(YEARS) / GENERATION_TIME) / (4 * float(Ne))
	MS_NAME = out+str(YEARS)+".ms_sim"
	print ms, N_CHRMS, N_REPS, "-p 8 -t", THETA, "-r", RECOMB, N_SITES, "-I 2 2 2 -ej", SPLIT, "2 1 >", MS_NAME, "&"
	return MS_NAME



### BODY ###

## make header ###

now = datetime.datetime.now()
print "#!/bin/bash/"
print ""
print "### File created:", now.strftime("%Y-%m-%d %H:%M")
print "### Using Command:", " ".join(sys.argv)
print ""


## build ms simulations ##
print "### Begin ms simulations ###"
step_size=(float(upper)-float(lower))/float(sims-1)
i=0
sim_names=[]   # I'm saving the simulation file names to avoid rounding errors
while i<sims:
	YEARS=lower+int(i*step_size)
	sim_names.append(ms_command(YEARS))
	i+=1
	if i%par==0:
		print "wait"
if i%par!=0:
	print "wait"
print ""


## convert ms to psmcfa ##
print "### Convert ms to psmcfa format ###"
i=0
while i<sims:
	sim_out="python "+hPSMC+"hPSMC_ms2psmcfa.py -b10 -d -c"+str(N_SITES)+" "+sim_names[i]+" > "+sim_names[i]+".psmcfa &"
	print sim_out
	i+=1
	if i%par==0:
		print "wait"
if i%par!=0:
	print "wait"
print ""

## run psmc ##
print "### Run PSMC ###"
i=0
while i<sims:
	sim_out=PSMC+" -N25 -t15 -r5 -p \"4+25*2+4+6\" -o "+sim_names[i]+".psmc "+ sim_names[i] +".psmcfa &"
	print sim_out
	i+=1
	if i%par==0:
		print "wait"
if i%par!=0:
	print "wait"
print ""

## Estimate Divergence time ## 
print "### Estimate Divergence time with hPSMC ###"
command = "ls " + out + "*psmc | python " + hPSMC + "hPSMC_compare_sims_to_data.py -i "+ sys.argv[-1] + " > " + out + "result.txt" 
print command
```

##### As you run above command following will be generated automatically as we run this command below.
- Below command runs within few seconds and you get below code.
- I ran them individually so that I could track all the intermediate stages.
```bash
python2 /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/hPSMC_quantify_split_time_Dama_gazelle.py -s 20 -l 1000 -u 100000 -p 20 -N 55000
```
#### Following will be generated here.
```bash
## Begin ms simulations ###
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.000777000777001 2 1 > ./hPSMC_sim_1000.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.00482517482517 2 1 > ./hPSMC_sim_6210.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.00887412587413 2 1 > ./hPSMC_sim_11421.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.0129222999223 2 1 > ./hPSMC_sim_16631.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.0169712509713 2 1 > ./hPSMC_sim_21842.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.0210194250194 2 1 > ./hPSMC_sim_27052.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.0250683760684 2 1 > ./hPSMC_sim_32263.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.0291165501166 2 1 > ./hPSMC_sim_37473.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.0331655011655 2 1 > ./hPSMC_sim_42684.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.0372136752137 2 1 > ./hPSMC_sim_47894.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.0412626262626 2 1 > ./hPSMC_sim_53105.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.0453108003108 2 1 > ./hPSMC_sim_58315.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.0493597513598 2 1 > ./hPSMC_sim_63526.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.0534079254079 2 1 > ./hPSMC_sim_68736.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.0574568764569 2 1 > ./hPSMC_sim_73947.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.0615050505051 2 1 > ./hPSMC_sim_79157.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.065554001554 2 1 > ./hPSMC_sim_84368.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.0696021756022 2 1 > ./hPSMC_sim_89578.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.0736511266511 2 1 > ./hPSMC_sim_94789.ms_sim &
mspms 4 40 -p 8 -t 19047.6 -r 11000.0 5000000 -I 2 2 2 -ej 0.0777000777001 2 1 > ./hPSMC_sim_100000.ms_sim &
wait
```

#### Convert ms to psmcfa format ###
```bash
# First, make sure the output folder exists
mkdir -p /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca

python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_1000.ms_sim > ./ms2psmca/hPSMC_sim_1000.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_6210.ms_sim > ./ms2psmca/hPSMC_sim_6210.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_11421.ms_sim > ./ms2psmca/hPSMC_sim_11421.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_16631.ms_sim > ./ms2psmca/hPSMC_sim_16631.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_21842.ms_sim > ./ms2psmca/hPSMC_sim_21842.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_27052.ms_sim > ./ms2psmca/hPSMC_sim_27052.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_32263.ms_sim > ./ms2psmca/hPSMC_sim_32263.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_37473.ms_sim > ./ms2psmca/hPSMC_sim_37473.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_42684.ms_sim > ./ms2psmca/hPSMC_sim_42684.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_47894.ms_sim > ./ms2psmca/hPSMC_sim_47894.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_53105.ms_sim > ./ms2psmca/hPSMC_sim_53105.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_58315.ms_sim > ./ms2psmca/hPSMC_sim_58315.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_63526.ms_sim > ./ms2psmca/hPSMC_sim_63526.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_68736.ms_sim > ./ms2psmca/hPSMC_sim_68736.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_73947.ms_sim > ./ms2psmca/hPSMC_sim_73947.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_79157.ms_sim > ./ms2psmca/hPSMC_sim_79157.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_84368.ms_sim > ./ms2psmca/hPSMC_sim_84368.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_89578.ms_sim > ./ms2psmca/hPSMC_sim_89578.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_94789.ms_sim > ./ms2psmca/hPSMC_sim_94789.ms_sim.psmcfa &
python /scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC/ms2psmcfa.py -b10 -d -c5000000 ./hPSMC_sim_100000.ms_sim > ./ms2psmca/hPSMC_sim_100000.ms_sim.psmcfa &
wait

```
### Run PSMC ###
```bash
/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_1000.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_1000.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_6210.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_6210.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_11421.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_11421.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_16631.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_16631.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_21842.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_21842.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_27052.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_27052.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_32263.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_32263.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_37473.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_37473.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_42684.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_42684.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_47894.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_47894.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_53105.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_53105.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_58315.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_58315.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_63526.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_63526.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_68736.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_68736.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_73947.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_73947.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_79157.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_79157.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_84368.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_84368.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_89578.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_89578.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_94789.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_94789.ms_sim.psmcfa

/scratch/bistbs/Population_Genomic_Analysis/hPSMC/Output_94_87/94_87simulation/ms2psmca/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
    -o /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/psmc_output/hPSMC_sim_100000.ms_sim.psmc \
    /scratch/bistbs/Population_Genomic_Analysis/hPSHC/Output_94_87/94_87simulation/ms2psmca/hPSMC_sim_100000.ms_sim.psmcfa

```
### Estimate Divergence time with hPSMC ###
```bash
ls ./hPSMC_sim_*psmc | python ./hPSMC_compare_sims_to_data.py -i 55000 > ./hPSMC_sim_result.txt
```

##### Do the same as you did for above two samples.
- For sample 85 and 88, I will use the same above simulated result with the empirical psmc file with the pairs: 94 and 87(Addra 1 vs. Mhorr1) and 85 and 88(Addra 2 vs Mhorr 2)
