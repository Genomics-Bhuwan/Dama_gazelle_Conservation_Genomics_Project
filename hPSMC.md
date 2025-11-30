##### hPSMC
---
#!/bin/bash -l
#SBATCH --job-name=hPSMC_all
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --partition=batch
#SBATCH --output=hpsmc_run_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
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
```bash
cd /scratch/bistbs/Population_Genomic_Analysis/PSMC/haploidized

python psmcfa_from_2_fastas.py -b10 -m5 SRR17129394_all.fa SRR17134087_all.fa > hPSMC.psmcfa
```
