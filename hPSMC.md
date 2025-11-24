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
- Step A. I am haplodizing for sample SRR17134085.bam only. It is to be done indivdiual sample wise in first few steps as also recommended by the author: https://github.com/Genomics-Bhuwan/Dama_gazelle_Conservation_Genomics_Project/edit/main/hPSMC.md
---

```bash
REF=/scratch/bistbs/Population_Genomic_Analysis/PSMC/Dama_gazelle_hifiasm-ULONT_primary.fasta
BAM=/scratch/bistbs/Population_Genomic_Analysis/PSMC/SRR17134085.bam
CHROM_FILE=/scratch/bistbs/Population_Genomic_Analysis/PSMC/chromosomes.txt
PU2FA=/scratch/bistbs/Population_Genomic_Analysis/PSMC/Chrom-Compare/pu2fa
OUTDIR=/scratch/bistbs/Population_Genomic_Analysis/hPSMC/85
MAX_COVERAGE=100

mkdir -p $OUTDIR

# -----------------------------
# Haploidize each chromosome/scaffold
# -----------------------------
for CHR in $(cat $CHROM_FILE); do
    echo "Processing $CHR ..."
    
    samtools mpileup -s -f $REF -q20 -Q20 -r $CHR $BAM | \
    $PU2FA -c $CHR -C $MAX_COVERAGE > $OUTDIR/SRR17134085_haploidized_${CHR}.fa

    if [ ! -s "$OUTDIR/SRR17134085_haploidized_${CHR}.fa" ]; then
        echo "Warning: $CHR output is empty"
    fi
done

# -----------------------------
# Concatenate all chromosomes into a single fasta
# -----------------------------
cat $OUTDIR/SRR17134085_haploidized_*.fa > $OUTDIR/SRR17134085_all.fa
echo "Haploidization complete: $OUTDIR/SRR17134085_all.fa"

```


##### Step 3: Generate hPSMC .psmcfa for all Addra × Mohrr pairs
```bash
ADDRA=("SRR17129394" "SRR17134085" "SRR17134086")
MOHR=("SRR17134087" "SRR17134088")

for A in "${ADDRA[@]}"; do
    for M in "${MOHR[@]}"; do
        echo "Generating hPSMC .psmcfa for $A × $M ..."
        python $HPSMC_DIR/psmcfa_from_2_fastas.py -b10 -m5 ${A}_all.fa ${M}_all.fa > hPSMC_${A}_${M}.psmcfa
    done
done
```
##### Step 4: Run PSMC on all pairwise .psmcfa
```bash
for A in "${ADDRA[@]}"; do
    for M in "${MOHR[@]}"; do
        echo "Running PSMC for $A × $M ..."
        $PSMC -N25 -t15 -r5 -p "4+25*2+4+6" -o hPSMC_${A}_${M}.psmc hPSMC_${A}_${M}.psmcfa
    done
done
```
##### Step 5: Plot hPSMC results
```bash
for A in "${ADDRA[@]}"; do
    for M in "${MOHR[@]}"; do
        echo "Plotting hPSMC for $A × $M ..."
        $PSMC_PLOT -g 5.85 -u 1.2e-8 -X 1000000 hPSMC_${A}_${M} hPSMC_${A}_${M}.psmc
    done
done
```
##### Step 6: Estimate pre-divergence Ne
```bash
for A in "${ADDRA[@]}"; do
    for M in "${MOHR[@]}"; do
        echo "Estimating pre-divergence Ne for $A × $M ..."
        python $HPSMC_DIR/PSMC_emit_last_iteration_coord.py \
          -s10 -g5.85 -m1.2e-8 hPSMC_${A}_${M}.psmc > hPSMC_${A}_${M}_preNe.txt
    done
done
```
##### Step 7 (optional): Run divergence simulations
```bash
Replace <PRE_DIV_NE>, <LOWER_TIME>, <UPPER_TIME> with values estimated from Step 6
for A in "${ADDRA[@]}"; do
  for M in "${MOHR[@]}"; do
        echo "Running divergence simulations for $A × $M ..."
     python $HPSMC_DIR/hPSMC_quantify_split_time.py \
          -o ./hPSMC_sim_${A}_${M} \
          -N <PRE_DIV_NE> -l <LOWER_TIME> -u <UPPER_TIME> \
         -s 10 -p 2 -P $PSMC -m $MS -H $HPSMC_DIR
    done
done

echo "hPSMC pipeline completed for all Addra × Mohrr pairs."
```

