# hPSMC

**Author:** Bhuwan Singh Bist

**Affiliation:** Jezkova lab

**Date:** 2025-10-11

#!/bin/bash -l
#SBATCH --job-name=hPSMC_all
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --partition=batch
#SBATCH --output=hpsmc_run_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# -----------------------------
# Paths to programs
# -----------------------------
PSMC=/scratch/bistbs/Population_Genomic_Analysis/PSMC/psmc/psmc
HPSMC_DIR=/scratch/bistbs/Population_Genomic_Analysis/hPSMC/hPSMC  # scripts folder
MS=/path/to/ms  # path to ms simulation program
REF=Dama_gazelle_hifiasm-ULONT_primary.fasta
PSMC_PLOT=/scratch/bistbs/Population_Genomic_Analysis/PSMC/psmc/utils/psmc_plot.pl

# -----------------------------
# Step 1 & 2: Haploidize BAMs and concatenate chromosome FASTAs
# -----------------------------
SAMPLES=("SRR17129394" "SRR17134085" "SRR17134086" "SRR17134087" "SRR17134088")

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Haploidizing $SAMPLE ..."
    for CHR in $(cat chromosomes.txt); do
        samtools mpileup -s -f $REF -q30 -Q60 -r $CHR ${SAMPLE}.bam | \
        pu2fa -c $CHR -C 100 > ${SAMPLE}_haploidized_${CHR}.fa
    done
    # Concatenate all chromosomes into one genome fasta
    cat ${SAMPLE}_haploidized_*.fa > ${SAMPLE}_all.fa
done

# -----------------------------
# Step 3: Generate hPSMC .psmcfa for all Addra × Mohrr pairs
# -----------------------------
ADDRA=("SRR17129394" "SRR17134085" "SRR17134086")
MOHR=("SRR17134087" "SRR17134088")

for A in "${ADDRA[@]}"; do
    for M in "${MOHR[@]}"; do
        echo "Generating hPSMC .psmcfa for $A × $M ..."
        python $HPSMC_DIR/psmcfa_from_2_fastas.py -b10 -m5 ${A}_all.fa ${M}_all.fa > hPSMC_${A}_${M}.psmcfa
    done
done

# -----------------------------
# Step 4: Run PSMC on all pairwise .psmcfa
# -----------------------------
for A in "${ADDRA[@]}"; do
    for M in "${MOHR[@]}"; do
        echo "Running PSMC for $A × $M ..."
        $PSMC -N25 -t15 -r5 -p "4+25*2+4+6" -o hPSMC_${A}_${M}.psmc hPSMC_${A}_${M}.psmcfa
    done
done

# -----------------------------
# Step 5: Plot hPSMC results
# -----------------------------
for A in "${ADDRA[@]}"; do
    for M in "${MOHR[@]}"; do
        echo "Plotting hPSMC for $A × $M ..."
        $PSMC_PLOT -g 5.85 -u 1.2e-8 -X 1000000 hPSMC_${A}_${M} hPSMC_${A}_${M}.psmc
    done
done

# -----------------------------
# Step 6: Estimate pre-divergence Ne
# -----------------------------
for A in "${ADDRA[@]}"; do
    for M in "${MOHR[@]}"; do
        echo "Estimating pre-divergence Ne for $A × $M ..."
        python $HPSMC_DIR/PSMC_emit_last_iteration_coord.py \
          -s10 -g5.85 -m1.2e-8 hPSMC_${A}_${M}.psmc > hPSMC_${A}_${M}_preNe.txt
    done
done

# -----------------------------
# Step 7 (optional): Run divergence simulations
# -----------------------------
# Replace <PRE_DIV_NE>, <LOWER_TIME>, <UPPER_TIME> with values estimated from Step 6
# for A in "${ADDRA[@]}"; do
#     for M in "${MOHR[@]}"; do
#         echo "Running divergence simulations for $A × $M ..."
#         python $HPSMC_DIR/hPSMC_quantify_split_time.py \
#           -o ./hPSMC_sim_${A}_${M} \
#           -N <PRE_DIV_NE> -l <LOWER_TIME> -u <UPPER_TIME> \
#           -s 10 -p 2 -P $PSMC -m $MS -H $HPSMC_DIR
#     done
# done

echo "hPSMC pipeline completed for all Addra × Mohrr pairs."


