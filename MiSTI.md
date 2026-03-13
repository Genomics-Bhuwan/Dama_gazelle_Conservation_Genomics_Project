#### Using MiSTI for divergence between two species.
- https://github.com/Genomics-HSE/MiSTI

#### Step 1. Run PSMC (Li, Durbin 2011) on two genomes.
#### Step 2. Run ANGSD and SFS to get it and then combine. Generate joint SFS
```bash
#!/bin/bash -l
#SBATCH --job-name=ANGSD_MiSTI
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --partition=batch
#SBATCH --output=angsd_misti_%A.log

# ----------------------------
# Paths & References
# ----------------------------
OUT_DIR="/home/bistbs/Dama_gazelle_MiSTI_divergence/Step2_ANGSD"
REF="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Dama_gazelle_hifiasm-ULONT_primary.fasta"

cd $OUT_DIR

# ----------------------------
# STEP 1: Generate SAF files for each population
# ----------------------------
echo "[$(date)] Generating SAF for Addra..."
angsd -P 24 -bam addra_bams.txt -anc $REF -ref $REF -out addra_pop \
      -dosaf 1 -gl 1 -C 50 -minQ 20 -minmapq 30

echo "[$(date)] Generating SAF for Mhorr..."
angsd -P 24 -bam mhorr_bams.txt -anc $REF -ref $REF -out mhorr_pop \
      -dosaf 1 -gl 1 -C 50 -minQ 20 -minmapq 30

# ----------------------------
# STEP 2: Generate Joint SFS (realSFS)
# ----------------------------
echo "[$(date)] Generating Joint SFS..."
realSFS addra_pop.saf.idx mhorr_pop.saf.idx -P 24 > addra_mhorr.real.sfs

echo "[$(date)] Process completed."
```
#### Step 3. Convert joint realSFS to MiSTI file format.
```bash
python ./utils/ANGSDSFS.py mhorr_addra.real.sfs Mhorr Addra > mhorr_addra.mi.sfs
```

#### Step 4. Align the Time Scales
- Since, we have 3 Addra vs 2 Mhorr. We assesed the index time interval for each pair.
- MiSTI merges the discrete time windows from your two separate PSMC runs into a single timeline.
```bash
python3 ~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/MiSTI/utils/calc_time.py \
~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/SRR17129394.psmc \
~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/SRR17134087.psmc \
> ~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/Index_files/time_index_29394_34087.txt

python3 ~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/MiSTI/utils/calc_time.py \
~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/SRR17129394.psmc \
~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/SRR17134088.psmc \
> ~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/Index_files/time_index_29394_34088.txt

python3 ~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/MiSTI/utils/calc_time.py \
~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/SRR17134085.psmc \
~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/SRR17134087.psmc \
> ~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/Index_files/time_index_34085_34087.txt

python3 ~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/MiSTI/utils/calc_time.py \
~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/SRR17134085.psmc \
~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/SRR17134088.psmc \
> ~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/Index_files/time_index_34085_34088.txt

python3 ~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/MiSTI/utils/calc_time.py \
~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/SRR17134086.psmc \
~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/SRR17134087.psmc \
> ~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/Index_files/time_index_34086_34087.txt

python3 ~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/MiSTI/utils/calc_time.py \
~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/SRR17134086.psmc \
~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/SRR17134088.psmc \
> ~/Dama_gazelle_MiSTI_divergence/Step2_ANGSD/Real_SFS/PSMC_files/Index_files/time_index_34086_34088.txt
```

#### Step 5. 
- Run the MiSTI Inference
- If you are not sure of excat split time, it is common practise to run several version with different split_time indices, and compare the likelihood.
-Choose the one with the lowest likelihood for runing the MiSTI.
- OR What I did is:- Looked into PSMC time of divergence. selected the range of 200k-300k.
- I did this ~250,000 years, generation time g=5.85: which is around 42,735 generations.
- I will look for time index that correspond to the 42725 generation.

  #############################################
  - 29394–34087- Index is 64
  - 29394_34088- 64
  - 34085_34087- 64
  - 34085_34088- 64
  - 34086_34087- 64
  - 34086_34088- 64

#### Step 6. I want to see the time of divergence of the indivudal.

```bash
python ./MiSTI.py mhorr.psmc addra.psmc mhorr_addra.mi.sfs <split_time_index> -o results.mi
```
#### Step 7. But, if you want to include "MIGRATION" or want to see the gene flow after divergence.
- If you suspect gene flow after the split; use the -mi flag.
```bash
python ./MiSTI.py mhorr.psmc addra.psmc mhorr_addra.mi.sfs 15 \
    -mi 2 2 10 0.5 1 -o results_with_mig.mi
```

#### Step 8. Plot the Ne and split time.
```bash
python ./MiSTIPlot.py mhorr.psmc addra.psmc results.mi -o gazelle_divergence_plot.pdf
```
