#### This is used for estimating the Effective Population Size using Linkage Disequilibrium method.
#### Load the modules.
module load R


# Run PLINK
plink --bfile /scratch/bistbs/Population_Genomic_Analysis/GONe/Dama_gazelle_biallelic_snps_filtered.recode \
      --keep /scratch/bistbs/Population_Genomic_Analysis/GONe/file_lists/Addra_gazelle.txt \
      --extract /scratch/bistbs/Population_Genomic_Analysis/GONe/temp_SNP_sub \
      --mac 1 \
      --recode \
      --not-chr NW_024070207.1 \
      --allow-extra-chr \
      --nonfounders \
      --debug \
      --out /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Addra/Dama_gazelle_biallelic_snps_filtered.recode_Addra

# Recode map file for GONe
Rscript /scratch/bistbs/Population_Genomic_Analysis/GONe/scripts/recode_map_gone.R \
        /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Addra/Dama_gazelle_biallelic_snps_filtered.recode_Addra.map

# Run GONe
bash /scratch/bistbs/Population_Genomic_Analysis/GONe/scripts/script_GONE.sh \
     /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Addra/Dama_gazelle_biallelic_snps_filtered.recode_Addra

# Move output files
mv /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Addra/outfileHWD \
   /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Addra/outfileLD_d2_sample \
   /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Addra/outfileLD_Ne_estimates \
   /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Addra/seedfile \
   /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Addra/timefile \
   /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Addra/

mv /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Addra/TEMPORARY_FILES \
   /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Addra/

echo "Addra analysis done!"


---
Do it for Mohrr gazelle
---


# Run PLINK
plink --bfile /scratch/bistbs/Population_Genomic_Analysis/GONe/Dama_gazelle_biallelic_snps_filtered.recode \
      --keep /scratch/bistbs/Population_Genomic_Analysis/GONe/file_lists/Mohrr_gazelle.txt \
      --extract /scratch/bistbs/Population_Genomic_Analysis/GONe/temp_SNP_sub \
      --mac 1 \
      --recode \
      --not-chr NW_024070207.1 \
      --allow-extra-chr \
      --nonfounders \
      --debug \
      --out /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Mohrr/Dama_gazelle_biallelic_snps_filtered.recode_Mohrr

# Recode map file for GONe
Rscript /scratch/bistbs/Population_Genomic_Analysis/GONe/scripts/recode_map_gone.R \
        /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Mohrr/Dama_gazelle_biallelic_snps_filtered.recode_Mohrr.map

# Run GONe
bash /scratch/bistbs/Population_Genomic_Analysis/GONe/scripts/script_GONE.sh \
     /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Mohrr/Dama_gazelle_biallelic_snps_filtered.recode_Mohrr

# Move output files
mv /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Mohrr/outfileHWD \
   /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Mohrr/outfileLD_d2_sample \
   /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Mohrr/outfileLD_Ne_estimates \
   /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Mohrr/seedfile \
   /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Mohrr/timefile \
   /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Mohrr/

mv /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Mohrr/TEMPORARY_FILES \
   /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Mohrr/

echo "Mohrr analysis done!"

