#### This is used for estimating the Effective Population Size using Linkage Disequilibrium method.
#### Load the modules.
module load R


# Run PLINK
```bash
plink --bfile /scratch/bistbs/Population_Genomic_Analysis/GONe/Dama_gazelle \
      --keep /scratch/bistbs/Population_Genomic_Analysis/GONe/Addra_gazelle.txt \
      --mac 1 \
      --recode \
      --allow-extra-chr \
      --nonfounders \
      --debug \
      --out /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Addra/Dama_gazelle_Addra
```


# Recode map file for GONe
### This is script from Emily: https://github.com/elhumble/SHO_reseq_2022/blob/master/scripts/recode_map_gone.R

```bash
library(data.table)
library(dplyr)

# Get input argument
args <- commandArgs(trailingOnly=TRUE)[[1]]

# Read map file
map <- fread(args)

# For your data, chromosomes are already 1-17, so no recoding needed
# Just write it back
write.table(map, args, col.names = FALSE, row.names = FALSE, quote = FALSE)
```
##### Run the recode_map_gone.R
```bash
Rscript /scratch/bistbs/Population_Genomic_Analysis/GONe/recode_map_gone.R \
        /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Addra/Dama_gazelle_Addra.map
```
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
plink --bfile /scratch/bistbs/Population_Genomic_Analysis/GONe/Dama_gazelle \
      --keep /scratch/bistbs/Population_Genomic_Analysis/GONe/Mohrr_gazelle.txt \
      --mac 1 \
      --recode \
      --not-chr h1tg000003l,h1tg000007l,h1tg000018l,h1tg000022l,h1tg000023l,h1tg000027l,h1tg000029l,h1tg000030l,h1tg000031l,h1tg000032l,h1tg000033l,h1tg000034l,h1tg000035l,h1tg000037l,h1tg000038l,h1tg000039l,h1tg000040l,h1tg000041l,h1tg000042l,h1tg000043l,h1tg000044l,h1tg000045l,h1tg000046l,h1tg000047l,h1tg000048l,h1tg000049l,h1tg000050l,h1tg000051l,h1tg000052l,h1tg000053l,h1tg000054l,h1tg000055l,h1tg000056l,h1tg000057l,h1tg000058l,h1tg000059l,h1tg000060l,h1tg000061l,h1tg000062l,h1tg000063l,h1tg000064l,h1tg000065l,h1tg000066l,h1tg000067l,h1tg000068l,h1tg000069l,h1tg000070l,h1tg000071l,h1tg000072l,h1tg000073l,h1tg000074l,h1tg000075l,h1tg000076l,h1tg000077l,h1tg000078l,h1tg000079l,h1tg000080l,h1tg000081l,h1tg000082l,h1tg000083l,h1tg000084l,h1tg000085l,h1tg000086l,h1tg000087l,h1tg000088l,h1tg000089l,h1tg000090l,h1tg000091l,h1tg000092l,h1tg000093l,h1tg000094l,h1tg000095l,h1tg000096l,h1tg000097l,h1tg000098l,h1tg000099l,h1tg000100l,h1tg000101l,h1tg000102l,h1tg000103l,h1tg000104l,h1tg000105l,h1tg000106l,h1tg000107l,h1tg000108l,h1tg000109l,h1tg000110l,h1tg000111l,h1tg000112l,h1tg000113l,h1tg000114l,h1tg000115l,h1tg000116l,h1tg000117l,h1tg000118l,h1tg000119l,h1tg000120l,h1tg000121l,h1tg000122l,h1tg000123l,h1tg000124l,h1tg000125l,h1tg000126l,h1tg000127l,h1tg000128l,h1tg000129l,h1tg000130l,h1tg000131l,h1tg000132l,h1tg000133l,h1tg000134l,h1tg000135l,h1tg000136l,h1tg000137l,h1tg000138l,h1tg000139l,h1tg000140l,h1tg000141l,h1tg000142l,h1tg000143l,h1tg000144l,h1tg000145l \
      --allow-extra-chr \
      --nonfounders \
      --debug \
      --out /scratch/bistbs/Population_Genomic_Analysis/GONe/output/Mohrr/Dama_gazelle_Mohrr


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

