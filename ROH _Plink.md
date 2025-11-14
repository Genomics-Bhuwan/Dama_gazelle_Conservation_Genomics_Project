#### Calculate the ROH using Plink

###### Convert VCF to PLINK binary format (BED/BIM/FAM)
```bash
plink --vcf /scratch/bistbs/Population_Genomic_Analysis/ROH/Plink/Dama_gazelle_biallelic_snps_filtered.recode.vcf \
      --make-bed \
 --allow-extra-chr \
      --out /scratch/bistbs/Population_Genomic_Analysis/ROH/Plink/Dama_gazelle
```

###### Calculate the ROH using Plink
```bash
plink --bfile /scratch/bistbs/Population_Genomic_Analysis/ROH/Plink/Dama_gazelle \
      --homozyg \
      --homozyg-window-snp 50 \
      --homozyg-snp 50 \
      --homozyg-kb 500 \
      --homozyg-gap 1000 \
      --homozyg-density 50 \
      --homozyg-window-missing 5 \
      --homozyg-window-het 3 \
--allow-extra-chr \
      --out /scratch/bistbs/Population_Genomic_Analysis/ROH/Plink/Dama_gazelle_ROH
```

