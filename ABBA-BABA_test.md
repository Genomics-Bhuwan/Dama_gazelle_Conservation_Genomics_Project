#### Hybridization among the two sub-species
- Firstly, We used five samples of Addra and Mhorr gazelle from this study including the Nanger granti as an outgroup NCBI Accession:SRR6878810: SRR6878810: https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR6878810).
- Link: https://github.com/millanek/Dsuite
- https://github.com/simonhmartin/tutorials/blob/master/ABBA_BABA_windows/README.md   (Recommended to use for analysis and for plotting)
- Fast calculation of Patterson's D (ABBA-BABA) and the f4-ratio statistics across many populations/species
- Input files: A VCF file: can contain multiallelic loci and indels but only biallelic SNPs will be used.

  ---
 - This excess can be expressed in terms of the D statistic, which ranges from -1 to 1,.
 -  and should equal 0 under the null hypothesis of no introgression.
 -  D > 1 indicates possible introgression between P3 and P2 (or other factors that would result in a deviation from a strict bifurcating species history).
  ---
- Population/species map SETS.txt
- contains a .txt file with one individual per row and a tab separating the individual's name from the name of the species/populations it belongs to
Ind1    Species1
Ind2    Species1
Ind3    Species2
Ind4    Species2
Ind5    Species3
Ind6    Outgroup
Ind7    Outgroup

#### Step 1. Calculate the Depth of the samples. Downsample them
#### Step 2. Run Dsuite
```bash
/home/bistbs/Dama_gazelle_Abba_babba/Combine_Gentopyes_Joint_Genotyping/Joint_Genotyping/Variant_Filtration/biallelic/ABBA_BABA_test/Dsuite/Build/Dsuite Dtrios \
/home/bistbs/Dama_gazelle_Abba_babba/Combine_Gentopyes_Joint_Genotyping/Joint_Genotyping/Variant_Filtration/biallelic/Dama_gazelle_genotypes_pass_biallelic.vcf.gz \
/home/bistbs/Dama_gazelle_Abba_babba/Combine_Gentopyes_Joint_Genotyping/Joint_Genotyping/Variant_Filtration/biallelic/ABBA_BABA_test/SETS.txt
```
