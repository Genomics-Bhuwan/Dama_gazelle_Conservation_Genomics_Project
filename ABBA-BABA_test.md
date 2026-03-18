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
#!/bin/bash -l
#SBATCH --job-name=Dsuite_Gazelle
#SBATCH --time=54:00:00        # 24 hours should be sufficient for ~25,000 SNPs
#SBATCH --cpus-per-task=14     # Dsuite can utilize multiple threads
#SBATCH --mem=40G              # 20 GB is usually enough, but keeping buffer
#SBATCH --partition=batch
#SBATCH --output=Dsuite_%J.out
#SBATCH --error=Dsuite_%J.err  # Use this file to check if the job is running or errors occur
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# 1. Load required modules (if needed on your server)
# module load gcc/9.3.0  

# 2. Navigate to the working directory (where VCF and SETS.txt are located)
cd ~/Dama_gazelle_Abba_babba/Combine_Gentopyes_Joint_Genotyping/Joint_Genotyping/Variant_Filtration/biallelic/ABBA_BABA_test/

# 3. Print a message to confirm the job has started
echo "Dsuite analysis started at: $(date)"

# 4. Run Dsuite Dtrios command
./Dsuite/Build/Dsuite Dtrios \
-n Gazelle_Total \
Dama_gazelles_final_biallelic.vcf.gz \
SETS.txt

# 5. Print a message when the job finishes
echo "Dsuite analysis finished at: $(date)"
```
#### Step 3. Dinvestigate
```bash
#!/bin/bash -l
#SBATCH --job-name=Dinvestigate_Gazelle
#SBATCH --time=48:00:00        # Sliding window analysis may take longer time
#SBATCH --cpus-per-task=12     # Use 12 CPUs
#SBATCH --mem=40G              # Allocate 40 GB memory
#SBATCH --partition=batch
#SBATCH --output=Dinvestigate_%J.out
#SBATCH --error=Dinvestigate_%J.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu

# 1. Move to the main working directory
cd ~/Dama_gazelle_Abba_babba/Combine_Gentopyes_Joint_Genotyping/Joint_Genotyping/Variant_Filtration/biallelic/ABBA_BABA_test/

# 2. Print a message indicating the job has started
echo "Dinvestigate sliding window analysis started at: $(date)"

# 3. Run the Dinvestigate command
# -w 50000,10000 means a window size of 50 kb and a step size of 10 kb
./Dsuite/Build/Dsuite Dinvestigate \
-w 50000,10000 \
-n Gazelle_Local \
Dama_gazelles_final_biallelic.vcf.gz \
SETS.txt \
test_trios.txt

# 4. Print a message when the analysis is finished
echo "Dinvestigate analysis finished at: $(date)"

```
