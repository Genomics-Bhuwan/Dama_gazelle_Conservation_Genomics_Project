#### Dxy and Pxy within and between the population of Addra and Mhorr gazelle.
- Link: GITHUB: https://github.com/ksamuk/pixy
- https://pixy.readthedocs.io/en/latest/; 
- Paper: https://doi.org/10.1111/1755-0998.13326

#### Step 1. Get the monomorphic sites or invariants site as well as variants site.
```bash
#!/bin/bash -l
#SBATCH --time=75:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --partition=batch
#SBATCH --job-name=Monomorphic_site
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bistbs@miamioh.edu

# Move to working directory
cd /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Pixy
module load gatk-4.1.2.0 

# Define your paths for clarity
REFERENCE="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Dama_gazelle_hifiasm-ULONT_primary.fasta"
INPUT_GVCF="all_samples_combined.g.vcf.gz"
OUTPUT_VCF="all_samples_allsites_genotyped.vcf.gz"

# Run GATK GenotypeGVCFs
# We increase memory to 90G given the T2T assembly size
gatk --java-options "-Xmx90g" GenotypeGVCFs \
  -R $REFERENCE \
  -V $INPUT_GVCF \
  -O $OUTPUT_VCF \
  --include-non-variant-sites true
```


#### Step 2. Go for variant filtration using vcf-tools
```bash
vcftools --gzvcf all_samples_allsites_genotyped.vcf.gz \
  --remove-indels \
  --max-missing 0.8 \
  --min-meanDP 10 \
  --max-meanDP 100 \
  --recode --stdout | gzip -c > filtered_allsites.vcf.gz
```

#### Step 3. Get the names of the samples of Addra and Mhorr
- Make a file and seperate them.
- Calcualte the pi dxy fst among these six pair of population using below code or revise it if needed.
  ---
Specifically, the dXY estimator can be used to estimate the absolute divergence between populations, along with Fst. I suggest we test the pairwise differences between each of the 3 addra and 2 mhorr gazelles (6x pairwise comparisons) as well as within each subspecies (between the 3 Addra and between the two Mhorr gazelles).
  ---


#### Step 4. Rather than running pixy on whole data file one at a time,
- Run it based on chromosomes. 

```bash
pixy --stats pi dxy fst \

  -vcf chr1_filtered.vcf.gz \

  --populations individual_pops.txt \

  --window_size 50000 \

  --n_cores 4 \

  --output_prefix chr1_results
```
