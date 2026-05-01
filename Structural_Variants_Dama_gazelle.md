#### Structural variants analysis using smoove pipeline
- https://github.com/brentp/smoove

#### Step 1 a. Download the smoove pipeline using wget
```bash
# Download the smoove binary
wget https://github.com/brentp/smoove/releases/download/v0.2.8/smoove

# Make it executable
chmod +x smoove

# Run the help command to see if it works
./smoove -h

```
#### Step 1 b. Install the dependencies
  ```bash
  mamba install -c bioconda -c conda-forge htslib gsort lumpy-sv samtools svtyper mosdepth
  ```

  #### Step 2. Running for small sample cohort for structural variants
```bash
  ./smoove call -x \
    --name Dama_Gazelle_Cohort \
    --fasta /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Structural_Variants_SaVOR/Dama_gazelle_hifiasm-ULONT_primary.fasta \
    -p 4 \
    --genotype \
    /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Structural_Variants_SaVOR/svArcher_Standalone/results/Dama_Gazelle/bams/*.bam
```
#### Step 3. Run the duphold to add the depth annotations(DHFFC).
# Optional: Run duphold to add depth annotations (DHFFC)
# Note: You would usually add -d to the initial 'call' to do this automatically
# If you didn't, you can run it now:

#### Step 3 a. Annotation included
```bash
./smoove duphold -x \
    --name Dama_Gazelle_Annotated \
    --vcf Dama_Gazelle_Cohort-smoove.genotyped.vcf.gz \
    --fasta /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Structural_Variants_SaVOR/Dama_gazelle_hifiasm-ULONT_primary.fasta \
    /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Structural_Variants_SaVOR/svArcher_Standalone/results/Dama_Gazelle/bams/*.bam
```

#### Step 3 b. Run annotate with a GFF file for the Dama gazelle.
```bash
./smoove annotate --gff your_gazelle.gff3.gz Dama_Gazelle_Annotated.vcf.gz | bgzip -c > Dama_Gazelle_Final.vcf.gz
```

#### Step 4. Variant filtration for Smoove for high-quality genotypes and real depth changes using LUMPY
- Filter Smoove for high-quality genotypes and real depth changes
- MSHQ > 3: High quality heterozygotes
- DHFFC: deletions with DHFFC < 0.7 and duplications with DHFFC > 1.25
- These are stringently requested by the smoove pipeline.
```bash
bcftools view -i 'MSHQ > 3 && (SVTYPE == "DEL" && DHFFC < 0.7 || SVTYPE == "DUP" && DHFFC > 1.3 || SVTYPE == "INV")' \
    Dama_Gazelle_Annotated.vcf.gz | \
bcftools view -i 'QUAL >= 100' -O z -o Dama_Gazelle_Smoove_Clean.vcf.gz

# Always index after filtering
tabix -p vcf Dama_Gazelle_Smoove_Clean.vcf.gz
```



#### Step 5. Running Manta and Delly for structural variant calling
- DELLY- https://github.com/dellytools/delly
- singularity pull delly.sif docker://dellytools/delly:latest    : Use this for installing the delly.

#### Step 5 a. Preparation of the files.
- Create the exclusion file cause you have unplaced scaffolds and variable sex chromosomes.
- Avoid wasting time on hundreds of h1tg scaffolds.

```bash
 # Set Paths
REF="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Structural_Variants_SaVOR/Dama_gazelle_hifiasm-ULONT_primary.fasta"
BAM_DIR="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Structural_Variants_SaVOR/svArcher_Standalone/results/Dama_Gazelle/bams"
DELLY_IMG="./delly.sif"
# Create exclusion list for anything not Chromosomes 1-17
grep ">" $REF | sed 's/>//g' | grep -vE "^([1-9]|1[0-7])$" > exclude.txt
```

#### Step 5 b: SV Discovery for each of the individual sample.
- Mine is Illumina paired-end short read sequence. Use -y ont for oxford nanoepore sequencing and -y pb for PacBio.
- Run variant discovery for the five samples.
```bash
SAMPLES=("SRR17129394" "SRR17134085" "SRR17134086" "SRR17134087" "SRR17134088")

for s in "${SAMPLES[@]}"; do
    echo "Short-Read Discovery: $s"
    singularity exec --bind /shared:/shared $DELLY_IMG delly call \
        -g $REF \
        -x exclude.txt \
        -o "${s}.bcf" \
        "${BAM_DIR}/${s}_mapped_sorted_RG_rmdup.bam"
done
```
#### Step 5 c.: Merge Sites
- We are now merging the findings from all 5 samples into a single unified list of candidate SVs.
```bash
singularity exec --bind /shared:/shared $DELLY_IMG delly merge -o sites.bcf \
    SRR17129394.bcf SRR17134085.bcf SRR17134086.bcf SRR17134087.bcf SRR17134088.bcf
```
#### Step 5 d.: Genotyping (Per Sample)
- We now go back to the original BAM files to see  which sample has which variant from the unified sites.bcf list.
```bash
for s in "${SAMPLES[@]}"; do
    echo "Genotyping: $s"
    singularity exec --bind /shared:/shared $DELLY_IMG delly call \
        -g $REF \
        -v sites.bcf \
        -x exclude.txt \
        -o "${s}.geno.bcf" \
        "${BAM_DIR}/${s}_mapped_sorted_RG_rmdup.bam"
done
```

#### Step 6 d.: Final Merge and Chromosome Filtering
- Finally, we merge the individual genotyped fiels and use bcftools to stringently keep only chromosomes 1-17.
```bash
# Merge all 5 samples into one multi-sample BCF
# Merge genotyped files
singularity exec --bind /shared:/shared $DELLY_IMG bcftools merge -m id -O b -o dama_merged.bcf \
    SRR17129394.geno.bcf SRR17134085.geno.bcf SRR17134086.geno.bcf SRR17134087.geno.bcf SRR17134088.geno.bcf

# Index the result
singularity exec --bind /shared:/shared $DELLY_IMG bcftools index dama_merged.bcf

# Filter for Chrs 1-17 and output VCF
singularity exec --bind /shared:/shared $DELLY_IMG bcftools view dama_merged.bcf \
    -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17 \
    -o dama_gazelle_final_shortread_1-17.vcf

```
- Summary of the above files created
- exclude.txt lists everuthing that you want to ignore.
- sites.bcf: master list of SVs found across the cohort of 5.
- *.geno.bcf: Individual files containing the presence/absence of SVs.
  
#### Step 6 e.: Final filtration
- Use bcftools and only keep the varaints that pass the quality. Rest eliminate.
- After filtration, look at the genotype field in the vcf file and see
- 0/0: Homozygous Reference (No SV).
- 0/1: Heterozygous SV.
- 1/1: Homozygous Alternative SV.
