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


