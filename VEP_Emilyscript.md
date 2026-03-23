#### Finding deleterious mutation in the  Dama gazelle.
- Usually kept three species but here I am keeping two outgroups.
- Grant’s gazellehttps://trace.ncbi.nlm.nih.gov/Traces/?run=SRR6878810)  and Thomson's gazelle (https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR6894844)
#### Step 1. Download the short-read sequences of Grant's gazelle and Thomson's gazelle.
- Map these with BWA-MEM using the reference genome assembly of Dama gazelle.
- Since, I had already mapped each outgroup individually to the ref. genome and then already added RG and dedeuplicated, now I am mapping. 
```bash
#!/bin/bash -l
#SBATCH --job-name=MergeOutgroups_Dama
#SBATCH --time=50:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=80G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/merge_%j.out
#SBATCH --error=/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/merge_%j.err

# 1. Load Java (Using your specific version)
module load java-20

# 2. Define Paths from your provided script
PICARD_JAR="/shared/jezkovt_bistbs_shared/Guam_Rail/Guam_Rail_Analysis/Final_data_analysis/Alignment_BWAmem/Add_RG/picard.jar"
IN_DIR="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/ABBA_BABA/Haplotype_Caller/Five_samples/rmdup_Dama"
OUT_DIR="/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# 3. Define the two specific outgroup files to merge
# Grant's Gazelle and Thompson's Gazelle
BAM1="$IN_DIR/SRR6878810_sorted_RG_rmdup10X.bam"
BAM2="$IN_DIR/SRR6894844_rmdup_10X.bam"

echo "Starting merge for Outgroups..."

# 4. Run Picard MergeSamFiles
java -Xmx30g -jar "$PICARD_JAR" MergeSamFiles \
    I="$BAM1" \
    I="$BAM2" \
    O="$OUT_DIR/gazelle_outgroup_merged.bam" \
    TMP_DIR="$OUT_DIR" \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    USE_THREADING=true

echo "Merge Complete. Output file: $OUT_DIR/gazelle_outgroup_merged.bam"
```

#### Step 2. Consensus generation for Dama gazelle outgroups.

```bash
#!/bin/bash -l
#SBATCH --job-name=ANGSD_Dama
#SBATCH --time=54:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/angsd_%j.out
#SBATCH --error=/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/angsd_%j.err

# 1. Load the ANGSD module 
module load angsd

# 2. Define Input and Output
# Your merged outgroup BAM file
INPUT_BAM="/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/gazelle_outgroup_merged.bam"
# Where the ancestral fasta will be saved
OUT_PREFIX="/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/gazelle_outgroup_consensus"

echo "Starting ANGSD Consensus Generation for Dama Gazelle Outgroups..."

# 3. Run ANGSD
# -i: Input BAM
# -P: Number of threads (CPUs to use)
# -doFasta 2: Generates a consensus sequence (The 'Ancestral' map)
# -doCounts 1: Counts the bases at each site (Needed for doFasta)
angsd \
    -i "$INPUT_BAM" \
    -P 12 \
    -doFasta 2 \
    -doCounts 1 \
    -out "$OUT_PREFIX"

echo "Process Complete. Your ancestral file is: ${OUT_PREFIX}.fa.gz"
```
#### Step 3. Conversion of vcf file to .bed format
- We need to take every mutation SNP that exists in your five gazelles and put their address into a simple list.
```bash
# Define your paths
IN_VCF="/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/Dama_Gazelle_Final_Filtered_biallelic.recode.vcf"
OUT_DIR="/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups"
PLINK_EXE="/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/plink"

# Run the conversion
$PLINK_EXE --vcf "$IN_VCF" \
--make-bed \
--allow-extra-chr \
--out "$OUT_DIR/Dama_Final_Binary"
```
- After running the plink, run below command
---
awk 'BEGIN{OFS="\t"} {print $1, $4-1, $4}' Dama_Final_Binary.bim > Dama_SNPs_Coordinates.bed
---

#### Step 4. Extract ancestral alleles using bedtools
- Since, I already have .bed file and ANGSD consensus fasta.
```bash
module load bedtools-2.28
# 2. Define your paths
FASTA="/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/gazelle_outgroup_consensus.fa"
BED="/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/BEDfiles/Dama_SNPs_Coordinates.bed"
OUT_FILE="/home/bistbs/Dama_gazelle_VEP/Merging_Outgroups/ancestral_alleles.out"

# 3. Extract the ancestral base
# Note: Ensure the FASTA is unzipped first
bedtools getfasta -fi $FASTA -bed $BED -fo $OUT_FILE
```
#### Step 5. Reformat the Ancestral Alleles.
```bash
# 1. Extract just the DNA letters from the bedtools output
grep -v ">" ancestral_alleles.out > alleles_only.txt

# 2. Get the SNP IDs (Chromosome:Position) from your FIXED bed file
# Note: We use the 3rd column (actual position) to match PLINK's @:# format
awk '{print $1":"$3}' Dama_SNPs_Coordinates.bed > positions_only.txt

# 3. Paste them together to make the reference file
paste positions_only.txt alleles_only.txt > ancestral_alleles.txt
```

#### Step 6. Prepare the VCF IDs.
- PLINK needs the IDs in your vcf to match the IDs in your ancestral_alleles.txt
- A file named temp_dama.vcf where every SNP is named Chromosome: Position.
```bash
./plink --vcf Dama_Gazelle_Final_Filtered_biallelic.recode.vcf \
--set-missing-var-ids @:# \
--recode vcf \
--allow-extra-chr \
--out temp_dama
```

#### Step 7. Polarize (Big Switch)
- Now, tell PLINK to force the ancestral letter from your list to be the Reference(0) allele.
```bash
 1. Create a list of positions to keep
awk '{print $1}' ancestral_alleles.txt > ancestral_positions.txt

# 2. Run the polarization
./plink --vcf temp_dama.vcf \
--extract ancestral_positions.txt \
--a2-allele ancestral_alleles.txt 2 1 \
--recode vcf \
--allow-extra-chr \
--out temp_dama_polarized
```    

#### Step 8. Final Cleanup(Removing Mismatches)
- Remove the SNPs where the ancestral letter doesn't match the Dama gazelle letters (e.g., Ancestro has "A", but Dama only has "C/G")
# 1. Identify the mismatches from the log file
grep 'Warning' temp_dama_polarized.log | awk '{print $7}' | sed 's/.//;s/.$//' > mismatches.txt

# 2. Create the final, clean, polarized VCF
./plink --vcf temp_dama_polarized.vcf \
--exclude mismatches.txt \
--allow-extra-chr \
--export vcf-4.2 \
--out Dama_Gazelle_POLARIZED_Final

# 3. Clean up the temporary files to save space
rm temp_dama.vcf temp_dama_polarized.vcf alleles_only.txt positions_only.txt
```

#### Step 9. Variant Annotation using VEP or SnpEff.
- The final vcf has the ALT allele (the '1') is officially the Derived Mutation.
- This is exactly what VEP needs to tell you if the mutation is harmful.
