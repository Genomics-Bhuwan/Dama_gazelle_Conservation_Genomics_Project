# Mitogenome Assembly & Phylogenomic Tutorial
  - Outgroups used is:
 ---
  Gazella dorcas – Dorcas gazelle		
     a. isolate 1 subspecies osiris	Mali- JN632637*
     b. isolate 2 subspecies pelzelnii AWWP 5722, Somalia- JN632638*
    - Downloaded from the paper: https://www.sciencedirect.com/science/article/pii/S1631069111002800
- Used sample:
Nanger soemmerringii isolate AWWP mitochondrion, complete genome
GenBank: JN632667.1
---

Heterozygosity is a simple and informative statistic that can be obtained by analyzing whole-genome data. You can calculate the average heterozygosity of an individual or assess the local heterozygosity throughout the genome to answer more specific questions.

There are many methods available to achieve this. In this tutorial, we will use two strategies. To obtain the average heterozygosity, we will manipulate the VCF file using bcftools and vcftools. To explore heterozygosity variation across the genome, we will use ANGSD and other associated tools.

### Mitogenome Assembly using MitoZ

First we want to subsample the whole genome reads for each individual we have. For the purposes of practicing we are only going to use our two resequenced individuals for the mitogenome assembly but for the whole phylogeny later we are going to have different individuals included

```
mkdir mitogenome_phylogeny/
mkdir mitogenome_phylogeny/sub_reads
cd /mitogenome_phylogeny/
```
##### Step 1. randomly subsample 20% of read pairs
- This will randomly sub-sample the 20% of the read pairs from reads1 or forward read cause the mtDNA retention is optimal enough from 20%.
```bash
#!/bin/bash

INPUT_DIR="/localscratch/bistbs/mitogenome_phylogeny/sub_reads"
OUTPUT_DIR="/localscratch/bistbs/mitogenome_phylogeny/sub_reads_sub20"
mkdir -p "$OUTPUT_DIR"

SAMPLES=("SRR17129394" "SRR17134085" "SRR17134086" "SRR17134087" "SRR17134088")
SEQTK="/localscratch/bistbs/mitogenome_phylogeny/seqtk/seqtk"
PIGZ="/localscratch/bistbs/mitogenome_phylogeny/pigz-2.8/pigz"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Subsampling read1 for $SAMPLE..."
    $SEQTK sample -s 100 "$INPUT_DIR/${SAMPLE}_1_val_1.fq" 0.2 | $PIGZ -p 20 > "$OUTPUT_DIR/${SAMPLE}_1.sub.fq.gz"
    echo "$SAMPLE read1 subsampled."
done

```

##### Step 2: Merge subsampled read1 with original read2 (also .fq, not .fq.gz)
- Interleave is like merging the 20% sub-sampled forward reads with the complete F2 reads.
```bash
   for SAMPLE in "${SAMPLES[@]}"; do
    echo "Merging/interleaving subsampled R1 with full R2 for $SAMPLE..."

    /localscratch/bistbs/mitogenome_phylogeny/seqtk/seqtk mergepe \
        /localscratch/bistbs/mitogenome_phylogeny/sub_reads_sub20/${SAMPLE}_1.sub.fq.gz \
        /localscratch/bistbs/mitogenome_phylogeny/sub_reads/${SAMPLE}_2_val_2.fq \
        > /localscratch/bistbs/mitogenome_phylogeny/sub_reads_sub20/${SAMPLE}.interleave.fq

    echo "Done: ${SAMPLE}.interleave.fq"
done

```

#### Step 3: Deinterleave using your local script
- Since, the 20% will only align with the resepctive region in the reverse reads.
- Since, paste command is faster than awk. Therefore, I am using awk.
```bash
  for f in *.interleave.fq; do
    base=${f%.interleave.fq}

    echo "Deinterleaving $base ..."

    paste - - - - - - - - < "$f" | \
        tee >(cut -f1-4 | tr '\t' '\n' | gzip > ${base}_1.sub20.fq.gz) \
            >(cut -f5-8 | tr '\t' '\n' | gzip > ${base}_2.sub20.fq.gz) >/dev/null

    echo "$base done."
done

```
##### Step 4. Our main target is to run MitoZ as a loop over all of our samples via their subreads.
---
- MitoZ generates the de novo mitogenome assembly and annotates teh resulting mitogenome.
---
- This is the step where I had alredy apptainer/singularity container.
- I am using MitoZ inside a Apptainer container(A singularity variant).
- I tried below command.  
Command you tried:
```bash
singularity exec MitoZ_v3.6.sif bash
```
- Now, Open the container in the interactive session and see the conda environment inside the container.
- MitoZ requires Python and some bioinformatics tools (MEGAHIT, SPAdes).
- Inside the container, these are installed in a Conda environment called mitoz3.6.

##### Check Conda and environments:
```bash
which conda
```
```bash
conda info --envs
```
- activate the environment:
```bash
source /app/anaconda/bin/activate mitoz3.6
```
- (mitoz3.6) appears in the prompt, meaning the environment is active.
- Now you can access MitoZ, MEGAHIT, SPAdes, etc.
```bash
 which megahit
```
```bash
which spades.py
```
---
- Running MitoZ for one sample. But will be replicating for rest of the four sample for this project.
- --
```bash
# Output directory
OUTDIR=/scratch/bistbs/Population_Genomic_Analysis/mitogenome_haplotype_phylogeny/MitoZ_output

# Correct sample
SAMPLES=(SRR17129394)

# ===============================================================
# Loop over each sample and run MitoZ with full paths to FASTQs
# ===============================================================
for SAMPLE in "${SAMPLES[@]}"
do
    echo "--------------------------------------------------"
    echo "Running MitoZ for $SAMPLE..."
    echo "Output will be saved in: ${OUTDIR}/${SAMPLE}"
    
    # Run MitoZ all module
    python /app/anaconda/envs/mitoz3.6/lib/python3.8/site-packages/mitoz/MitoZ.py all \
        --outprefix $SAMPLE \
        --thread_number 24 \
        --fq1 /scratch/bistbs/Population_Genomic_Analysis/mitogenome_haplotype_phylogeny/${SAMPLE}_1.sub20.fq.gz \
        --fq2 /scratch/bistbs/Population_Genomic_Analysis/mitogenome_haplotype_phylogeny/${SAMPLE}_2.sub20.fq.gz \
        --clade Chordata \
        --genetic_code 2 \
        --requiring_taxa Mammalia \
        --species_name 'Dama gazelle' \
        --insert_size 350 \
        --skip_filter \
        --workdir ${OUTDIR}/${SAMPLE}

    echo "$SAMPLE done."
done

```
##### Run the loop for rest of the samples.
- Step 1: Enter the Apptainer container
- Launch the container that has MitoZ installed.
- This gives you a controlled environment with all dependencies.
 ```bash
singularity exec MitoZ_v3.6.sif bash
```
##### Step 2: Check and activate conda environment
- Locate conda in the container
```bash
which conda
```
##### Activate the MitoZ environment (contains megahit, spades, python, etc.)
```bash
source /app/anaconda/bin/activate mitoz3.6
```
- Prompt should show: (mitoz3.6) Apptainer>
- Verify main tools are available
```bash
which megahit     # /app/anaconda/envs/mitoz3.6/bin/megahit
which spades.py   # /app/anaconda/envs/mitoz3.6/bin/spades.py
which python      # /app/anaconda/envs/mitoz3.6/bin/python
```

##### Step 3: Organize output
##### Define a base directory for all sample outputs
```bash
OUTDIR="/scratch/bistbs/Population_Genomic_Analysis/mitogenome_haplotype_phylogeny/MitoZ_Output_Final"
```
##### Make sure the directory exists
```bash
mkdir -p "$OUTDIR"
```
##### Step 4: Prepare sample list
- List all the SRR samples you want to process
```bash
SAMPLES=(SRR17134085 SRR17134086 SRR17134087 SRR17134088)

# ===============================================================
# Step 5: Loop over each sample and run MitoZ
# ===============================================================
for SAMPLE in "${SAMPLES[@]}"
do
    echo "--------------------------------------------------"
    echo "Running MitoZ for $SAMPLE..."
    echo "Output will be saved in: ${OUTDIR}/${SAMPLE}"
    
    # Run MitoZ all module
    python /app/anaconda/envs/mitoz3.6/lib/python3.8/site-packages/mitoz/MitoZ.py all \
        --outprefix $SAMPLE \
        --thread_number 24 \
        --fq1 ./${SAMPLE}_1.sub20.fq.gz \
        --fq2 ./${SAMPLE}_2.sub20.fq.gz \
        --clade Chordata \
        --genetic_code 2 \
        --requiring_taxa Mammalia \
        --species_name 'Dama gazella' \
        --fastq_read_length 125 \
        --insert_size 350 \
        --skip_filter \
        --workdir ${OUTDIR}/${SAMPLE}

    echo "$SAMPLE done."
done
```
```bash
# ===============================================================
# Step 6: Check final outputs
# ===============================================================
# After completion, each sample will have its own folder:
# /scratch/bistbs/Population_Genomic_Analysis/mitogenome_haplotype_phylogeny/MitoZ_Output_Final/SRR17134085/
# Each folder contains:
# - assembled mitochondrial contigs (FASTA)
# - annotation files
# - log files
# - intermediate files from Megahit/Spades

# ===============================================================
# Notes / Comments
# ===============================================================
# 1. We activate the correct conda environment inside Apptainer to ensure MitoZ, megahit, spades, and Python versions are correct.
# 2. Output directory structure keeps all samples separate to avoid overwriting.
# 3. Thread number (--thread_number) can be adjusted depending on available CPUs.
# 4. skip_filter can be used if you don't want MitoZ to filter contigs by taxonomy.
# 5. You can parallelize this loop if running on HPC with job arrays for faster processing.
```
###### Step 5. Finding a mitogenone in an assembly using BLASTN
- I downloaded the mitochondrial reference genome assembly of Addra gazelle from NCBI(https://www.ncbi.nlm.nih.gov/nuccore/NC_020724.1) in the folder given below.

```bash
cd /scratch/bistbs/Population_Genomic_Analysis/mitogenome_haplotype_phylogeny/mitogenome

```


##### Step 6: Make a BLAST database from the mitochondrial reference genome
```bash
module load blast-2.13.0+
makeblastdb -in /scratch/bistbs/Population_Genomic_Analysis/mitogenome_haplotype_phylogeny/Nanger_dama_mitochondrial_reference_genome.fasta -dbtype nucl
```
##### Step 7: BLAST the  genome assembly against the mitochondrial reference
```bash
blastn -db /scratch/bistbs/Population_Genomic_Analysis/mitogenome_haplotype_phylogeny/Nanger_dama_mitochondrial_reference_genome.fasta \
-query /scratch/bistbs/Population_Genomic_Analysis/mitogenome_haplotype_phylogeny/Dama_gazelle_hifiasm-ULONT_primary.fasta \
-outfmt 6 -max_target_seqs 1 -max_hsps 1 -num_threads 12 \
| sort -V -k12,12r \
| head -n1 \
| cut -f1 > Dama_gazelle_mitoscaf.txt
```
##### Step 8: Extract the mitochondrial scaffold
```bash
seqtk subseq /scratch/bistbs/Population_Genomic_Analysis/mitogenome_haplotype_phylogeny/Dama_gazelle_hifiasm-ULONT_primary.fasta Dama_gazelle_mitoscaf.txt > Dama_gazelle_mitogenome.fasta
```

##### 9. Mitogenome Phylogenomics
---
- Now, I have a mitogenome for our samples.
- I will do whole genome phylogeny using MAFFT and IQTREE.
- Typically, We use the annotations to only include the 12 protein-coding genes(PCGs).
- But, here we are only going to align the whole mitogenome.
- How to concatenate the PCGs can be found here https://github.com/rtfcoimbra/Coimbra-et-al-2021_CurrBiol/blob/main/mitogenomes_workflow.txt
---

##### Step 1. For whole genome phylogeny.
- *Step A*: Organize your mitogenomes
- Make a new working folder for alignment/phylogeny:
```bash
mkdir ~/mitogenomes/Phylogeny
cd ~/mitogenomes/Phylogeny
```
- Copy your MitoZ output FASTAs:
```bash
cp /scratch/bistbs/Population_Genomic_Analysis/mitogenome_haplotype_phylogeny/MitoZ_Output_Final/*.fasta to phylogeny folder.
```
- *Step B*: Align the whole mitogenome (optional) If you want whole mitogenome phylogeny:
- You have to combine all the mitogenomes in one file using cat command
- combine all mitogenomes in one file

```bash
- Renamed the samples with the sample names.
mv SRR17129394_SRR17129394.megahit.mitogenome.fa_mitoscaf.fa.gbf.fasta SRR17129394.fasta
mv SRR17134085_SRR17134085.megahit.mitogenome.fa_mitoscaf.fa.gbf.fasta SRR17134085.fasta 
mv SRR17134086_SRR17134086.megahit.mitogenome.fa_mitoscaf.fa.gbf.fasta SRR17134086.fasta
mv SRR17134087_SRR17134087.megahit.mitogenome.fa_mitoscaf.fa.gbf.fasta SRR17134087.fasta
mv SRR17134088_SRR17134088.megahit.mitogenome.fa_mitoscaf.fa.gbf.fasta SRR17134088.fasta

- Combine all fastas into one fasta
- Downloaded the Nanger granti as an outgroup for these individuals.
- Download it using wget
```bash
wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_020725.1&db=nuccore&report=fasta&extrafeat=null&conwithfeat=on&hide-cdd=on" -O Nanger_granti_NC_020725.1.fasta

cat *.fasta > Dama_whole_mitogenomes.fa
```
- **Align using MAFFT* 
```bash
mafft --auto Dama_whole_mitogenomes.fa > Dama_whole_mitogenomes_aligned.fa
```
##### Run whole genome phylogenetic tree
```bash
#!/bin/bash -l
#SBATCH --time=150:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=90G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --job-name=Whole_genome_mtDNA

### Load modules if needed ###
module load gcc/10.2.0   # or your system’s default
module load mafft         # optional if using MAFFT for alignment

### Change to your working directory ###
cd /scratch/bistbs/Population_Genomic_Analysis/mitogenome_haplotype_phylogeny/mitogenome/Phylogeny/Whole_genome_mtDNA_phylogeny/

### Download IQ-TREE if not already present ###
wget https://github.com/iqtree/iqtree2/releases/download/v2.2.0/iqtree-2.2.0-Linux.tar.gz
tar -zxvf iqtree-2.2.0-Linux.tar.gz
export PATH=$(pwd)/iqtree-2.2.0-Linux/bin:$PATH

### Run IQ-TREE on the whole mitogenome alignment ###

iqtree2 -s Dama_with_outgroup_aligned.fasta -m TESTNEW -bb 1000 -alrt 1000 -nt 24 -pre Dama_whole_mtDNA 
```

##### Phylogeny using Protein-Coding Genes(Recommended Method).

- Extract protein-coding genes (PCGs) (recommended)
- Mitochondrial phylogenies are usually more robust using PCGs.
- Create a folder for PCGs:

##### Step 1: Prepare working folder for PCGs
```bash
mkdir -p mtdna_cds_work
cp *.cds.fasta mtdna_cds_work/
cd mtdna_cds_work
```

##### Step 2: Rename headers for simplicity
- Goal: make each header Sample.Gene (avoiding long messy headers).
```bash
mkdir -p genes_split

for f in *.cds.fasta; do
    sample=$(basename "$f" .cds.fasta)
    awk -v s="$sample" '
    /^>/ {
        # Extract gene name from header (works for both formats)
        if ($0 ~ /\[gene=([A-Za-z0-9]+)\]/) {
            match($0, /\[gene=([A-Za-z0-9]+)\]/, a)
            gene = a[1]
        } else if ($0 ~ /\.([A-Za-z0-9_-]+)$/) {
            match($0, /\.([A-Za-z0-9_-]+)$/, a)
            gene = a[1]
        } else {
            gene = "UNKNOWN"
        }
        print ">" s "." gene
        next
    }
    {print}
    ' "$f" > tmp && mv tmp "$f"
done
```
- This works for both Nanger_granti_CDS.fasta and SRR*.cds.fasta.
- Ensures headers are sample name + gene only.

##### Step 3: Split all genes into per-gene FASTA files
```bash
pcgs=(ATP6 ATP8 COX1 COX2 COX3 CYTB ND1 ND2 ND3 ND4 ND4L ND5 ND6)

for g in "${pcgs[@]}"; do
    grep -h -A1 --no-group-separator "\.${g}$" *.cds.fasta > genes_split/${g}.fasta
done
```

##### Creates genes_split/ATP6.fasta, genes_split/ND1.fasta, etc.
##### Step 4: Align each gene separately
```bash
cd genes_split

for g in "${pcgs[@]}"; do
    mafft --auto ${g}.fasta > ${g}_aligned.fasta
done
```

- Produces files like ATP6_aligned.fasta, ND1_aligned.fasta etc.

##### Step 5: Concatenate aligned genes into a supermatrix
- Make sure all sequences are in same order of samples across genes.
- Simple concatenation:
```bash
cat CYTB_aligned.fasta ND6_aligned.fasta ND5_aligned.fasta ND4_aligned.fasta ND4L_aligned.fasta ND3_aligned.fasta COX3_aligned.fasta ATP6_aligned.fasta ATP8_aligned.fasta COX2_aligned.fasta COX1_aligned.fasta ND2_aligned.fasta ND1_aligned.fasta > Dama_gazelle_PCGs_aligned.fasta
```

##### Step 6: Optional: generate partition file for IQ-TREE or RAxML:

```bash
# Define genes in the order they appear in your concatenated file
genes=(CYTB ND6 ND5 ND4 ND4L ND3 COX3 ATP6 ATP8 COX2 COX1 ND2 ND1)

# Initialize variables
start=1
> partitions.txt  # empty file to store partitions

# Loop through each gene
for gene in "${genes[@]}"; do
    # get length of the gene alignment (number of characters excluding headers)
    length=$(grep -v ">" "${gene}_aligned.fasta" | tr -d '\n' | wc -c)
    
    # calculate end position
    end=$((start + length - 1))
    
    # write to partition file in IQ-TREE format
    echo "DNA, $gene = $start-$end" >> partitions.txt
    
    # update start for next gene
    start=$((end + 1))
done

echo "Partition file created: partitions.txt"
```


# RUN IQTREE2 for phylogenetics for 13 PCGs.
```bash
iqtree2 -s Dama_gazelle_PCGs_aligned_with_Nanger.fasta -q partitions.txt -m MFP -bb 1000 -alrt 1000
```


Label your samples clearly (Dama_1 … Dama_5).
