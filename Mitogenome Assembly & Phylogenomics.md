# Mitogenome Assembly & Phylogenomic Tutorial

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
```
#!/bin/bash

# List of samples
for sample in SRR17129394 SRR17134085 SRR17134086 SRR17134087 SRR17134088
do
    # Run MitoZ
    python3 /localscratch/bistbs/mitogenome_phylogeny/sub_reads_sub20/mitoz-3.6/mitoz/MitoZ.py all \
        --genetic_code 2 \
        --clade 'Chordata' \
        --outprefix ${sample} \
        --thread_number 24 \
        --fastq1 ../../sub_reads_sub20/sub_reads/${sample}_1.sub.fq.gz \
        --fastq2 ../../sub_reads_sub20/sub_reads/${sample}_2.sub.fq.gz \
        --fastq_quality_shift \
        --fastq_read_length 125 \
        --duplication \
        --insert_size 350 \
        --run_mode 2 \
        --filter_taxa_method 1 \
        --requiring_taxa 'Mammalia' \
        --species_name 'Dama gazella' &> ${sample}.mitoz.log
done


```

##### Step 5. Finding a mitogenone in an assembly using BLASTN
- I downloaded the mitochondrial reference genome assembly of Addra gazelle from NCBI(https://www.ncbi.nlm.nih.gov/nuccore/NC_020724.1)
```bash
cp /data/genomics/workshops/smsc_2024/mitogenomes/*.fa . mitogenomes/
```


##### Step 6: Make a BLAST database from the mitochondrial reference genome
```bash
makeblastdb -in /scratch/bistbs/Population_Genomic_Analysis/mitogenome_haplotype_phylogeny/Nanger_dama_mitochondrial_reference_genome.fasta -dbtype nucl
```
##### Step 7: BLAST the nuclear genome against the mitochondrial reference
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
- Now, I have a mitogenome for our samples.
- I will do whole genome phylogeny using MAFFT and IQTREE.
- Typically, We use the annotations to only include the 12 protein-coding genes(PCGs).
- But, here we are only going to align the whole mitogenome.
- How to concatenate the PCGs can be found here https://github.com/rtfcoimbra/Coimbra-et-al-2021_CurrBiol/blob/main/mitogenomes_workflow.txt

##### 10.A Combine all the mitogenomes
- Assuming your mitochondrial fasta files are named something like *.fa:
```bash
cat *.fa > Dama_gazelle_mitogenomes.fasta
```

##### 10.B Align sequences with MAFFT
##### --adjustdirection ensures reverse-complement sequences are corrected.
```bash
mafft --thread 10 --adjustdirection Dama_gazelle_mitogenomes.fasta > Dama_gazelle_mitogenomes.aln
```


##### 10.C Build phylogenetic tree with IQ-TREE
-B 1000 → ultrafast bootstrap with 1000 replicates
--pre Dama_gazelle_mito → output prefix for all IQ-TREE results
```bash
iqtree -nt 10 --pre Dama_gazelle_mito -s Dama_gazelle_mitogenomes.aln -B 1000
```


##################################################################################################################################################
##################################################################################################################################################
### This will create a partition.nexus file for each gene type and get the Phylogenetic tree for all the five individuals.
#!/bin/bash
#============================================#
# Mitogenome Assembly & Phylogenomics Workflow
# Production-Ready Version
# - Subsample reads
# - Run MitoZ assembly
# - Extract PCGs
# - Align sequences
# - Build phylogenetic tree with IQ-TREE
#============================================#

#---------------Directories------------------#
INPUT_DIR="/localscratch/bistbs/mitogenome_phylogeny/sub_reads"
OUTPUT_DIR="/localscratch/bistbs/mitogenome_phylogeny/sub_reads_sub20"
MITOZ_OUT="/localscratch/bistbs/mitogenome_phylogeny/mitoz_out"
PCG_OUT="/localscratch/bistbs/mitogenome_phylogeny/mitoz_PCGs"
ALIGN_OUT="/localscratch/bistbs/mitogenome_phylogeny/alignment"
TREE_OUT="/localscratch/bistbs/mitogenome_phylogeny/tree"

# Create directories if not exist
mkdir -p $OUTPUT_DIR $MITOZ_OUT $PCG_OUT $ALIGN_OUT $TREE_OUT

#---------------Samples----------------------#
samples=("SRR17129394" "SRR17134085" "SRR17134086" "SRR17134087" "SRR17134088")

#---------------Subsample Reads---------------#
echo "### Subsampling reads (20% read1) ###"
for sample in "${samples[@]}"; do
    echo "Processing $sample..."

    # Subsample 20% of read1
    seqtk sample -s 100 $INPUT_DIR/${sample}_1_val_1.fq.gz 0.2 | pigz -p 8 > $OUTPUT_DIR/${sample}_1.sub.fq.gz

    # Merge subsampled read1 with original read2
    seqtk mergepe $OUTPUT_DIR/${sample}_1.sub.fq.gz $INPUT_DIR/${sample}_2_val_2.fq.gz > $OUTPUT_DIR/${sample}.interleave.fq

    # Deinterleave to get final read1 and read2
    sh /shared/jezkovt_bistbs_shared/scripts/deinterleave_fastq.sh < $OUTPUT_DIR/${sample}.interleave.fq \
        $OUTPUT_DIR/${sample}_1.sub20.fq.gz $OUTPUT_DIR/${sample}_2.sub20.fq.gz compress

    echo "$sample subsampling done."
done

#---------------Run MitoZ Assembly---------------#
echo "### Running MitoZ Assembly ###"
for sample in "${samples[@]}"; do
    echo "Assembling mitogenome for $sample..."

    MitoZ.py all \
        --genetic_code 2 \
        --clade 'Chordata' \
        --outprefix $MITOZ_OUT/${sample} \
        --thread_number 24 \
        --fastq1 $OUTPUT_DIR/${sample}_1.sub20.fq.gz \
        --fastq2 $OUTPUT_DIR/${sample}_2.sub20.fq.gz \
        --fastq_quality_shift \
        --fastq_read_length 125 \
        --duplication \
        --insert_size 350 \
        --run_mode 2 \
        --filter_taxa_method 1 \
        --requiring_taxa 'Mammalia' \
        --species_name 'Dama gazella' &> $MITOZ_OUT/${sample}.mitoz.log

    echo "MitoZ assembly for $sample completed."
done

#---------------Extract Protein-Coding Genes (PCGs)---------------#
# MitoZ output directory structure: ${MITOZ_OUT}/${sample}/Genome_Annotation/MT_genes.fasta
echo "### Extracting PCGs from MitoZ output ###"
for sample in "${samples[@]}"; do
    cp $MITOZ_OUT/${sample}/Genome_Annotation/MT_genes.fasta $PCG_OUT/${sample}_PCGs.fasta
done

#---------------Combine All PCGs for Alignment---------------#
echo "### Combining all PCGs ###"
cat $PCG_OUT/*_PCGs.fasta > $ALIGN_OUT/Dama_gazelle_PCGs.fasta

#---------------Align Sequences with MAFFT---------------#
echo "### Running MAFFT alignment ###"
mafft --thread 10 --adjustdirection $ALIGN_OUT/Dama_gazelle_PCGs.fasta > $ALIGN_OUT/Dama_gazelle_PCGs.aln

#---------------Build Phylogenetic Tree with IQ-TREE---------------#
echo "### Building IQ-TREE phylogenetic tree ###"
iqtree -nt 10 \
       --pre $TREE_OUT/Dama_gazelle_PCGs \
       -s $ALIGN_OUT/Dama_gazelle_PCGs.aln \
       -B 1000 \
       -m MFP  # ModelFinder to select best substitution model automatically

echo "### Workflow Complete ###"

