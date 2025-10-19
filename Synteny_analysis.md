# Synteny analysis for Dama gazelle(Addra gazelle & Mohrr gazelle)

**Author:** Bhuwan Singh Bist

**Date:** 2025-10-11


> **Note:** All code blocks are for **demonstration purposes only** and are not executed in this document.

---

# Dama Gazelle(Addra vs Mohrr gazelle sub-species) WGS synteny analysis

This tutorial demonstrates the workflow for whole-genome sequencing (WGS) synteny analysis in the Dama gazelle. Each step includes the commands in a copyable code block.

---

#### 1. The code is using EGAPx for gene annotation.
#### I have Addra gazelle T2T reference genome assembly sequenced (Dama_gazelle_hifiasm-ULONT_primary.fasta.gz) using Oxford Nanopore Technology(ONT).I used RNAseq (NCBI SRA (SRR5647654) data for gene annotation using EGAPx pipeline(I have my pipeline for this in the same github page from where you are seeing it). 
#### I have Mohrr gazelle scaffol-level genome assemlby (GCA_917880005.1.fna) downloaded from NCBI sequenced using ONT. I am also used RNAseq (NCBI SRA (SRR5647654) data for gene annotation using EGAPx pipeline for this sub-species.

#### I generated the output files from these two sub-species from Genome Annotation pipepline to be used in the below analysis of "Synteny Analysis".
#### I am using the pipleine of MCScanX: Multiple Collinearity Scan toolkit X version. The most popular synteny analysis tool in the world: https://github.com/wyp1125/MCScanX

```bash
#!/bin/bash -l
# To be submitted by: sbatch gemoma_Dama_gazelle.slurm

#SBATCH --time=300:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=256G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --job-name=GeMoMa_Dama_gazelle
#SBATCH --output=GeMoMa_Dama_gazelle_%j.out
#SBATCH --error=GeMoMa_Dama_gazelle_%j.err

```bash
# Step 1. Downloading the MCScanX
# Go to your working directory
cd /scratch/bistbs/Synteny_Analysis
# Download it
# Unzip and compile
unzip MCScanX.zip
cd MCScanX
module load gcc-9.2.0 
make
```

# Step 2. Input files
a. .GFF file for both species
The gff file of each species has alot of heading and data. MCScanX needs only four columns: Scafold, GeneID, Start and End. Since, the .gff file consits of many headers. Remove them for both sub-species.
```bash
i. For Addra:
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene"{
  split($9,a,";"); 
  for(i in a) if(a[i] ~ /^ID=/){split(a[i],b,"="); id=b[2]}
  print $1, id, $4, $5
}' Addra_complete.genomic.gff > Addra_genes_minimal.gff

ii. For Mohrr
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene"{
  split($9,a,";"); 
  for(i in a) if(a[i] ~ /^ID=/){split(a[i],b,"="); id=b[2]}
  print $1, id, $4, $5
}' Mohrr_complete.genomic.gff > Mohrr_genes_minimal.gff
```
# Step 3. Input files
It needs two input files:
a. Addra_Mohrr.blast and Addra_Mohrr.gff or .bed

# Step 3. Create a combined BLAST database and run BLASTP
a. Combined both protein files and make a BLAST database.
cat Addra_complete.proteins.faa Mohrr_complete.proteins.faa > Addra_Mohrr_combined.faa

b. Make BLAST database
makeblastdb -in Addra_Mohrr_combined.faa -dbtype prot

c.Run BLASTP (Limit to top 5 hits per gene as MCScanX recommends);
## This produces the BLAST output file in the required m8 format.
blastp -query Addra_Mohrr_combined.faa \
       -db Addra_Mohrr_combined.faa \
       -evalue 1e-10 -num_threads 24 -outfmt 6 -max_target_seqs 5 \
       -out Addra_Mohrr.blast


# Step 4. Generate the .gff or .bed file for MCScanX.
# The software needs a tab-delimited file with gene locations as shown below.
chrID    start_position    end_position    gene_id
# Generate this from your GFF files using awk
You can generate this from your GFF files using awk:

# 4.a. For Addra gazelle
awk '$3=="gene" {split($9,a,";"); for(i in a){if(a[i]~/ID=/){split(a[i],b,"="); print "Ad" $1, $4, $5, b[2]}}}' Addra_complete.genomic.gff > Addra.gff

# 4.b. For Mohrr gazelle
awk '$3=="gene" {split($9,a,";"); for(i in a){if(a[i]~/ID=/){split(a[i],b,"="); print "Mh" $1, $4, $5, b[2]}}}' Mohrr_complete.genomic.gff > Mohrr.gff

# 4.c. Combined both files
cat Addra.gff Mohrr.gff > Addra_Mohrr.gff

#Make sure the chromosome IDs have unique prefixes like Ad1, Ad2, Mh1, Mh2, etc. (MCScanX uses these prefixes to distinguish species).

# 5. Run MCScanX
# This will use Addra_Mohrr.blast and Addra_Mohrr.gff and will produce teh Addra_Mohrr.collinearity(syntenic block file) as well as the Addra_Mohrr.html(visualization files).

./MCScanX Addra_Mohrr
---


#########################################################################################################################
#########################################################################################################################
## Gene Annotation Using "UniProt or SwissProt) for both sub-species.
#########################################################################################################################
Step 1. Module load
module load blast-2.13.0+
which blastp
which makeblastdb

## Step 2. Download Swiss Prot Database. Download in the working directory.
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

# Make BLAST database
makeblastdb -in uniprot_sprot.fasta -dbtype prot -out swissprot_db

## Step 3. Run Blast for both sub-species.
# Addra gazelle
blastp -query Addra_complete.proteins.faa \
       -db swissprot_db \
       -evalue 1e-5 \
       -out Addra_vs_swissprot.tsv \
       -outfmt 6 \
       -max_target_seqs 5 \
       -num_threads 24

# Mohrr gazelle
blastp -query Mohrr_complete.proteins.faa \
       -db swissprot_db \
       -evalue 1e-5 \
       -out Mohrr_vs_swissprot.tsv \
       -outfmt 6 \
       -max_target_seqs 5 \
       -num_threads 24


#########################################################################################################################
#########################################################################################################################
## Gene Annotation Using "TrEMBL) for both sub-species.
#########################################################################################################################

## Download the TrEMBL database
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
gunzip uniprot_trembl.fasta.gz

## Make blast database using TrEMBL
makeblastdb -in uniprot_trembl.fasta -dbtype prot -out trembl_db

## Run BLASTP for Addra gazelle
```bash
blastp -query Addra_complete.proteins.faa \
       -db trembl_db \
       -evalue 1e-5 \
       -out Addra_vs_trembl.tsv \
       -outfmt 6 \
       -max_target_seqs 5 \
       -num_threads 24
```
 ## Run BLASTP for Mohrr gazelle.
```bash
blastp -query Mohrr_complete.proteins.faa \
       -db trembl_db \
       -evalue 1e-5 \
       -out Mohrr_vs_trembl.tsv \
       -outfmt 6 \
       -max_target_seqs 5 \
       -num_threads 24
```
 


#########################################################################################################################
#########################################################################################################################
## Gene Annotation Using "InterProScan for both sub-species.
#########################################################################################################################
## It is a database integrating the predictive information about proteins function from a number of partner resources proving overview of the families that a protein belongs to and the domains and sites it contains. 
## If you have a protein fasta file and want to functionally characterize can use this software to run the scanning algorithms from the InterPro database.
## Since, our protein sequences data is in fasta from and the matches  are then calcualted against all of the required member databases signatures and the results are then output in variety of format. 


## Step 1. Download InterProScan
```bash
cd /scratch/bistbs/Synteny_Analysis/InterProScan_Annotation
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.50-84.0/interproscan-5.50-84.0-64-bit.tar.gz
## Extract the package.
tar -zxvf interproscan-5.50-84.0-64-bit.tar.gz
```

## Step 2. Change the directory
```bash
cd interproscan-5.50-84.0
module load hmmer-3.3.2
cd /scratch/bistbs/Synteny_Analysis/InterProScan_Annotation/interproscan-5.50-84.0/data/panther/15.0
hmmpress panther.hmm
```

## Step 3. Run the interproscan for Addra
```bash
#!/bin/bash -l
# Submit with: sbatch slurm_Addra.sh
#SBATCH --time=300:00:00          # Maximum runtime (5 days)
#SBATCH --nodes=1                 
#SBATCH --ntasks-per-node=24      
#SBATCH --mem=128G                
#SBATCH --partition=batch         
#SBATCH --mail-type=BEGIN,END     
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --job-name=Addra_InterProScan

# Load Java 
module load java-20

# Set InterProScan directory
IPR_DIR=/scratch/bistbs/Synteny_Analysis/InterProScan_Annotation/interproscan-5.50-84.0

# Input protein file
ADDRA_PROT=/scratch/bistbs/Synteny_Analysis/InterProScan_Annotation/Addra_complete.proteins.faa

# Output directory
ADDRA_OUT=/scratch/bistbs/Synteny_Analysis/InterProScan_Annotation/Addra_output
mkdir -p $ADDRA_OUT

# Run InterProScan
$IPR_DIR/interproscan.sh -i $ADDRA_PROT \
                         -f tsv \
                         -dp \
                         -goterms \
                         -pa \
                         -cpu 24 \
                         -o $ADDRA_OUT/Addra_interproscan.out

```

 ## Run the interproscan for Mohrr gazelle
```bash 
#!/bin/bash -l
# Submit with: sbatch slurm_Mohrr.sh
#SBATCH --time=300:00:00          # Maximum runtime (5 days)
#SBATCH --nodes=1                 
#SBATCH --ntasks-per-node=24      
#SBATCH --mem=128G                
#SBATCH --partition=batch         
#SBATCH --mail-type=BEGIN,END     
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --job-name=Mohrr_InterProScan

# Load Java if needed (InterProScan requires Java 11+)
# module load java-20

# Set InterProScan directory
IPR_DIR=/scratch/bistbs/Synteny_Analysis/InterProScan_Annotation/interproscan-5.50-84.0

# Input protein file
MOHRR_PROT=/scratch/bistbs/Synteny_Analysis/InterProScan_Annotation/Mohrr_complete.proteins.faa

# Output directory
MOHRR_OUT=/scratch/bistbs/Synteny_Analysis/InterProScan_Annotation/Mohrr_output
mkdir -p $MOHRR_OUT

# Run InterProScan
$IPR_DIR/interproscan.sh -i $MOHRR_PROT \
                         -f tsv \
                         -dp \
                         -goterms \
                         -pa \
                         -cpu 24 \
                         -o $MOHRR_OUT/Mohrr_interproscan.out

```






###################################################################################
###################################################################################
#### Genome Annotation Summary for both species using .gff file as input
###################################################################################
# 1. Addra gazelle
```bash
   grep -v "^#" Addra_complete.genomic.gff | awk '{print $3}' | sort | uniq -c
```

# 2. Mohrr gazelle
   ```bash
grep -v "^#" Mohrr_complete.genomic.gff | awk '{print $3}' | sort | uniq -c
```

