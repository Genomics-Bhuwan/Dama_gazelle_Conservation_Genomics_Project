# Synteny analysis for Dama gazelle(Addra gazelle & Mohrr gazelle)

**Author:** Bhuwan Singh Bist

**Affiliation:** Jezkova lab

**Date:** 2025-10-11


> **Note:** All code blocks are for **demonstration purposes only** and are not executed in this document.

---

# Dama Gazelle WGS synteny analysis

This tutorial demonstrates the workflow for whole-genome sequencing (WGS) synteny analysis in the Dama gazelle. Each step includes the commands in a copyable code block.

---

## 1. The code is using GeMoMa for gene annotation.

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

# Masking and annotating repetitive elements with Repeatmodeler and RepeatMasker

# Repeatmodeler: 
a. used for identifying the repeates in the genome. It provides a list of repeat family sequence to mask repeats in the genome with RepeatMasker. 
b. screens DNA sequences for interspersed repeats and low complexity DNA sequences.
c. output would be: detailed annotation of the repeats that are present in the query sequence as well as modified version of the query sequence in which all the annotated repeats have been masked.
d. software takes way long time with large genomes. 
e. Set correct parameters in repeatemodeler so that we could get repeats that are not only grouped by family, but are also annotated.

# Link to the software
Repeatmodeler http://www.repeatmasker.org/RepeatModeler/
RepeatMasker http://www.repeatmasker.org/RMDownload.html

#usage:BuildDatabase -name {database_name} {genome_file-in_fasta_format}
# # Since, dama gazelle doesnot have the database in repbase, I am creating the denovo repeat library and use it.
## This step is building the database for downstream analysis.
```bash
BuildDatabase -name Dama_gazelle /scratch/bistbs_new/Genome_annotation_Homology_based/Dama_gazelle_hifiasm-ULONT_primary.fasta.gz
```


# RepeatModeler is a de novo transposable element (TE) family identification and modeling package. At the heart of RepeatModeler are three de-novo repeat finding programs ( RECON, RepeatScout and LtrHarvest/Ltr_retriever ) which employ complementary computational methods for identifying repeat element boundaries and family relationships from sequence data.

# RepeatModeler assists in automating the runs of the various algorithms given a genomic database, clustering redundant results, refining and classifying the families and producing a high quality library of TE families suitable for use with RepeatMasker and ultimately for submission to the Dfam database ( http://dfam.org ).
```bash
# Usage: RepeatModeler -database {database_name} -pa {number of cores} -LTRStruct > out.log
RepeatModeler -database Dama_gazelle -threads 24 -engine ncbi -LTRStruct  > repeatmodeler_Dama_gazelle_out.log
```


#Repeat Masker: It is used for masking the repetitive elements. Soft mask with lower case letter and hard mask with N.
# usage: RepeatMasker -pa 30 -gff -lib {consensi_classified} -dir {dir_name} {genome_in_fasta}

RepeatMasker -pa $NSLOTS -xsmall -gff -lib consensi.fa.classified -dir ../repeatmasker /path/to_assembly/bHypOws1_hifiasm.bp.p_ctg.fasta



