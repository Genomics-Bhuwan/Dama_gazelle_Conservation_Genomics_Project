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

```
---
