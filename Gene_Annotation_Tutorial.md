# Gene Annotation and Prediction for RNAseq Data with the Ruminant T2T Genome Assembly of Dama Gazelle

## Prerequisites

```bash
a. Download Python 3.9+ version
b. Load Singularity or Docker (this example uses Singularity v3.7.1)
```

### Load Singularity

```bash
module load singularity-3.7.1
```

---

## Step 1: Change to Working Directory

```bash
cd /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Genome_Annotation
```

---

## Step 2: Clone EGAPx and Set Up Python Environment

```bash
git clone https://github.com/ncbi/egapx.git        # Clone EGAPx
cd egapx                                          # Move into egapx folder
pip install -r requirements.txt                   # Install dependencies
```

---

## Step 3: Prepare Your Input YAML (`input_Dama.yaml`)

```yaml
genome: /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Genome_Annotation/Dama_gazelle_hifiasm-ULONT_primary.fasta.gz
taxid: 9805
short_reads:
  - /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Genome_Annotation/SRR5647654.fastq.gz
locus_tag_prefix: damaga
```

---

## Step 4: Run EGAPx (Local Execution or Batch Script)

```bash
python3 ui/egapx.py input_Dama.yaml \
  -o dama_output \
  -w /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Genome_Annotation/workdir
```

---

## Step 5: Re-Run with SLURM (After config file is created)

```bash
python3 ui/egapx.py input_Dama.yaml \
  -o dama_output \
  -w /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Genome_Annotation/workdir \
  -e slurm
```

---

## Step 6: Output Files

- **dama_output/complete.genomic.gff** → Annotated genome
- **dama_output/proteins/** → Predicted proteins
- **dama_output/transcripts/** → Predicted transcripts
- **dama_output/logs/** → Pipeline logs

---

## Explanation of Output Flags

- `-o dama_output` : Output folder
- `-w .../workdir` : Working directory for temporary files

---

### Author: Bhuwan Singh Bist, Jezkova lab
