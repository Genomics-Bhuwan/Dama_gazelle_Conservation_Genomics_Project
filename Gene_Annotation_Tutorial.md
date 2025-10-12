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
cd /scratch/bistbs_new/Genome_Annotation
```

---

## Step 2: Clone EGAPx and Set Up Python Environment

```bash
git clone https://github.com/ncbi/egapx.git
cd egapx

python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

---

## Step 3: Prepare Your Input YAML (`input_Dama.yaml`)

```yaml
genome: /scratch/bistbs_new/Genome_Annotation/Dama_gazelle_hifiasm-ULONT_primary.fasta.gz
taxid:  67940
short_reads:
  - /scratch/bistbs_new/Genome_Annotation/SRR5647654.fastq.gz
locus_tag_prefix: damaga
```

---

## Step 4: Run EGAPx (Local Execution or Batch Script)

```bash
python3 ui/egapx.py input_Dama.yaml \
  -o dama_output \
  -w /scratch/bistbs_new/Genome_Annotation/workdir
```
---
## Step 5. Preparing the workflow. Most important step. Fixing the slurm.org file. I used the aptainer container, not the Singularity.
## I am using aptainer which is pretty easy. Based on what you have such as singularity or docker you may want to use. But aptainer worked for me.
1.Change directory to working directory

```{batch}
// ===============================
// SLURM CONFIG for EGAPx Pipeline (Apptainer) - bigmem
// ===============================
apptainer {
    enabled = true
    autoMounts = true
    path = "/usr/bin/apptainer"
    runOptions = ''
    envWhitelist = "https_proxy,http_proxy,ftp_proxy,DISPLAY,SLURM_JOB_ID,APPTAINER_BINDPATH"
}

env {
    APPTAINER_CACHEDIR = "/scratch/bistbs_new/Genome_Annotation/egapx/apptainer_cache"
    APPTAINER_TMPDIR   = "/scratch/bistbs_new/Genome_Annotation/egapx/tmp"
    DEBUG_STACK_TRACE_LEVEL = "Warning"
    EXCEPTION_STACK_TRACE_LEVEL = "Warning"
    DIAG_POST_LEVEL = "Trace"
    DEBUG_CATCH_UNHANDLED_EXCEPTIONS = "0"
}

process {
    executor = 'slurm'
    queue = 'batch'
    queueSize = 50
    pollInterval = '2 min'
    queueStatInterval = '5 min'
    submitRateLimit = '10/1min'
    retry.maxAttempts = 3
    
    // Apptainer container image
    container = '/software/egapx/ncbi-egapx-0.4.1-alpha.img'
    
    // Specific resource settings for miniprot - using bigmem partition
    withName: 'egapx:target_proteins_plane:miniprot:run_miniprot' {
        memory = '128.GB'  // Plenty available on mualhpcp27
        cpus = 12          // Max CPUs on bigmem nodes
        time = '300.h'
        clusterOptions = '--partition=batch'
        errorStrategy = 'retry'
        maxRetries = 2
    }
}
```


---

## Step 6: Re-Run with SLURM (After config file is created)

```bash
python3 ui/egapx.py input_Dama.yaml \
  -o /scratch/bistbs_new/Genome_Annotation/dama_output \
  -w /scratch/bistbs_new/Genome_Annotation/egapx/workdir \
  -e slurm
```

---

## Step 7: Output Files

- **dama_output/complete.genomic.gff** → Annotated genome
- **dama_output/proteins/** → Predicted proteins
- **dama_output/transcripts/** → Predicted transcripts
- **dama_output/logs/** → Pipeline logs

---

## Explanation of Output Flags

- `-o dama_output` : Output folder
- `-w .../workdir` : Working directory for temporary files
