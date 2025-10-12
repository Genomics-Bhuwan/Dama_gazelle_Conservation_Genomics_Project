---
title: "Genome_Annotation_Dama_gazelle Using EGAPx pipeline"
Author: "Bhuwan Singh Bist"

---

#### Pre-requisite for using EGAPx



```{batch}
cd /scratch/bistbs_new/Genome_Annotation
```

2.  Clone EGAPx and set up Python environment

```{bash}
git clone https://github.com/ncbi/egapx.git
cd egapx

python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

```

3.  Make your input YAML (input_Dama.yaml)

```{bash}
genome: /scratch/bistbs_new/Genome_Annotation/Dama_gazelle_hifiasm-ULONT_primary.fasta.gz
taxid:  67940
short_reads:
  - /scratch/bistbs_new/Genome_Annotation/SRR5647654.fastq.gz
locus_tag_prefix: damaga
```

4.  Run the local execution or the batch script

```{batch}
python3 ui/egapx.py input_Dama.yaml \
  -o dama_output \
  -w /scratch/bistbs_new/Genome_Annotation/workdir
```

5. Most important step. Fixing the slurm.org file. I used the aptainer container, not the Singularity.
#### I am using aptainer which is pretty easy. Based on what you have such as singularity or docker you may want to use. But aptainer worked for me.

```{batch}
a.	Download Python 3.9+ version
b.	Load the singularity or docker whichever you have. I am using aptainer as recommended by my admin.
```

#### Work flow

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

6. Re-run. Because the first run above created the file named "EGAPX_config".
```{batch}
python3 ui/egapx.py input_Dama.yaml \
  -o /scratch/bistbs_new/Genome_Annotation/dama_output \
  -w /scratch/bistbs_new/Genome_Annotation/egapx/workdir \
  -e slurm

```
6.  Output

```{batch}
Explanation:

-o dama_output → output folder

-w .../workdir → working directory for temporary files

#### After it finishes, you should find:

dama_output/complete.genomic.gff → annotated genome

dama_output/proteins/ → predicted proteins

dama_output/transcripts/ → predicted transcripts

dama_output/logs/ → pipeline logs
```
