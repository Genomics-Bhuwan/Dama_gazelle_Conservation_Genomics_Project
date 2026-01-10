#### There are many method to make the SNAILPLOT.
- I am using https://github.com/hanwnetao/snailplot-assembly-stats

#### Change the directory to assemlby-stats, use the perl secript to convert fasta into JSON for both assemlby as given below.
#### Installation of the repository.
```bash
cd /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Assemlby_Statistics/SNAILPLOT

#### Install the blobtoolkit using apptainer or singularity.
apptainer pull docker://genomehubs/blobtoolkit:latest
```

#### These will run blobtool snail analysis for producing the SNAILPLOT.
```bash
# Process Hifiasm
apptainer exec -B /shared:/shared blobtoolkit_latest.sif blobtools create \
     --fasta Dama_gazelle_hifiasm-ULONT_primary.fasta \
     Dama_Hifiasm_Dir

# Process SCBI Reference
apptainer exec -B /shared:/shared blobtoolkit_latest.sif blobtools create \
     --fasta GCA_019969365.1_SCBI_Ndam_1.0_genomic.fna \
     SCBI_Ndam_Dir

# Process ORGONE Reference
apptainer exec -B /shared:/shared blobtoolkit_latest.sif blobtools create \
     --fasta GCA_917880005.1_ORGONE_02_genomic.fna \
     ORGONE_Dir
```


#### Visualize the plot for all the three reference genome assmlby.
- Generate PNG for Hifiasm
  ```bash
# 1. Create an images folder to keep things clean
mkdir -p images

# 2. Run the command pointing --out to that folder
apptainer exec blobtoolkit_latest.sif blobtools add \
  --busco /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Assemlby_Statistics/SNAILPLOT/T2TBUSCO_output/run_mammalia_odb10/full_table.tsv \
  /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Assemlby_Statistics/SNAILPLOT/Dama_Hifiasm_Dir

```
- Generate PNG for SCBI
```
# Generate PNG for SCBI Reference
```bash
apptainer exec /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Assemlby_Statistics/SNAILPLOT/blobtoolkit_latest.sif \
  blobtools view \
    --plot \
    --view snail \
    --driver chromium \
    --out /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Assemlby_Statistics/SNAILPLOT/images \
    /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Assemlby_Statistics/SNAILPLOT/SCBI_Ndam_Dir
```
# Generate PNG for ORGONE Reference

```bash
Add BUSCO data to ORGONE
Run this to link the BUSCO results to the ORGONE dataset:

Bash

apptainer exec /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Assemlby_Statistics/SNAILPLOT/blobtoolkit_latest.sif \
  blobtools add \
    --busco /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Assemlby_Statistics/SNAILPLOT/GCA_917880005.1/full_table.tsv \
    /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Assemlby_Statistics/SNAILPLOT/ORGONE_Dir
2. Generate the Snail Plot PNG
Run this to render the final image using the Chromium driver:

Bash

apptainer exec /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Assemlby_Statistics/SNAILPLOT/blobtoolkit_latest.sif \
  blobtools view \
    --plot \
    --view snail \
    --driver chromium \
    --out /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Assemlby_Statistics/SNAILPLOT/images \
    /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Assemlby_Statistics/SNAILPLOT/ORGONE_Dir
```
