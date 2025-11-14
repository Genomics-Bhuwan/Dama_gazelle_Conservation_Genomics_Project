# Synteny analysis of Addra and Bos taurus
## 1. The code is using https://github.com/PombertLab/SYNY

```bash
## Download the synteny tool
wget https://github.com/Eco-Flow/synteny/archive/refs/tags/v4.0.0.tar.gz -O - | tar -xzvf -


# Download Nextflow 24.10.0 as version 25 fails bad.
wget -O nextflow https://github.com/nextflow-io/nextflow/releases/download/v24.10.0/nextflow

# Make it executable
chmod +x nextflow

# Verify version
./nextflow -v

# change directory
cd /scratch/bistbs_new/Synteny_analysis/synteny-4.0.0


# Run the pipeline offline.
./nextflow run main.nf -profile apptainer -offline \
  --input /scratch/bistbs_new/Synteny_analysis/2_Addra_vs_Bos_taurus/Addra_Bos_taurus.csv \
  --outdir /scratch/bistbs_new/Synteny_analysis/2_Addra_vs_Bos_taurus/Results \
  --clean true



```