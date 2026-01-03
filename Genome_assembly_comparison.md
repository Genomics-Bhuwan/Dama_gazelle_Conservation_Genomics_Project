#### Our project is comparing three reference genome assembly.
#### Two methods: 1. BUSCO summary statistics and another using 
- NCBI GCA_019969365.1: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_019969365.1/ (scaffold level).
- GCA_917880005.1 https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_917880005.1/ (scaffold level).
- Telomere 2 Telomere (T2T) reference genome assembly from this project: Dama_gazelle_hifiasm-ULONT_primary.fasta with T2T reference genome assembly.

#### Location of directories
```bash
PROJECT_DIR="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Assemlby_Statistics"
ASSEMBLY_STATS_DIR="$PROJECT_DIR/assembly-stats"
OUTPUT_DIR="$PROJECT_DIR/assembly_stats_output"

GENOME1="$PROJECT_DIR/GCA_019969365.1_SCBI_Ndam_1.0_genomic.fna"
GENOME2="$PROJECT_DIR/GCA_917880005.1_ORGONE_02_genomic.fna"
GENOME3="$PROJECT_DIR/Dama_gazelle_hifiasm-ULONT_primary.fasta"

ASSEMBLY_STATS_EXEC="$ASSEMBLY_STATS_DIR/build/assembly-stats"
```

#### Load the directories
```bash
cd "$PROJECT_DIR" || { echo "Project directory not found! Exiting."; exit 1; }

module load gcc-12.2.0
module load cmake/3.23.0
```

#### build the assembly stats
```bash
mkdir -p "$ASSEMBLY_STATS_DIR/build"
cd "$ASSEMBLY_STATS_DIR/build" || { echo "Failed to enter build folder! Exiting."; exit 1; }
```

#### Configure build with CMake
```bash
cmake ..
```
#### Compile the program
```bash
make
```
#### Make the output directory
```bash
mkdir -p "$OUTPUT_DIR"
```
#### Run assembly-stats on all three genomes
```bash
"$ASSEMBLY_STATS_EXEC" \
"$GENOME1" \
"$GENOME2" \
"$GENOME3" \
> "$OUTPUT_DIR/all_genomes_stats.txt"
```



#### Summary statistics using BUSCO genes

- For GCA_917880005.1
```bash
#!/bin/bash -l
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=90G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --job-name=5_BUSCO_Analysis


# Load Java and activate BUSCO environment
module load java-20
module load hmmer-3.3.2
conda activate busco-env

# Make sure BLAST+ and Miniprot are in PATH
export PATH=/software/blast+/ncbi-blast-2.13.0+/bin:$PATH
export PATH=/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/BUSCO/miniprot:$PATH
export PATH=/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/BUSCO/bbmap:$PATH

# Set Java heap space for BBTools to 90 GB
export BBMAP_JAVA_OPTS="-Xmx80g"

# Run BUSCO using localscratch paths
busco \
  -i /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/BUSCO/GCA_917880005.1_ORGONE_02_genomic.fna \
  -l mammalia_odb10 \
  -m genome \
  -o /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/BUSCO/GCA_917880005.1_Output \
  -c 20 \
  --miniprot \
  -f
```

- For GCA_019969365.1

```bash
#!/bin/bash -l
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=90G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --job-name=65_BUSCO_Analysis


# Load Java and activate BUSCO environment
module load java-20
module load hmmer-3.3.2
conda activate busco-env

# Make sure BLAST+ and Miniprot are in PATH
export PATH=/software/blast+/ncbi-blast-2.13.0+/bin:$PATH
export PATH=/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/BUSCO/miniprot:$PATH
export PATH=/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/BUSCO/bbmap:$PATH

# Set Java heap space for BBTools to 90 GB
export BBMAP_JAVA_OPTS="-Xmx80g"

# Run BUSCO using localscratch paths
busco \
  -i /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/BUSCO/GCA_019969365.1_SCBI_Ndam_1.0_genomic.fna \
  -l mammalia_odb10 \
  -m genome \
  -o /shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/BUSCO/GCA_019969365.1_Output \
  -c 20 \
  --miniprot \
  -f

```

