#### Our project is comparing three reference genome assembly.
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

