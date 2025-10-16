## Creating a sequence dictionary for a reference genome using Picard

----
```bash
#!/bin/bash -l
#SBATCH --job-name=Seq_Dict
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=logs/Seq_Dict_%A.out
#SBATCH --error=logs/Seq_Dict_%A.err

# Load Java
module load java

# Paths
PICARD_PATH=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/picard.jar
REFERENCE=/localscratch/bistbs/4_aligning_with_BWA_Mem_Final_1/5_Sorted_BAMs/6_ReadGroups/7_MergeSam/8_MarkDuplicates/Dama_gazelle_hifiasm-ULONT_primary.fasta.gz

# Create output folder for sequence dictionary
SEQ_DICT_DIR=9_Sequence_Dictionary
mkdir -p $SEQ_DICT_DIR

# Output dictionary path
DICT_OUTPUT=$SEQ_DICT_DIR/$(basename ${REFERENCE%.fasta.gz}.dict)

# Run Picard to create sequence dictionary
java -Xmx8G -jar $PICARD_PATH CreateSequenceDictionary \
      R=$REFERENCE \
      O=$DICT_OUTPUT

echo "✅ Sequence dictionary created at: $DICT_OUTPUT"

   
```
---
