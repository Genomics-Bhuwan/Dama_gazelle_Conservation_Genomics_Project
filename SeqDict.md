## Creating a sequence dictionary for a reference genome using Picard

----
```bash
#!/bin/bash -l
#SBATCH --job-name=Seq_Dict
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=50G
#SBATCH --partition=batch
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --output=logs/Coverage_%A.out
#SBATCH --error=logs/Coverage_%A.err

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Load dependencies and set paths   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load Dependencies and setup env variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Jobscript to create sequence dictionary


# Load java
module java

# Specifiy some paths

PICARD_PATH=/exports/cmvm/eddie/eb/groups/ogden_grp/software/picard
REFERENCE=/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2022/data/reference/GCF_014754425.2_SCBI_Odam_1.1_genomic.fna

# Run Picard

java -Xmx4G -jar $PICARD_PATH/picard.jar CreateSequenceDictionary \
      R=$REFERENCE \
      O=${REFERENCE%.fasta}.dict
```
---