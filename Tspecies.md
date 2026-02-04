####  I am using Tspeceis R code to get the divergence between the two sub-species.
- Link: [https://github.com/arzwa/wgd/](https://github.com/OrthoFinder/OrthoFinder)
- TSPECIES:https://github.com/limj0987/Tspecies
  - First of all, I am using OrthoFinder to get the Single copy orthologs using the link given above.
  - 

  #### Step 1. Install the Orthofinder using Conda 
```bash
conda create -n of3_env python=3.12
conda activate of3_env
conda install orthofinder
```
#### Step 2.  Next, you can run the following commands to install OrthoFinder inside the of3_env virtural environment.
```bash
cd OrthoFinder
python3 -m venv of3_env # Create an virtural environment named of3_env
. of3_env/bin/activate # Activate of3_env
pip install .
```
#### Step 3. Conda install the OrthoFinder
```bash
conda install -c bioconda orthofinder
```

#### Step 4. Run the Orthofiner

```bash
- Extract the longest transcripts (as recommended by the manual). Cause I have multiple transcripts and when we have that, we stick to the longest.
# Note: Since you installed via conda, the script is in your environment path
python /home/bistbs/miniforge3/envs/of3_env/bin/scripts_of/primary_transcript.py Addra.pep.clean.faa -o of_input/
python /home/bistbs/miniforge3/envs/of3_env/bin/scripts_of/primary_transcript.py Mohrr.pep.clean.faa -o of_input/
```

#### Step 5. Run the orthofinder
```bash
orthofinder -f of_input/ -t 24
```

#### Step 6. Input files for Calculateion of Ks once Orthofinder finishes.
- Perform every 1-to-1 gene pair: protein alignment
- Align the two protein sequences using tools like MAFFT.back translation.
- Use pal2nal to convert the protein alignment into a DNA codon alignment using my .cds.clean.fna files.
- Now, Run the KsKs_calculator on the codon alignment.
#### Step 6. A. Alignment using MAFFT
```bash
# 1. Create the specific output directory
# Define the output directory
OUT_DIR="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Ks_calculation/new_try/MAFFT_New_Final_try"

# Create the directory if it doesn't exist
mkdir -p $OUT_DIR

# Loop through all orthogroup fasta files and align them
for f in *.fa; do
    echo "Aligning $f..."
    # --auto: automatically selects strategy (FFT-NS-1, FFT-NS-2, etc.)
    # --quiet: suppresses progress output for a cleaner terminal
    mafft --auto --quiet "$f" > "${OUT_DIR}/${f%.fa}.aln"
done

echo "Done! All alignments are in $OUT_DIR"
```
#### Step 6.B. Match the protein IDs with the DNA pairs
```bash
# 1. Set the specific Directory for DNA pairs
DNA_OUT="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Ks_calculation/of_input/OrthoFinder/DNA_pairs"
SOURCE_DIR="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Ks_calculation"
ALN_DIR="$SOURCE_DIR/of_input/OrthoFinder/mafftalignment"

mkdir -p $DNA_OUT

# 2. Loop through your Protein Alignment files
cd $ALN_DIR

for aln in *.aln; do
    sample=$(basename $aln .aln)
    
    # Get the IDs from the protein alignment (removing the '>' for clean searching)
    ID1=$(grep ">" "$aln" | sed -n '1p' | sed 's/>//')
    ID2=$(grep ">" "$aln" | sed -n '2p' | sed 's/>//')

    # Pull matching DNA from the clean source files and save to your new folder
    # We use -A 1 to get the sequence line immediately following the header
    grep -A 1 -h "$ID1" "$SOURCE_DIR/Addra.cds.clean.fna" "$SOURCE_DIR/Mohrr.cds.clean.fna" | grep -v "\--" > "$DNA_OUT/${sample}.cds.fa"
    grep -A 1 -h "$ID2" "$SOURCE_DIR/Addra.cds.clean.fna" "$SOURCE_DIR/Mohrr.cds.clean.fna" | grep -v "\--" >> "$DNA_OUT/${sample}.cds.fa"
    
    echo "Successfully paired DNA for $sample in DNA_pairs/"
done
```
#### Step 6C. Converson of the protein alignment into a DNA codon alignment using the clean .fna files using pal2nal.

```bash
# 1. Define Directories
BASE_DIR="/shared/jezkovt_bistbs_shared/Dama_Gazelle_Project/Ks_calculation/of_input/OrthoFinder"
ALN_DIR="$BASE_DIR/mafftalignment"
DNA_DIR="$BASE_DIR/DNA_pairs"
PAL_OUT="$BASE_DIR/pal2nal_results"

mkdir -p $PAL_OUT

# 2. Loop for PAL2NAL
for aln in $ALN_DIR/*.aln; do
    sample=$(basename $aln .aln)
    
    # Path to pal2nal.pl - make sure this path is correct for your server
    # -nogap removes columns with any gaps
    # -nomismatch ensures the DNA matches the Protein code
    perl pal2nal.pl \
        "$aln" \
        "$DNA_DIR/${sample}.cds.fa" \
        -nogap -nomismatch \
        -output clustal > "$PAL_OUT/cds.${sample}.aln.clustal.fa"
    
    echo "PAL2NAL finished for $sample"
done
```
#### Step 6 D. 
```bash
# 3. Define AXT Output Directory
AXT_OUT="$BASE_DIR/axt_results"
mkdir -p $AXT_OUT

# 4. Loop for AXTConvertor
for f in $PAL_OUT/*.clustal.fa; do
    sample=$(basename $f .aln.clustal.fa)
    
    # AXTConvertor usage: [Input Clustal] [Output AXT] [Output AUX]
    # Note: We name it .axt for clarity
    AXTConvertor "$f" "$AXT_OUT/${sample}.axt" "$AXT_OUT/${sample}.aux"
    
    echo "Converted $sample to AXT format"
done
```
