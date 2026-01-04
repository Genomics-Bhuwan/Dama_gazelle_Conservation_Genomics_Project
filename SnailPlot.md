#### I am using Blobtools 2 to assess the reference genome assemly quality as well as detect the contaminant in the reference genome assembly. 
---
This is what is to be done ultimately.
    # 1. Create a new BlobDir from a FASTA file:
    blobtools create --fasta examples/assembly.fasta BlobDir

    # 2. Create a new BlobDir from a BlobDB:
    blobtools create --blobdb examples/blobDB.json BlobDir

    # 3. Add Coverage data from a BAM file:
    blobtools add --cov examples/assembly.reads.bam BlobDir

    # 4. Assign taxonomy from BLAST hits:
    blobtools add add --hits examples/blast.out --taxdump ../taxdump BlobDir

    # 5. Add BUSCO results:
    blobtools add --busco examples/busco.tsv BlobDir

    # 6. Host an interactive viewer:
    blobtools host BlobDir

    # 7. Filter a BlobDir:
    blobtools filter --param length--Min=5000 --output BlobDir_len_gt_5000 BlobDir

---



#### Step 1. Installation
- https://blobtoolkit.genomehubs.org/install/

#### Create conda environment and also install python 3.9
```bash
  conda create -y -n btk -c conda-forge python=3.9
conda activate btk
```
#### Install blobtoolkit:
```bash
pip install "blobtoolkit[full]"

 blobtools -h  ## used to see parameters
```
#### Database
- Local copy of the NCBI taxdump newere format
- Copies of NCBI nucleotide(nt)
- Uniprot database
- aLL of these can be fetched automatically when running the BlobToolkit pipeline.
- Alternatively use the commands below to fetch copies for standalone use.

#### 1. Fetch the NCBI Taxdump
```bash
mkdir -p taxdump;
cd taxdump;
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -;
cd ..;
```

#### 2. Fetch the nt database
```bash
mkdir -p nt
cd nt
wget -r -nd -np -A "nt.*.tar.gz" ftp://ftp.ncbi.nlm.nih.gov/blast/db/

for f in nt.*.tar.gz; do
    tar -xvf "$f"
    rm "$f"
done

```

#### 3. Fetch and format the UniProt database
- Formatting the UniProt database requires the NCBI taxdump to be downloaded and uncompressed in a sister directory, as described in step 2 above.
```bash
mkdir -p uniprot
wget -q -O uniprot/reference_proteomes.tar.gz \
 ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$(curl \
     -vs ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ 2>&1 | \
     awk '/tar.gz/ {print $9}')
cd uniprot
tar xf reference_proteomes.tar.gz

touch reference_proteomes.fasta.gz
find . -mindepth 2 | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional' | xargs cat >> reference_proteomes.fasta.gz

echo -e "accession\taccession.version\ttaxid\tgi" > reference_proteomes.taxid_map
zcat */*/*.idmapping.gz | grep "NCBI_TaxID" | awk '{print $1 "\t" $1 "\t" $3 "\t" 0}' >> reference_proteomes.taxid_map
```
#### Step 4. Install Diamond and run the diamaond database hit
```bash
wget http://github.com/bbuchfink/diamond/releases/download/v2.1.17/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz

./diamond makedb -p 16 \
  --in uniprot/reference_proteomes.fasta.gz \
  --taxonmap uniprot/reference_proteomes.taxid_map \
  --taxonnodes taxdump/nodes.dmp \
  -d uniprot/reference_proteomes.dmnd

```

#### Fetch any BUSCO lineages that you plan to use
```bash
mkdir -p busco
wget -q -O eukaryota_odb10.gz "https://busco-data.ezlab.org/v4/data/lineages/eukaryota_odb10.2020-09-10.tar.gz" \
        && tar xf eukaryota_odb10.gz -C busco
```

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
####################### Final run with Bloobtools2 ################################################################################
###################################################################################################################################

#### STEP 1. Create BlobDir (correct for unpublished T2T genome)
```bash
blobtools create \
  --fasta Dama_gazelle_hifiasm-ULONT_primary.fasta \
  --taxdump taxdump \
  --taxid 59553 \
  BlobDir
```

#### STEP 2. Add taxonomic hits (MOST IMPORTANT STEP)
- For UL-ONT assemblies, DIAMOND blastx is strongly recommended.
```bash
diamond blastx \
  -d uniprot/reference_proteomes.dmnd \
  -q Dama_gazelle_hifiasm-ULONT_primary.fasta \
  -o diamond.out \
  -f 6 qseqid staxids bitscore evalue \
  --max-target-seqs 1 \
  --evalue 1e-25 \
  -p 32
```

- Then add hits:
```bash
blobtools add \
  --hits diamond.out \
  --taxdump taxdump \
  BlobDir
```
#### STEP 3. Add taxonomic hits (MOST IMPORTANT STEP)
- For mammalian T2T assemblies, bacterial contaminants will be very obvious here.

- Add BUSCO (optional but recommended)
- Even for T2T, BUSCO is useful for sanity checking.
```bash
blobtools add \
  --busco BUSCO_output/*/full_table.tsv \
  BlobDir
```

#### STEP 4. Coverage: what to do for UL-ONT ?
#### Option A — Skip coverage (acceptable for T2T)
- If: You only have UL-ONT reads; Assembly is already manually curated
- You may skip this step.
- BlobTools plots will still work (GC vs length).

#### Option B — Add UL-ONT coverage (recommended if you have reads)
```bash
minimap2 -ax map-ont \
  -t 32 \
  Dama_gazelle_hifiasm-ULONT_primary.fasta \
  ultralong_ONT_reads.fastq.gz \
| samtools sort -@16 -O BAM -o assembly.ont.bam -

samtools index assembly.ont.bam
```

- Then:
```bash
blobtools add \
  --cov assembly.ont.bam \
  BlobDir
```
- Expect very high, uneven coverage — this is normal for UL-ONT.

#### STEP 5. Launch the viewer
```bash
blobtools host BlobDir
```

#### Open the URL in your browser (via SSH tunnel if on HPC).
- What contamination looks like in T2T UL-ONT assemblies?
- Typical findings:
- Bacterial contigs: very short, extreme GC, bacterial taxonomy
- Mitochondrial genome: small (~16–17 kb), very high coverage


#### STEP 6. Filtering (example)
- Remove non-Chordata contigs:
```bash
blobtools filter \
  --param taxon--phylum=Chordata \
  --output BlobDir_chordata_only \
  BlobDir
```

- Remove tiny contigs:
```bash
blobtools filter \
  --param length--Min=50000 \
  --output BlobDir_len_gt_50kb \
  BlobDir
```
- Recommended final output for a T2T genome

- One contig per chromosome

- One mitochondrial contig

- Zero bacterial/fungal contigs
