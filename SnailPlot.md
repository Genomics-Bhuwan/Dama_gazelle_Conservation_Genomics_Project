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

diamond makedb -p 16 --in reference_proteomes.fasta.gz --taxonmap reference_proteomes.taxid_map --taxonnodes ../taxdump/nodes.dmp -d reference_proteomes.dmnd
cd -
```

#### Fetch any BUSCO lineages that you plan to use
```bash
mkdir -p busco
wget -q -O eukaryota_odb10.gz "https://busco-data.ezlab.org/v4/data/lineages/eukaryota_odb10.2020-09-10.tar.gz" \
        && tar xf eukaryota_odb10.gz -C busco
```
