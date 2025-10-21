## InterProScan
   Change directory to the output file of InterProScan.
```bash
cd /scratch/bistbs/Synteny_Analysis/InterProScan_Annotation/InterProScan_output_for_both
```
----
### Step 1: Extract GO, KEGG, MetaCyc, and Reactome terms for each species
#### This will give you 4 summary files per species showing how many proteins are linked to each functional term.
```bash
🦌 For Addra
# GO terms
grep -o 'GO:[0-9]*' Addra_interproscan_combined.out | sort | uniq -c | sort -nr > Addra_GO_counts.txt

# KEGG pathways
grep -o 'KEGG:[^|]*' Addra_interproscan_combined.out | sort | uniq -c | sort -nr > Addra_KEGG_counts.txt

# MetaCyc pathways
grep -o 'MetaCyc:[^|]*' Addra_interproscan_combined.out | sort | uniq -c | sort -nr > Addra_MetaCyc_counts.txt

# Reactome pathways
grep -o 'Reactome:[^|]*' Addra_interproscan_combined.out | sort | uniq -c | sort -nr > Addra_Reactome_counts.txt
```

```bash
🐐 For Mohrr
# GO terms
grep -o 'GO:[0-9]*' Mohrr_interproscan.out | sort | uniq -c | sort -nr > Mohrr_GO_counts.txt

# KEGG pathways
grep -o 'KEGG:[^|]*' Mohrr_interproscan.out | sort | uniq -c | sort -nr > Mohrr_KEGG_counts.txt

# MetaCyc pathways
grep -o 'MetaCyc:[^|]*' Mohrr_interproscan.out | sort | uniq -c | sort -nr > Mohrr_MetaCyc_counts.txt

# Reactome pathways
grep -o 'Reactome:[^|]*' Mohrr_interproscan.out | sort | uniq -c | sort -nr > Mohrr_Reactome_counts.txt
```

### 🧠 Step 3: Compare the results between species
#### A. Compare the counts between Addra and Mohrr to see which reactome pathways are unique or shared.
```bash
cut -d' ' -f2 Addra_Reactome_counts.txt > Addra_Reactome_terms.txt
cut -d' ' -f2 Mohrr_Reactome_counts.txt > Mohrr_Reactome_terms.txt

##### Shared Reactome pathways
comm -12 <(sort Addra_Reactome_terms.txt) <(sort Mohrr_Reactome_terms.txt) > Shared_Reactome.txt

##### Unique to Addra
comm -23 <(sort Addra_Reactome_terms.txt) <(sort Mohrr_Reactome_terms.txt) > Addra_unique_Reactome.txt

##### Unique to Mohrr
comm -13 <(sort Addra_Reactome_terms.txt) <(sort Mohrr_Reactome_terms.txt) > Mohrr_unique_Reactome.txt
```

### B: GO pathways
cut -d' ' -f2 Addra_GO_counts.txt > Addra_GO_terms.txt
cut -d' ' -f2 Mohrr_GO_counts.txt > Mohrr_GO_terms.txt

#### Shared GO terms
comm -12 <(sort Addra_GO_terms.txt) <(sort Mohrr_GO_terms.txt) > Shared_GO.txt

#### Unique to Addra
comm -23 <(sort Addra_GO_terms.txt) <(sort Mohrr_GO_terms.txt) > Addra_unique_GO.txt

#### Unique to Mohrr
comm -13 <(sort Addra_GO_terms.txt) <(sort Mohrr_GO_terms.txt) > Mohrr_unique_GO.txt


#### KEGG pathways
cut -d' ' -f2 Addra_KEGG_counts.txt > Addra_KEGG_terms.txt
cut -d' ' -f2 Mohrr_KEGG_counts.txt > Mohrr_KEGG_terms.txt

#### Shared KEGG pathways
comm -12 <(sort Addra_KEGG_terms.txt) <(sort Mohrr_KEGG_terms.txt) > Shared_KEGG.txt

#### Unique to Addra
comm -23 <(sort Addra_KEGG_terms.txt) <(sort Mohrr_KEGG_terms.txt) > Addra_unique_KEGG.txt

#### Unique to Mohrr
comm -13 <(sort Addra_KEGG_terms.txt) <(sort Mohrr_KEGG_terms.txt) > Mohrr_unique_KEGG.txt


```
---
