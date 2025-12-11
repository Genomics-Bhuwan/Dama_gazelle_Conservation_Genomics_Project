##### SVbyEye Visualization: Addra vs Mohrr
##### Step 1. Install the packages..
```bash
required_packages <- c(
  "ggnewscale", "gggenes", "wesanderson", "randomcoloR",
  "ggplot2", "dplyr", "tibble", "magrittr", "scales",
  "stringr", "data.table", "ggforce", "devtools", "Biostrings"
)
```
- Install any missing packages
```bash
installed_packages <- rownames(installed.packages())
for(pkg in required_packages) {
  if(!pkg %in% installed_packages) {
    install.packages(pkg, repos="https://cloud.r-project.org")
  }
}
```
- Load all required packages
  ```bash
lapply(required_packages, library, character.only = TRUE)
```
- Install SVbyEye from GitHub if not already installed
```bash
if(!"SVbyEye" %in% installed_packages) {
  devtools::install_github("daewoooo/SVbyEye", branch = "master")
}
library(SVbyEye)
```
##### 2️. Load genome sequences (FASTA)
library(Biostrings)
- Use the full path
```bash
addra_fasta <- "/scratch/bistbs/Synteny_Analysis/SVbyEye/Addra_complete.genomic.fna"
mhorr_fasta <- "/scratch/bistbs/Synteny_Analysis/SVbyEye/Mohrr_complete.genomic.fna"
```
##### Load genome sequences (may take time and memory)
```bash
addra_genome <- readDNAStringSet(addra_fasta)
mhorr_genome <- readDNAStringSet(mhorr_fasta)
```
##### Check loaded sequences
```bash
addra_genome
mhorr_genome
```

##### Step 3. Load PAF alignment
##### A. Generate with minimap2 in terminal if not done:
```bash
./minimap2-2.30_x64-linux/minimap2 -x asm5 -c Addra_complete.genomic.fna Mohrr_complete.genomic.fna > Addra_vs_Mohrr.paf
```
##### B. Index the reference for faster repeated runs.
```bash
./minimap2-2.30_x64-linux/minimap2 -d Addra.mmi Addra_complete.genomic.fna
./minimap2-2.30_x64-linux/minimap2 -x asm5 -c Addra.mmi Mohrr_complete.genomic.fna > Addra_vs_Mohrr.paf
```

##### Step 4. Read the .paf file
paf <- readPaf("Addra_vs_Mohrr.paf")

##### Step 5. Filter and Flip the alignment
paf_filtered<-SVbyEye:::filterPaf(paf, min.align.len = 1000) # Keep only alignments >=1 kb
paf_filtered

#### Step 5. Load teh annotation.
```bash
library(rtracklayer)
```
- Import Addra and Mohrr annotation
```bash
addra_gff <- import.gff("/scratch/bistbs/Synteny_Analysis/SVbyEye/Addra_complete.genomic.gff")
mhorr_gff <- import.gff("/scratch/bistbs/Synteny_Analysis/SVbyEye/Mohrr_complete.genomic.gff")
```
##### Step 6. Pairwise Genome Alignment for Miropeats.
```bash
pairwise_plot <- SVbyEye:::plotMiro(
  paf.table = paf_filtered,
  min.deletion.size = 1000,   # highlight deletions ≥1kb
  min.insertion.size = 1000,  # highlight insertions ≥1kb
  highlight.sv = TRUE,        # turn on SV highlighting
  color.by = "strand"         # color by alignment direction
)

ggsave("Addra_vs_Mohrr_pairwise.jpeg", pairwise_plot, width = 12, height = 6, dpi = 300)

ggsave("Addra_vs_Mohrr_pairwise.pdf", pairwise_plot, width = 12, height = 6, dpi = 300)
```

##### Self-alignment of Addra (horizontal dotplot)
###Preparing the self alignment file for Addra and Mohrr
```bash
/scratch/bistbs/Synteny_Analysis/SVbyEye/minimap2-2.30_x64-linux

./minimap2 -x asm20 -c -t 24 \
  ../Addra_complete.genomic.fna \
  ../Addra_complete.genomic.fna \
  > ../Addra_self.paf

```
```bash
./minimap2 -x asm20 -c -t 24 \
  ../Mohrr_complete.genomic.fna \
  Mohrr_complete.genomic.fna \
  > Mohrr_self.paf
```

- Load the self PAF inside R
```bash
addra_self_paf <- SVbyEye:::readPaf("Addra_self.paf")
mhorr_self_paf <- SVbyEye:::readPaf("Mohrr_self.paf")
```
- Filter the self paf
```bash
addra_self_filtered <- SVbyEye:::filterPaf(addra_self_paf, min.align.len = 1000)
mhorr_self_filtered <- SVbyEye:::filterPaf(mhorr_self_paf, min.align.len = 1000)

```
- Visualization
```bash
library(dplyr)    # for pipes and data manipulation
library(magrittr) # also provides %>% (optional if dplyr is loaded)
```
##### I am not removing or filtering cause a lot less.
###### Addra self-dotplot
```bash
self_plot <- plotSelf(
  paf.table = addra_self_paf,  # raw self-PAF
  color.by = "identity",
  shape = "segment",
  sort.by = "position"
)

ggsave("Addra_self_dotplot.jpeg", self_plot, width = 12, height = 6, dpi = 300)

ggsave("Addra_self_dotplot.pdf", self_plot, width = 12, height = 6)
```

##### Mohrr self-dotplot

mhorr_self_plot <- plotSelf(
  paf.table = mhorr_self_paf,  # raw self-PAF
  color.by = "identity",
  shape = "segment",
  sort.by = "position"
)

# Save high-res JPEG
ggsave("Mohrr_self_dotplot.jpeg", mhorr_self_plot, width = 12, height = 6, dpi = 300)

# Save PDF
ggsave("Mohrr_self_dotplot.pdf", mhorr_self_plot, width = 12, height = 6)

```

###### Step 8. Plot the stacked or all-versus-all OR AVA plot.
```bash
# I am using paf without filtering.
# Use the raw PAF tibble directly
ava_plot <- SVbyEye:::plotAVA(
  paf.table = paf_filtered,   # use your PAF tibble directly
  binsize = 5000,             # 5 kb bins
  color.by = "identity",      # color by % identity
  highlight.sv = TRUE         # optionally highlight SVs
)

ggsave("Addra_vs_Mohrr_AVA_stacked.jpeg", ava_plot, width = 12, height = 8)

```
##### Step 9. Plot the Whole-genome overview plot
```bash
genome_overview_plot <- SVbyEye:::plotGenome(
  paf.table = paf_filtered, 
  binsize = 5000,
  color.by = "identity",
  highlight.sv = TRUE,
  min.deletion.size = 1000,   # 1 kb
  min.insertion.size = 1000   # 1 kb
)

ggsave("Addra_vs_Mohrr_genome_overview.jpeg", genome_overview_plot, width = 12, height = 6)
```

##### Step 10. Extract the structural variants-SVs and plot the size distribution.
```bash
sv_breaks <-SVbyEye:::breakPaf(paf_filtered, minSize = 1000) # SVs >=1 kb
write.csv(sv_breaks, "Addra_vs_Mohrr_SVs.csv", row.names = FALSE)
```
##### Create SV size column
```bash
sv_breaks$SV_size <- abs(sv_breaks$end_qry - sv_breaks$start_qry)
```
##### Plot SV size distribution
```bash
sv_size_plot <- ggplot(sv_breaks, aes(x = SV_size)) +
  geom_histogram(binwidth = 1000, fill = "steelblue", color = "black") +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title = "SV Size Distribution: Addra vs Mohrr",
    x = "SV size (bp)", y = "Count"
  ) +
  theme_minimal()
ggsave("Addra_vs_Mohrr_SV_size_distribution.png", sv_size_plot, width = 10, height = 6)
```
# ==============================================

### All done combined targeting only 1-17 autosomes only.
################################################################################
# SVbyEye Visualization Pipeline
# Addra vs Mohrr genome comparison (AUTOSOMES 1–17 ONLY)
#
# This script:
# ✔ Loads genomes, annotations, PAF files
# ✔ Restricts analysis to chromosomes chr1–chr17
# ✔ Generates:
#     - Pairwise Miropeats plot
#     - Self-dotplots for Addra & Mohrr
#     - All-versus-all (AVA) plot
#     - Whole-genome overview
#     - Structural Variant (SV) extraction
################################################################################


##############################
### 1. Load Required Packages
##############################
```bash
required_packages <- c(
  "ggnewscale", "gggenes", "wesanderson", "randomcoloR",
  "ggplot2", "dplyr", "tibble", "magrittr", "scales",
  "stringr", "data.table", "ggforce", "devtools", "Biostrings",
  "rtracklayer"
)

installed_packages <- rownames(installed.packages())
for(pkg in required_packages) {
  if(!pkg %in% installed_packages) {
    install.packages(pkg, repos="https://cloud.r-project.org")
  }
}

lapply(required_packages, library, character.only = TRUE)

# Install SVbyEye if missing
if(!"SVbyEye" %in% installed_packages) {
  devtools::install_github("daewoooo/SVbyEye", branch = "master")
}
library(SVbyEye)


#########################################
### 2. Define autosomes to keep (1–17)
#########################################

keep_chr <- paste0("chr", 1:17)


#########################################
### 3. Load Genome FASTA Files
#########################################

addra_fasta <- "/scratch/bistbs/Synteny_Analysis/SVbyEye/Addra_complete.genomic.fna"
mhorr_fasta <- "/scratch/bistbs/Synteny_Analysis/SVbyEye/Mohrr_complete.genomic.fna"

addra_genome <- readDNAStringSet(addra_fasta)
mhorr_genome <- readDNAStringSet(mhorr_fasta)

# Keep only chr1–chr17
addra_genome <- addra_genome[names(addra_genome) %in% keep_chr]
mhorr_genome <- mhorr_genome[names(mhorr_genome) %in% keep_chr]


#########################################
### 4. Load PAF Alignment (Addra vs Mohrr)
#########################################

paf <- SVbyEye:::readPaf("Addra_vs_Mohrr.paf")

# Filter to retain alignments only on chr1–chr17
paf <- paf %>% filter(qry.name %in% keep_chr & ref.name %in% keep_chr)

# Filter alignment length ≥1kb
paf_filtered <- SVbyEye:::filterPaf(paf, min.align.len = 1000)


#########################################
### 5. Load & Filter Annotation (GFF)
#########################################

addra_gff <- import.gff("/scratch/bistbs/Synteny_Analysis/SVbyEye/Addra_complete.genomic.gff")
mhorr_gff <- import.gff("/scratch/bistbs/Synteny_Analysis/SVbyEye/Mohrr_complete.genomic.gff")

# Keep only chr1–chr17
addra_gff <- addra_gff[seqnames(addra_gff) %in% keep_chr]
mhorr_gff <- mhorr_gff[seqnames(mhorr_gff) %in% keep_chr]


#########################################
### 6. Pairwise Genome Alignment Plot
#########################################

pairwise_plot <- SVbyEye:::plotMiro(
  paf.table = paf_filtered,
  min.deletion.size = 1000,
  min.insertion.size = 1000,
  highlight.sv = TRUE,
  color.by = "strand"
)

ggsave("Addra_vs_Mohrr_pairwise.jpeg", pairwise_plot, width = 12, height = 6, dpi = 300)
ggsave("Addra_vs_Mohrr_pairwise.pdf",  pairwise_plot, width = 12, height = 6, dpi = 300)


#########################################
### 7. Self-Alignment of Addra & Mohrr
#########################################

# Load self-PAF
addra_self <- SVbyEye:::readPaf("Addra_self.paf")
mhorr_self <- SVbyEye:::readPaf("Mohrr_self.paf")

# Filter to chr1–17
addra_self <- addra_self %>% filter(qry.name %in% keep_chr & ref.name %in% keep_chr)
mhorr_self <- mhorr_self %>% filter(qry.name %in% keep_chr & ref.name %in% keep_chr)

# Addra self-dotplot
addra_self_plot <- plotSelf(
  paf.table = addra_self,
  color.by = "identity",
  shape = "segment",
  sort.by = "position"
)

ggsave("Addra_self_dotplot.jpeg", addra_self_plot, width = 12, height = 6, dpi = 300)
ggsave("Addra_self_dotplot.pdf",  addra_self_plot, width = 12, height = 6)

# Mohrr self-dotplot
mhorr_self_plot <- plotSelf(
  paf.table = mhorr_self,
  color.by = "identity",
  shape = "segment",
  sort.by = "position"
)

ggsave("Mohrr_self_dotplot.jpeg", mhorr_self_plot, width = 12, height = 6, dpi = 300)
ggsave("Mohrr_self_dotplot.pdf",  mhorr_self_plot, width = 12, height = 6)


#########################################
### 8. AVA (All-vs-All) Plot
#########################################

ava_plot <- SVbyEye:::plotAVA(
  paf.table = paf_filtered,
  binsize = 5000,
  color.by = "identity",
  highlight.sv = TRUE
)

ggsave("Addra_vs_Mohrr_AVA_stacked.jpeg", ava_plot, width = 12, height = 8, dpi = 300)


#########################################
### 9. Whole-genome Overview Plot
#########################################

genome_overview_plot <- SVbyEye:::plotGenome(
  paf.table = paf_filtered,
  binsize = 5000,
  color.by = "identity",
  highlight.sv = TRUE,
  min.deletion.size = 1000,
  min.insertion.size = 1000
)

ggsave("Addra_vs_Mohrr_genome_overview.jpeg", genome_overview_plot, width = 12, height = 6, dpi = 300)


#########################################
### 10. SV Extraction + Size Distribution
#########################################

sv_breaks <- SVbyEye:::breakPaf(paf_filtered, minSize = 1000)
write.csv(sv_breaks, "Addra_vs_Mohrr_SVs.csv", row.names = FALSE)

# Add SV size column
sv_breaks$SV_size <- abs(sv_breaks$end_qry - sv_breaks$start_qry)

sv_size_plot <- ggplot(sv_breaks, aes(x = SV_size)) +
  geom_histogram(binwidth = 1000, fill = "steelblue", color = "black") +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title = "SV Size Distribution: Addra vs Mohrr",
    x = "SV size (bp)",
    y = "Count"
  ) +
  theme_minimal()

ggsave("Addra_vs_Mohrr_SV_size_distribution.jpeg", sv_size_plot, width = 10, height = 6, dpi = 300)
```



