#### Calculate the ROH using Plink

###### Convert VCF to PLINK binary format (BED/BIM/FAM)
```bash
plink --vcf /scratch/bistbs/Population_Genomic_Analysis/ROH/Plink/Dama_gazelle_biallelic_snps_filtered.recode.vcf \
      --make-bed \
 --allow-extra-chr \
      --out /scratch/bistbs/Population_Genomic_Analysis/ROH/Plink/Dama_gazelle
```

###### Calculate the ROH using Plink
```bash
plink --bfile /scratch/bistbs/Population_Genomic_Analysis/ROH/Plink/Dama_gazelle \
      --homozyg \
      --homozyg-window-snp 50 \
      --homozyg-snp 50 \
      --homozyg-kb 500 \
      --homozyg-gap 1000 \
      --homozyg-density 50 \
      --homozyg-window-missing 5 \
      --homozyg-window-het 3 \
--allow-extra-chr \
      --out /scratch/bistbs/Population_Genomic_Analysis/ROH/Plink/Dama_gazelle_ROH
```


#### Visualization of ROH 

```bash
# Load required packages
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

# -----------------------------
# Parameters
# -----------------------------
cd F:/Collaborative_Projects/Dama_Gazelle_Project/ROH/Plink_Final
genome_size <- 3e9  # 3 Gb
samples <- c("SRR17134085","SRR17134086","SRR17129394", # Addra
             "SRR17134087","SRR17134088")               # Mohrr
output_dir <- "F:/Collaborative_Projects/Dama_Gazelle_Project/ROH/Plink_Final"

# -----------------------------
# 1️⃣ Percent Genome in ROH per Subspecies
# -----------------------------
roh_indiv <- fread(file.path(output_dir, "Dama_gazelle_ROH.hom.indiv"))

# Filter samples and assign species
roh_sub <- roh_indiv %>%
  filter(IID %in% samples) %>%
  mutate(Species = ifelse(IID %in% c("SRR17134085","SRR17134086","SRR17129394"), "Addra", "Mohrr"),
         Percent_ROH = (KB*1000 / genome_size)*100)  # Convert KB to bp then %

# Plot
p1 <- ggplot(roh_sub, aes(x=Species, y=Percent_ROH, fill=Species)) +
  geom_boxplot() +
  geom_jitter(width=0.1, size=3) +
  theme_classic() +
  labs(title="Percent Genome in ROH per Subspecies",
       x="Subspecies",
       y="Percent genome in ROH") +
  scale_fill_manual(values=c("Addra"="skyblue","Mohrr"="orange"))

print(p1)
ggsave(file.path(output_dir,"Percent_ROH_by_sub_species.jpeg"), plot=p1, width=6, height=4)

# -----------------------------
# 2️⃣ ROH Landscape across Chromosomes
# -----------------------------
roh_segments <- fread(file.path(output_dir, "Dama_gazelle_ROH.hom"))

# Filter samples and assign species
roh_segments <- roh_segments %>%
  filter(IID %in% samples) %>%
  mutate(Species = ifelse(IID %in% c("SRR17134085","SRR17134086","SRR17129394"), "Addra", "Mohrr"),
         ROH_Mb = KB / 1000)  # Convert KB to Mb

# Ensure chromosome order
roh_segments$CHR <- factor(roh_segments$CHR, levels=unique(roh_segments$CHR))

# Plot
p2 <- ggplot(roh_segments) +
  geom_segment(aes(x=POS1, xend=POS2, y=IID, yend=IID, color=Species), linewidth=3) +
  facet_wrap(~ CHR, scales="free_x") +
  theme_minimal() +
  labs(title="ROH Landscape across Individuals",
       x="Chromosome Position (bp)",
       y="Individual") +
  scale_color_manual(values=c("Addra"="skyblue","Mohrr"="orange"))

print(p2)
ggsave(file.path(output_dir,"ROH_landscape.jpeg"), plot=p2, width=12, height=6)

# -----------------------------
# 3️⃣ Percent Genome in ROH by ROH Size Category
# -----------------------------
# Filter ROH > 100 kb
roh_segments_filtered <- roh_segments %>% filter(KB*1000 > 100000)

# Define ROH categories
roh_segments_filtered <- roh_segments_filtered %>%
  mutate(ROH_category = case_when(
    KB >= 0.1 & KB < 1000 ~ "0.1-1Mb",
    KB >= 1000 & KB < 5000 ~ "1-5Mb",
    KB >= 5000 & KB < 10000 ~ "5-10Mb",
    KB >= 10000 ~ "10-100Mb"
  ))

# Aggregate per individual per category
roh_by_cat <- roh_segments_filtered %>%
  group_by(IID, Species, ROH_category) %>%
  summarise(Total_ROH_bp = sum(KB*1000), .groups="drop") %>%
  mutate(Percent_Genome_ROH = (Total_ROH_bp / genome_size)*100)

# Plot stacked bar
p3 <- ggplot(roh_by_cat, aes(x=IID, y=Percent_Genome_ROH, fill=ROH_category)) +
  geom_bar(stat="identity", color="black") +
  facet_wrap(~Species, scales="free_x") +
  theme_classic() +
  labs(title="Percent Genome in ROH by Size Category",
       x="Individual",
       y="Percent of Genome in ROH",
       fill="ROH Size Category") +
  scale_fill_manual(values=c("0.1-1Mb"="#e6ab02", "1-5Mb"="#d95f02",
                             "5-10Mb"="#7570b3", "10-100Mb"="#1b9e77"))

print(p3)
ggsave(file.path(output_dir,"Percent_ROH_by_size_category.jpeg"), plot=p3, width=8, height=5)

# -----------------------------
# 4️⃣ Cumulative Percent Genome in ROH vs ROH Size
# -----------------------------
roh_segments_filtered <- roh_segments_filtered[order(IID, ROH_Mb)]
roh_segments_filtered <- roh_segments_filtered %>%
  group_by(IID) %>%
  mutate(Cumulative_bp = cumsum(KB*1000),
         Percent_Cumulative = (Cumulative_bp / genome_size)*100) %>%
  ungroup()

# Plot cumulative
p4 <- ggplot(roh_segments_filtered, aes(x=ROH_Mb, y=Percent_Cumulative, color=Species)) +
  geom_line(size=1.2) +
  geom_point(alpha=0.7, size=2) +
  theme_classic() +
  labs(title="Cumulative Percent Genome in ROH vs ROH Size",
       x="ROH Size (Mb)",
       y="Cumulative Percent of Genome in ROH") +
  scale_color_manual(values=c("Addra"="skyblue","Mohrr"="orange"))

print(p4)
ggsave(file.path(output_dir,"Cumulative_ROH_vs_size.jpeg"), plot=p4, width=8, height=5)
```
