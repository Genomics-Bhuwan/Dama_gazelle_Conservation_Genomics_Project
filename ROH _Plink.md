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

##### ---------------------------------------
##### ROH Analysis Script for Dama Gazelle
##### Working Directory: F:/Collaborative_Projects/Dama_Gazelle_Project/ROH/Plink_Final
##### ---------------------------------------

##### Load required packages
```bash
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
```
##### -----------------------------
##### Set Working Directory & Parameters
##### -----------------------------
```bash
setwd("F:/Collaborative_Projects/Dama_Gazelle_Project/ROH/Plink_Final")

genome_size <- 3e9  # 3 Gb genome
samples <- c("SRR17134085","SRR17134086","SRR17129394", # Addra
             "SRR17134087","SRR17134088")               # Mohrr
output_dir <- getwd()  # use current directory for outputs
```
##### -----------------------------
##### 1️⃣ Percent Genome in ROH per Subspecies
##### -----------------------------
```bash
roh_indiv <- fread(file.path(output_dir, "Dama_gazelle_ROH.hom.indiv"))

roh_sub <- roh_indiv %>%
  filter(IID %in% samples) %>%
  mutate(Species = ifelse(IID %in% c("SRR17134085","SRR17134086","SRR17129394"), 
                          "Addra", "Mhorr"),
         Percent_ROH = (KB * 1000 / genome_size) * 100)

p1 <- ggplot(roh_sub, aes(x=Species, y=Percent_ROH, fill=Species)) +
  geom_boxplot() +
  geom_jitter(width=0.1, size=3) +
  labs(#title="Percent Genome in ROH per Subspecies",
       x="Sub-species",
       y="% Genome in ROH",
       fill="Sub-species") +
  scale_fill_manual(values=c("Addra"="skyblue","Mhorr"="orange")) +
  scale_y_continuous(breaks = seq(0, max(roh_sub$Percent_ROH, na.rm=TRUE), by = 2)) +
  theme_classic() +
  theme(
    panel.border = element_rect(color="black", fill=NA, size=1),
    
    # Legend bottom-right
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill="white", colour="black", size=0.6),
    
    # Brighter Y-axis tick values
    axis.text.y = element_text(color="black", size=13, face="plain"),
    
    # Brighter X-axis tick values
    axis.text.x = element_text(color="black", size=13, face="plain"),
    
    # Brighter axis labels
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )

print(p1)
ggsave(file.path(output_dir, "Percent_ROH_by_subspecies.jpeg"), p1, width=6, height=4, dpi=300)
ggsave(file.path(output_dir, "Percent_ROH_by_subspecies.pdf"), p1, width=6, height=4)
```

##### -----------------------------
##### 2️⃣ ROH Landscape across Chromosomes
##### -----------------------------
```bash
roh_segments <- fread(file.path(output_dir, "Dama_gazelle_ROH.hom"))

roh_segments <- roh_segments %>%
  filter(IID %in% samples) %>%
  mutate(Species = ifelse(IID %in% c("SRR17134085","SRR17134086","SRR17129394"),
                          "Addra", "Mhorr"),
         ROH_Mb = KB / 1000)

roh_segments$CHR <- factor(roh_segments$CHR, levels=unique(roh_segments$CHR))

p2 <- ggplot(roh_segments) +
  geom_segment(aes(x=POS1, xend=POS2, y=IID, yend=IID, color=Species), linewidth=3) +
  facet_wrap(~CHR, scales="free_x") +
  theme_minimal() +
  labs(title="ROH Landscape across Individuals",
       x="Chromosome Position (bp)",
       y="Individual") +
  scale_color_manual(values=c("Addra"="skyblue","Mhorr"="orange"))

print(p2)
ggsave(file.path(output_dir, "ROH_landscape.jpeg"), p2, width=12, height=6)
ggsave(file.path(output_dir, "ROH_landscape.pdf"), p2, width=12, height=6)

```

#####-----------------------------
##### 3️⃣ Percent Genome in ROH by Size Category (3 groups)
##### -----------------------------
```bash
roh_segments_filtered <- roh_segments %>%
  filter(KB * 1000 > 100000) %>%  # keep ROH > 100 kb
  mutate(ROH_category = case_when(
    KB >= 100 & KB < 1000 ~ "0.1-1Mb",       # 0.1–1 Mb
    KB >= 1000 & KB < 5000 ~ "1-5Mb",        # 1–5 Mb
    KB >= 5000 & KB < 10000 ~ "5-10Mb"       # 5–10 Mb
  )) %>%
  filter(!is.na(ROH_category))  # remove any NA categories

roh_by_cat <- roh_segments_filtered %>%
  group_by(IID, Species, ROH_category) %>%
  summarise(Total_ROH_bp = sum(KB * 1000), .groups="drop") %>%
  mutate(Percent_Genome_ROH = (Total_ROH_bp / genome_size) * 100)

# Fixed colors to match new categories
p3 <- ggplot(roh_by_cat, aes(x=IID, y=Percent_Genome_ROH, fill=ROH_category)) +
  geom_bar(stat="identity", color="black") +
  facet_wrap(~Species, scales="free_x") +
  labs(
    x="Individual",
    y="% Genome in ROH",
    fill="ROH size category") +
  scale_fill_manual(values=c(
    "0.1-1Mb"="#e6ab02", 
    "1-5Mb"="#d95f02",
    "5-10Mb"="#1b9e77"
  )) +
  theme_classic() +
  theme(
    panel.border = element_rect(color="black", fill=NA, size=1),
    axis.text.y = element_text(color="black", size=12),
    axis.text.x = element_text(color="black", size=12, angle=0, hjust=0.5),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    
    # Legend top-left corner
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill="white", colour="black", size=0.6)
  )
p3 <- p3 + theme(
  strip.text = element_text(color="black", size=14, face="bold")
)

print(p3)

ggsave(file.path(output_dir, "Percent_ROH_by_size_category.jpeg"), p3, width=8, height=5, dpi=300)
ggsave(file.path(output_dir, "Percent_ROH_by_size_category.pdf"), p3, width=8, height=5)
```

##### -----------------------------
##### 4️⃣ Cumulative % Genome in ROH vs ROH Size
##### -----------------------------
```bash
roh_segments_filtered <- roh_segments_filtered[order(IID, ROH_Mb)]

roh_segments_filtered <- roh_segments_filtered %>%
  group_by(IID) %>%
  mutate(Cumulative_bp = cumsum(KB * 1000),
         Percent_Cumulative = (Cumulative_bp / genome_size) * 100) %>%
  ungroup()

p4 <- ggplot(roh_segments_filtered, aes(x=ROH_Mb, y=Percent_Cumulative, color=Species)) +
  geom_line(size=1.2) +
  geom_point(alpha=0.7, size=2) +
  labs(
    x="ROH Size (Mb)",
    y="Cumulative % of Genome in ROH",
    color="Sub-species"
  ) +
  scale_color_manual(values=c("Addra"="skyblue", "Mhorr"="orange")) +
  theme_classic() +
  theme(
    panel.border = element_rect(color="black", fill=NA, size=1),
    
    # Brighter axis titles
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    
    # Brighter axis tick labels
    axis.text.x = element_text(color="black", size=12),
    axis.text.y = element_text(color="black", size=12),
    
    # Legend slightly lower in bottom-right
    legend.position = c(0.95, 0.02),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill="white", colour="black", size=0.6)
  )

print(p4)

ggsave(file.path(output_dir, "Cumulative_ROH_vs_size.jpeg"), p4, width=8, height=5, dpi=300)
ggsave(file.path(output_dir, "Cumulative_ROH_vs_size.pdf"), p4, width=8, height=5)
```
##### -----------------------------
##### 5️⃣ NEW: Number of ROH vs Total ROH (Mb) per Size Category (3 categories)
##### -----------------------------
##### Assign Sub-species colors
```bash
subspecies_colors <- c("Addra"="skyblue", "Mhorr"="orange")

# Assign unique shapes to each sample
samples_shapes <- c("SRR17129394"=16,  # Circle
                    "SRR17134085"=17,  # Triangle
                    "SRR17134086"=15,  # Square
                    "SRR17134087"=18,  # Diamond
                    "SRR17134088"=8)   # Star

roh_threecat <- roh_segments %>%
  filter(KB * 1000 > 100000) %>%
  mutate(ROH_category = case_when(
    KB >= 100 & KB < 1000   ~ "0.1-1Mb",
    KB >= 1000 & KB < 5000  ~ "1-5Mb",
    KB >= 5000 & KB < 10000 ~ "5-10Mb",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(ROH_category))

roh_count_sum <- roh_threecat %>%
  group_by(IID, Species, ROH_category) %>%
  summarise(
    Num_ROH = n(),
    Sum_ROH_Mb = sum(ROH_Mb),
    .groups="drop"
  )

p5 <- ggplot(roh_count_sum, aes(x=Sum_ROH_Mb, y=Num_ROH)) +
  geom_point(aes(shape=IID, color=Species), size=4, alpha=0.9) +
  facet_wrap(~ROH_category) +
  labs(title="Number of ROH vs Total ROH Length (Mb) per Size Category",
       x="Total ROH Length (Mb)",
       y="Number of ROH",
       shape="Individual Sample",
       color="Sub-species") +
  scale_color_manual(values=subspecies_colors) +
  scale_shape_manual(values=samples_shapes) +
  theme_classic() +
  theme(
    panel.border = element_rect(color="black", fill=NA, size=1),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    axis.text.x = element_text(color="black", size=12),
    axis.text.y = element_text(color="black", size=12),
    
    # Legends in top-left corner
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill="white", colour="black", size=0.6)
  )

print(p5)

ggsave(file.path(output_dir, "ROH_Count_vs_Sum_shapes.jpeg"), p5, width=9, height=5, dpi=300)
ggsave(file.path(output_dir, "ROH_Count_vs_Sum_shapes.pdf"), p5, width=9, height=5)

```
