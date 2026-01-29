#### If you want to detect and annotate gene within ROH islands: FOLLOW THIS PAPER:https://link.springer.com/article/10.1186/s12864-025-11616-8

#### Calculate the ROH using Plink

###### Convert VCF to PLINK binary format (BED/BIM/FAM)
```bash
plink --vcf /scratch/bistbs/Population_Genomic_Analysis/ROH/Plink/Dama_Gazelle_biallelic_Autosomes_Only_withID.vcf.gz \
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
```bash
#######################################################
#######################################################
######################################################
# ----------------------------------------------------
#  ROH Analysis Script for Dama Gazelle
# ----------------------------------------------------
#######################################################
#######################################################
######################################################
# ----------------------------------------------------
#  ROH Analysis Script for Dama Gazelle
# ----------------------------------------------------

# -----------------------------
# Load required packages
# -----------------------------
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(IRanges) # for merging overlapping ROHs
library(scales) # Needed for better log-scale breaks

# -----------------------------
# Set Working Directory & Parameters
# -----------------------------
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)

# -----------------------------
# Set Working Directory & Parameters
# -----------------------------
setwd("F:/Collaborative_Projects/Dama_Gazelle_Project/ROH/ROH_Final_Second_time")

# UPDATED: Using the actual autosomal size (2.56 Gb) for accurate F_ROH calculation
genome_size <- 2557520234 
samples <- c("SRR17134085","SRR17134086","SRR17129394", # Addra
             "SRR17134087","SRR17134088") # Mhorr
output_dir <- getwd()

# Define shared aesthetics
subspecies_colors <- c("Addra"="#56B4E9", "Mhorr"="#E69F00") # Colorblind friendly hex codes

# -----------------------------
# Load and Process Data
# -----------------------------
roh_indiv <- fread(file.path(output_dir, "Dama_gazelle_ROH.hom.indiv"))

roh_sub <- roh_indiv %>%
  filter(IID %in% samples) %>%
  mutate(Species = ifelse(IID %in% c("SRR17134085","SRR17134086","SRR17129394"),
                          "Addra", "Mhorr"),
         Percent_ROH = (KB * 1000 / genome_size) * 100)

# -----------------------------
# Create Publication Plot
# -----------------------------
p1 <- ggplot(roh_sub, aes(x=Species, y=Percent_ROH, fill=Species)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, size = 0.4) + 
  geom_jitter(width=0.1, size=2, shape=21, color="black", stroke=0.3) + 
  labs(x="Subspecies",
       y="Total ROH (% of genome)",
       fill="Subspecies") +
  scale_fill_manual(values=subspecies_colors) +
  scale_y_continuous(breaks = pretty_breaks(n = 6)) + 
  theme_classic(base_size = 10) + # Molecular Ecology prefers ~8-10pt font for reduction
  theme(
    # Positioning legend top-left as requested
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill="white", color="black", size=0.3),
    legend.title = element_text(size=8, face="bold"),
    legend.text = element_text(size=8),
    # Styling axes for print
    panel.border = element_rect(color="black", fill=NA, size=0.5),
    axis.line = element_blank(), # border replaces lines
    axis.text = element_text(color="black", size=9),
    axis.title = element_text(color="black", size=10, face="bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

# -----------------------------
# Exporting as Vector Graphics
# -----------------------------
# Molecular Ecology Single Column Width: 80 mm (approx 3.15 inches)
# We export at a height that preserves the aspect ratio, e.g., 80 mm width x 80 mm height
# 1. Save as Vector PDF (Recommended for Molecular Ecology)
ggsave(filename = file.path(output_dir, "Fig1_ROH_Percentage.pdf"), 
       plot = p1, 
       width = 80, 
       height = 80, 
       units = "mm", 
       device = "pdf",
       useDingbats = FALSE)

# 2. Save as Vector EPS (Standard for high-end print publishing)
# Note: ggsave uses the postscript device for .eps
ggsave(filename = file.path(output_dir, "Fig1_ROH_Percentage.eps"), 
       plot = p1, 
       width = 80, 
       height = 80, 
       units = "mm", 
       device = "eps")

# 3. Save as High-Resolution JPEG (Min 300 DPI as per Wiley guidelines)
# We set quality to 100 to ensure no compression artifacts
ggsave(filename = file.path(output_dir, "Fig1_ROH_Percentage.jpeg"), 
       plot = p1, 
       width = 80, 
       height = 80, 
       units = "mm", 
       dpi = 300, 
       device = "jpeg",
       quality = 100)

print("Figure successfully exported in PDF, EPS, and JPEG formats.")


# --------------------------------------------------------------------------------------
# 2.  Percent Genome in ROH by Size Category (Barplot) - FILTER FIXED
# --------------------------------------------------------------------------------------
# Use the loaded roh_segments data
roh_segments_3 <- roh_segments %>%
  filter(IID %in% samples) %>%
  mutate(Species = ifelse(IID %in% c("SRR17134085","SRR17134086","SRR17129394"),
                          "Addra", "Mhorr"),
         ROH_Mb = KB / 1000)

roh_segments_3$CHR <- factor(roh_segments_3$CHR, levels=unique(roh_segments_3$CHR))

roh_segments_filtered <- roh_segments_3 %>%
  # --- CRITICAL FIX: Removed the line: filter(KB * 1000 > 100000) ---
  mutate(ROH_category = case_when(
    KB >= 100 & KB < 1000 ~ "0.1-1Mb",
    KB >= 1000 & KB < 5000 ~ "1-5Mb",
    KB >= 5000 & KB < 10000 ~ "5-10Mb"
  )) %>%
  filter(!is.na(ROH_category)) # Filter out ROHs < 100 KB or > 10 Mb


roh_by_cat <- roh_segments_filtered %>%
  group_by(IID, Species, ROH_category) %>%
  summarise(Total_ROH_Mb = sum(KB / 1000), .groups="drop") %>%
  mutate(Percent_Genome_ROH = (Total_ROH_Mb * 1e6 / genome_size) * 100)

# --------------------------------------------------------------------------------------
# Final Plotting: Horizontal Labels, Full Bounding Box, and Facet Label Boxes
# --------------------------------------------------------------------------------------
p3 <- ggplot(roh_by_cat, aes(x=IID, y=Percent_Genome_ROH, fill=ROH_category)) +
  geom_bar(stat="identity", color="black", linewidth = 0.3) +
  facet_wrap(~Species, scales="free_x") +
  labs(x="Individual",
       y="Total ROH (% of genome)",
       fill="ROH size category") +
  scale_fill_manual(values=c("0.1-1Mb"="#e6ab02", "1-5Mb"="#d95f02", "5-10Mb"="#1b9e77")) +
  # Expand y-axis to 40% and ensure bars sit on the x-axis line
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10), expand = c(0,0)) +
  theme_classic(base_size = 10) +
  theme(
    # --- FULL EXTERNAL BOUNDING BOX ---
    panel.border = element_rect(color="black", fill=NA, linewidth=0.8),
    axis.line = element_blank(), 
    
    # --- FACET LABEL BOXES (Addra and Mhorr boxes) ---
    strip.background = element_rect(color="black", fill="gray95", linewidth=0.8),
    strip.text = element_text(color="black", size=10, face="bold"),
    
    # LEGEND: Top-left corner inside the box
    legend.position = c(0.05, 0.95),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill="white", colour="black", linewidth=0.3),
    
    # X-AXIS: Horizontal and readable
    axis.text.x = element_text(color="black", size=7, angle=0, hjust=0.5),
    axis.text.y = element_text(color="black", size=9),
    axis.title = element_text(color="black", size=10, face="bold")
  )

# --- Save for Molecular Ecology Resources ---
ggsave(
  filename = file.path(output_dir, "Percent_ROH_by_size_category_Final_FullBoxes.pdf"), 
  plot = p3, 
  width = 169, 
  height = 100, 
  units = "mm", 
  device = "pdf",
  useDingbats = FALSE
)

# --- Save as High-Res JPEG (300 DPI for standard image viewing) ---
ggsave(
  filename = file.path(output_dir, "Percent_ROH_by_size_category_Final.jpeg"), 
  plot = p3, 
  width = 169, 
  height = 100, 
  units = "mm", 
  dpi = 300,
  device = "jpeg",
  quality = 100
)

print("Final plots saved in PDF and 300 DPI JPEG formats.")


# --------------------------------------------------------------------------------------
# 3. Number of ROH vs Total ROH (Mb) per Size Category (Scatterplot)
# --------------------------------------------------------------------------------------

# Process Data: Ensure Species and ROH_Mb are created before grouping
roh_threecat <- roh_segments %>%
  filter(IID %in% samples) %>%
  mutate(
    Species = ifelse(IID %in% c("SRR17134085","SRR17134086","SRR17129394"), "Addra", "Mhorr"),
    ROH_Mb = KB / 1000
  ) %>%
  mutate(ROH_category = case_when(
    KB >= 100 & KB < 1000 ~ "0.1-1Mb",
    KB >= 1000 & KB < 5000 ~ "1-5Mb",
    KB >= 5000 & KB < 10000 ~ "5-10Mb",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(ROH_category)) %>%
  mutate(ROH_category = factor(ROH_category, levels=c("0.1-1Mb", "1-5Mb", "5-10Mb")))

# Summarize
roh_count_sum <- roh_threecat %>%
  group_by(IID, Species, ROH_category) %>%
  summarise(Num_ROH = n(),
            Sum_ROH_Mb = sum(ROH_Mb),
            .groups="drop")

# Create Plot
p5 <- ggplot(roh_count_sum, aes(x=Sum_ROH_Mb, y=Num_ROH)) +
  # Using shape for individuals and color for subspecies
  geom_point(aes(shape=IID, color=Species), size=3.5, alpha=0.85, stroke = 0.8) +
  facet_wrap(~ROH_category, scales = "free") +
  labs(x="Total ROH Length (Mb)",
       y="Number of ROH segments",
       shape="Individual",
       color="Subspecies") +
  scale_color_manual(values=subspecies_colors) +
  scale_shape_manual(values=individual_shapes) +
  theme_classic(base_size = 10) +
  theme(
    # --- FULL BOUNDING BOX ---
    panel.border = element_rect(color="black", fill=NA, linewidth=0.8),
    axis.line = element_blank(),
    
    # --- FACET LABEL BOXES ---
    strip.background = element_rect(color="black", fill="gray95", linewidth=0.8),
    strip.text = element_text(color="black", size=10, face="bold"),
    
    # --- LEGEND SETTINGS ---
    # Placing legend bottom-right often avoids data overlap in scatterplots
    legend.position = "right", 
    legend.background = element_rect(fill="white", color="black", linewidth=0.3),
    legend.title = element_text(size=9, face="bold"),
    legend.text = element_text(size=8),
    
    # --- AXIS SETTINGS ---
    axis.title = element_text(color="black", size=10, face="bold"),
    axis.text = element_text(color="black", size=9)
  )

# --- Save for Molecular Ecology Resources (169mm Full Width) ---
ggsave(file.path(output_dir, "ROH_Count_vs_Sum_Final.pdf"), p5, 
       width=169, height=100, units="mm", device="pdf", useDingbats=FALSE)

ggsave(file.path(output_dir, "ROH_Count_vs_Sum_Final.jpeg"), p5, 
       width=169, height=100, units="mm", dpi=300)

print("Scatterplot p5 exported successfully.")

  # ----------------------------------------------------
# 4.  Cumulative % Genome in ROH vs ROH Size (PLAIN LINE CDF) - CORRECTED
# ----------------------------------------------------
# ----------------------------------------------------
# ROH Cumulative % Genome vs ROH Size (PLAIN LINE CDF)
# ----------------------------------------------------
library(IRanges)
library(dplyr)
library(ggplot2)
library(scales)

# -----------------------------
# 1. Merge ROHs per chromosome
# -----------------------------
merge_overlaps_chr <- function(df) {
  df %>%
    group_by(CHR) %>%
    group_modify(~{
      ir <- IRanges(start = .x$POS1, end = .x$POS2)
      merged <- reduce(ir)
      tibble(ROH_Mb = width(merged) / 1e6)
    }) %>% ungroup()
}

# -----------------------------
# 2. Sample and species info
# -----------------------------
sample_info <- tibble(
  IID = c("SRR17129394","SRR17134085","SRR17134086","SRR17134087","SRR17134088"),
  Species = c("Addra","Addra","Addra","Mhorr","Mhorr")
)

# Optional: assign colors per individual (blue shades = Addra, orange shades = Mhorr)
addra_colors <- c("#1f78b4", "#6baed6", "#a6cee3")  # 3 Addra
mhorr_colors <- c("#ff7f00", "#fdbf6f")            # 2 Mhorr

individual_colors <- c(
  setNames(addra_colors, sample_info$IID[sample_info$Species=="Addra"]),
  setNames(mhorr_colors, sample_info$IID[sample_info$Species=="Mhorr"])
)

# -----------------------------
# 3. Merge ROHs genome-wide and attach species
# -----------------------------
roh_merged <- roh_segments %>%
  filter(IID %in% sample_info$IID) %>%
  group_by(IID) %>%
  group_modify(~ merge_overlaps_chr(.x)) %>%
  ungroup() %>%
  left_join(sample_info, by="IID")

# -----------------------------
# 4. Compute CDF
# -----------------------------
roh_cdf <- roh_merged %>%
  arrange(IID, ROH_Mb) %>%
  group_by(IID, Species) %>%
  mutate(
    Total_ROH_bp = sum(ROH_Mb * 1e6),
    Cumulative_bp = cumsum(ROH_Mb * 1e6),
    Percent_Cumulative = (Cumulative_bp / genome_size) * 100,
    Label = paste(IID, "(", Species, ")", sep="")   # For legend
  ) %>%
  ungroup()

# -----------------------------
# 5. Plot CDF with legend slightly left
# -----------------------------
# -----------------------------
# 4. Cumulative Distribution Plot (CDF) - Compact Legend Box
# -----------------------------

p4 <- ggplot(roh_cdf, aes(x = ROH_Mb, y = Percent_Cumulative, 
                          group = IID, color = IID)) +
  geom_line(linewidth = 1.2) +  
  scale_x_log10(labels = label_comma()) +
  scale_color_manual(
    values = individual_colors,
    labels = setNames(roh_cdf$Label[match(names(individual_colors), roh_cdf$IID)],
                      names(individual_colors)),
    name = "Individual (Sub-species)"
  ) +
  labs(
    x = "ROH Size (Mb) (log10 scale)",
    y = "Cumulative % Genome in ROH"
  ) +
  theme_classic(base_size = 12) + 
  theme(
    # --- FULL BOUNDING BOX ---
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    
    # --- AXIS TEXT ---
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    
    # --- LEGEND POSITION (Top-Left) ---
    legend.position = c(0.01, 0.99), 
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    
    # --- COMPACT LEGEND BOX SETTINGS ---
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(0.8, "cm"),   # Shorter line samples
    legend.key.height = unit(0.4, "cm"),  # Tightens vertical space
    legend.margin = margin(t = 2, r = 5, b = 2, l = 5) # Tightens padding around text
  )

# -----------------------------
# Save for Molecular Ecology Resources
# -----------------------------

ggsave(
  filename = file.path(output_dir, "ROH_CDF_Compact_Legend.pdf"), 
  plot = p4, 
  width = 169, 
  height = 110, 
  units = "mm", 
  device = "pdf",
  useDingbats = FALSE
)

ggsave(
  filename = file.path(output_dir, "ROH_CDF_Compact_Legend.jpeg"), 
  plot = p4, 
  width = 169, 
  height = 110, 
  units = "mm", 
  dpi = 300
)

print("CDF plot saved: Legend box size reduced for a compact look.")



# --------------------------------------------------------------------------------------
# 5. Calculate NROH, SROH, LROH, and correlation between NROH and SROH
# --------------------------------------------------------------------------------------
#######################################################
# 5. Calculate NROH, SROH, LROH, and Correlation
#######################################################
library(dplyr)
library(ggplot2)
library(Cairo)
library(ragg)

# -----------------------------
# 1. Calculate Summary Statistics
# -----------------------------
roh_stats <- roh_segments_3 %>%
  group_by(IID, Species) %>%
  summarise(
    NROH = n(),                  # Number of ROH
    SROH = sum(ROH_Mb),          # Sum of ROH in Mb
    LROH = mean(ROH_Mb),         # Mean length of ROH in Mb
    .groups = "drop"
  )

# Calculate Pearson Correlation
cor_result <- cor.test(roh_stats$NROH, roh_stats$SROH, method = "pearson")
r_val <- round(cor_result$estimate, 3)
p_val_formatted <- format.pval(cor_result$p.value, digits = 2)

# -----------------------------
# 2. Define Clean Ultra-Tight Theme
# -----------------------------
clean_tight_theme <- theme_classic() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "none",         # Remove legend entirely
    plot.title = element_blank(),     # Remove Title area
    plot.subtitle = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# -----------------------------
# 3. Create Plots (Sub-species Labels)
# -----------------------------

# Plot A: Correlation NROH vs SROH
p_corr_clean <- ggplot(roh_stats, aes(x = NROH, y = SROH)) +
  geom_point(aes(color = Species), size = 4) + 
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
  labs(x = "NROH", y = "SROH") +
  annotate("text", x = -Inf, y = Inf, 
           label = paste0("r = ", r_val, "\np = ", p_val_formatted), 
           hjust = -0.2, vjust = 1.2, fontface = "bold", size = 4.5) +
  scale_color_manual(values = subspecies_colors) +
  clean_tight_theme

# Plot B: NROH Boxplot
p_nroh <- ggplot(roh_stats, aes(x = Species, y = NROH, fill = Species)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3) +
  labs(x = "Sub-species", y = "NROH") +
  scale_fill_manual(values = subspecies_colors) +
  clean_tight_theme

# Plot C: SROH Boxplot
p_sroh <- ggplot(roh_stats, aes(x = Species, y = SROH, fill = Species)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3) +
  labs(x = "Sub-species", y = "SROH") +
  scale_fill_manual(values = subspecies_colors) +
  clean_tight_theme

# Plot D: LROH Boxplot
p_lroh <- ggplot(roh_stats, aes(x = Species, y = LROH, fill = Species)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3) +
  labs(x = "Sub-species", y = "LROH") +
  scale_fill_manual(values = subspecies_colors) +
  clean_tight_theme

# -----------------------------
# 4. Save Function (PDF + JPEG)
# -----------------------------
save_ultra_tight <- function(plot_obj, name, w = 5, h = 4) {
  # Save PDF using Cairo device
  ggsave(paste0(name, ".pdf"), plot_obj, width = w, height = h, device = cairo_pdf)
  
  # Save JPEG using ragg
  ragg::agg_jpeg(paste0(name, ".jpeg"), width = w, height = h, units = "in", res = 300)
  print(plot_obj)
  invisible(dev.off())
}

# Execute Saving
save_ultra_tight(p_corr_clean, "Correlation_NROH_SROH_clean")
save_ultra_tight(p_nroh, "Boxplot_NROH_clean")
save_ultra_tight(p_sroh, "Boxplot_SROH_clean")
save_ultra_tight(p_lroh, "Boxplot_LROH_clean")
```

