#### Dxy vs Fst plane

```bash
#### Fst vs. Dxy plane for species delimitation.


# Load libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

setwd("F:/Collaborative_Projects/Dama_Gazelle_Project/Pixy/Outgroup_pxy_dxy/Final/Klaus_Plot")

# --- 1. PARAMETERS & CURVE GENERATION ---
# Mutation rate based on Dama Gazelle data
mu <- 2.96e-9

# Function to generate theoretical curves for speciation plane
calc_theory <- function(Ne, mu, max_gen = 2.5e6) {
  gens <- seq(0, max_gen, by = 5000)
  pi_anc <- 4 * Ne * mu
  dxy <- pi_anc + (2 * mu * gens)
  fst <- (2 * mu * gens) / dxy
  data.frame(dxy = dxy * 100, fst = fst, mga = gens / 1e6)
}

# Generate the theoretical "envelope" using Ne bounds (10k and 70k)
# Based on parameters in de Jong et al. preprint
curve_low <- calc_theory(Ne = 30000, mu = mu)
curve_high <- calc_theory(Ne = 150000, mu = mu)

# --- 2. DATA PROCESSING ---
# Assumes files dama_fst.txt and dama_dxy.txt are in your working directory
# Check file paths before running
fst_raw <- read.delim("dama_fst.txt")
dxy_raw <- read.delim("dama_dxy.txt")

plot_data <- fst_raw %>%
  group_by(pop1, pop2) %>%
  summarise(mean_fst = max(0, mean(avg_wc_fst, na.rm = TRUE)), .groups = 'drop') %>%
  inner_join(
    dxy_raw %>%
      group_by(pop1, pop2) %>%
      summarise(mean_dxy = mean(avg_dxy, na.rm = TRUE) * 100, .groups = 'drop'),
    by = c("pop1", "pop2")
  ) %>%
  mutate(
    label = paste(pop1, "vs", pop2),
    # Assigning categories for coloring based on taxonomic ranks
    rank = case_when(
      pop1 == "Addra" & pop2 == "Mhorr" ~ "Within Sub-Species",
      grepl("Grant", label) & !grepl("Thompson", label) ~ "Between Species",
      grepl("Thompson", label) ~ "Between Genera",
      TRUE ~ "Other"
    )
  )

# --- 3. CREATE THE FINAL ULTRA-BOLD, LEGEND-LOWERED PLOT ---
speciation_plane <- ggplot() +
  # 1. Theoretical Curves: Ultra-Thick (size = 3)
  geom_line(data = curve_low, aes(x = dxy, y = fst), linetype = "dashed", color = "grey60", size = 3) +
  geom_line(data = curve_high, aes(x = dxy, y = fst), linetype = "dashed", color = "grey60", size = 3) +
  
  # 2. Time Labels (Mga): Large numbers along curves
  geom_text(data = filter(curve_low, mga %in% c(0.02,0.05, 0.1, 0.2, 0.4, 0.8, 1.5)), 
            aes(x = dxy, y = fst, label = mga), 
            size = 9, vjust = -1.8, color = "black", fontface = "bold") +
  
  # 3. Decision Thresholds: Red Dotted lines (size = 3)
  # Based on de Jong et al.: Dxy = 0.225% and Fst = 0.26
  geom_vline(xintercept = 0.225, linetype = "dotted", color = "#E41A1C", size = 3) +
  geom_hline(yintercept = 0.26, linetype = "dotted", color = "#E41A1C", size = 3) +
  
  # 4. Data Points: Massive (size = 16) with heavy stroke
  geom_point(data = plot_data, aes(x = mean_dxy, y = mean_fst, fill = rank), 
             shape = 21, size = 16, stroke = 3.5, color = "black") +
  
  # 5. Internal Labels: Bold italic and large, repelled to avoid overlap
  geom_text_repel(data = plot_data, aes(x = mean_dxy, y = mean_fst, label = label), 
                  size = 10, fontface = "bold.italic", point.padding = 1.2, 
                  box.padding = 1.5, segment.size = 2, segment.color = 'black') +
  
  # Aesthetics and Scales
  scale_fill_manual(values = c("Within Sub-Species" = "#377EB8", 
                               "Between Species" = "#4DAF4A", 
                               "Between Genera" = "#E41A1C")) +
  scale_x_continuous(limits = c(0, 2.1), breaks = seq(0, 2, 0.5), expand = c(0.01, 0.01)) +
  # Y-axis capped at 1.0
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
  
  # --- THEME: MAXIMUM IMPACT ---
  theme_bw() + 
  theme(
    # Bounding Box: Extra Thick
    panel.border = element_rect(colour = "black", fill = NA, size = 6), 
    panel.grid = element_blank(),
    
    # Axis Ticks: Extra Thick
    axis.ticks = element_line(size = 3, color = "black"),
    axis.ticks.length = unit(0.6, "cm"),
    
    # Axis Text (Numbers): Extremely Large
    axis.text.x = element_text(size = 32, face = "bold", color = "black"),
    axis.text.y = element_text(size = 32, face = "bold", color = "black"),
    
    # Axis Titles: Huge, Bold (Not Italic)
    axis.title.x = element_text(size = 40, face = "bold", margin = margin(t = 25)),
    axis.title.y = element_text(size = 40, face = "bold", margin = margin(r = 25)),
    
    # --- FINAL LEGEND POSITION ---
    # Moved down even further, centered around Fst 0.60
    legend.position = c(0.82, 0.60), 
    legend.background = element_blank(),
    legend.title = element_text(size = 26, face = "bold"),
    legend.text = element_text(size = 24, face = "bold"),
    legend.key.size = unit(1.8, "cm"),
    
    # Main Title
    plot.title = element_text(size = 44, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 22, face = "italic", hjust = 0.5, color = "grey20", margin = margin(b = 25))
  ) +
  labs(
    # Expression for mathematical formatting (Dxy % and Fst)
    x = expression(bold(D)[xy]~("%")), 
    y = expression(bold(F)[ST]), 
    title = "DAMA GAZELLE SPECIATION PLANE",
    fill = "Taxonomic Status:"
  )

# --- 4. EXPORT TO HIGH RES ---
# Large output dimensions prevent label overlap
ggsave("Dama_Gazelle_Speciation_Plane_FINAL_ULTRA_V5.jpeg", 
       plot = speciation_plane, 
       width = 20, 
       height = 16, 
       dpi = 300)
```

