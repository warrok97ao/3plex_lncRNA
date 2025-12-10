library(tidyverse)
library(ggpubr)
library(data.table)
library(hexbin)  
library(ggpointdensity) 
# --- CONFIGURATION ---
#  three different datasets to test the hypothesis against increasing levels of stringency:
# 1. Random: Baseline comparison against random genome.
# 2. cCRE: Controls for general open chromatin (rules out accessibility bias).
# 3. Biosample: Controls for cell-type specific chromatin (strictest control).

files <- list(
  "Random_Negatives" = "Random_Negatives/ALL_shape.3plex_stability.matrix",
  "cCRE_Balanced"    = "cCRE_Balanced/ALL_shape.3plex_stability.matrix",
  "Biosample_Specific" = "Biosample_Specific/ALL_shape.3plex_stability.matrix"
)

shape_features <- c("HelT", "MGW", "ProT", "Roll")
score_col <- "Stability_best" 
label_col <- "pos_neg"

# --- PART 1: AGGREGATED ANALYSIS (Mean values per lncRNA) ---
run_analysis <- function(dataset_name, file_path) {
  message(paste("--- Processing:", dataset_name, "---"))
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
 
  # Use fread for fast reading of large genomic matrix files
  df <- fread(file_path)
  
  print(paste("Columns in file:", paste(colnames(df), collapse = ", ")))

  # 1. AGGREGATION: Calculate the MEAN Score and MEAN Shape for each lncRNA.
  # This treats each lncRNA as a biological unit (centroid), reducing noise from individual peaks.
  df_summary <- df %>%
    group_by(lncRNA, get(label_col)) %>% 
    summarise(
      mean_score = mean(get(score_col), na.rm = TRUE),
      mean_HelT = mean(HelT, na.rm = TRUE),
      mean_MGW  = mean(MGW, na.rm = TRUE),
      mean_ProT = mean(ProT, na.rm = TRUE),
      mean_Roll = mean(Roll, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    rename(label = `get(label_col)`) %>%
    mutate(
      Type = ifelse(label == 1 | label == "pos", "Positive", "Negative"),
      Type = factor(Type, levels = c("Negative", "Positive")) 
    )

  df_long <- df_summary %>%
    pivot_longer(
      cols = c(mean_HelT, mean_MGW, mean_ProT, mean_Roll),
      names_to = "Feature",
      values_to = "Shape_Value"
    ) %>%
    mutate(Feature = str_remove(Feature, "mean_")) 

  
# --- 3. PLOTTING: Shape (X) vs Score (Y) ---
 p <- ggplot(df_long, aes(x = Shape_Value, y = mean_score)) +
  
  geom_point(aes(color = Type, shape = Type), size = 3, alpha = 0.7) +
  
  geom_smooth(
    data = subset(df_long, Type == "Positive"),
    method = "lm",
    color = "black",
    linetype = "dashed",
    se = FALSE
  ) +

  stat_cor(
    data = subset(df_long, Type == "Positive"),
    method = "pearson",
    label.x.npc = 1,     # Align to far left edge
      label.y.npc = 0.98,     # Align to very top edge
      hjust = 1,           # Left-justify text
      vjust = 0.25,           # Top-justify box (makes it hang DOWN)
    size = 3,
    label.sep = "\n",
    geom = "label",
    fill = "white",
    alpha = 0.95,
    label.padding = unit(0.5, "lines"),
    label.r = unit(0.15, "lines")
  ) +

  # ---- ADD THIS LINE ----
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +

  facet_wrap(~Feature, scales = "free") +

  scale_color_manual(values = c("Negative" = "blue", "Positive" = "red")) +
  
  labs(
    title = paste("Correlation Analysis:", dataset_name),
    subtitle = "Mean DNA Shape vs. Mean Triplex Stability (Positives Only)",
    x = "Normalized Shape Feature Value",
    y = "Mean 3plex Stability Score"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.text = element_text(size = 10)
  )

  # Save as high-resolution PDF for publication
  output_filename_pdf <- paste0("Scatterplot_", dataset_name, ".pdf")
  ggsave(output_filename_pdf, plot = p, width = 10, height = 8, device = "pdf")
  message(paste("Saved PDF plot to:", output_filename_pdf))
  
  # Also save PNG version (optional, for presentations/slides)
  output_filename_png <- paste0("Scatterplot_", dataset_name, ".png")
  ggsave(output_filename_png, plot = p, width = 10, height = 8, dpi = 600)
  message(paste("Saved PNG plot to:", output_filename_png))
}

# Execute the analysis for all 3 control datasets 
for (name in names(files)) {
  run_analysis(name, files[[name]])
}


# --- PART 2: REGION-LEVEL ANALYSIS (Individual Peaks) ---
run_region_analysis <- function(dataset_name, file_path, bins = 25, top_percent = 0.1) {
  message(paste0("--- Processing Region Plot for: ", dataset_name, " (Top ", top_percent * 100, "%) ---"))
  if (!file.exists(file_path)) return(NULL)
  df <- fread(file_path)

  # FILTER: Select ONLY Positive peaks
  df_pos <- df %>% filter(get(label_col) == "pos")

  # (Optional) GLOBAL CUTOFF for PLOT A/B (unchanged)
  score_cutoff <- quantile(df_pos[[score_col]], probs = 1 - top_percent, na.rm = TRUE)
  df_pos <- df_pos %>%
    mutate(
      Group_global = ifelse(.data[[score_col]] >= score_cutoff, "Top", "Rest"),
      Group_global = factor(Group_global, levels = c("Rest", "Top"))
    )

  # LONG FORMAT for PLOT A/B (global top x%)
  df_long <- df_pos %>%
    pivot_longer(
      cols = c(HelT, MGW, ProT, Roll),
      names_to = "Feature",
      values_to = "Shape_Value"
    )

  # --- PLOT A: SCATTER WITH ALL POINTS, COLORED BY GLOBAL TOP X% vs REST ---
  p_overall_points <- ggplot(df_long, aes(x = Shape_Value, y = .data[[score_col]], color = Group_global)) +
    geom_point(alpha = 0.1, size = 1) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed") +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", color = "black") +
    facet_wrap(~Feature, scales = "free") +
    scale_color_manual(
      values = c("Rest" = "gray80", "Top" = "gold"),
      labels = c(
        Rest = paste0("Bottom ", 100 - top_percent * 100, "%"),
        Top  = paste0("Top ", top_percent * 100, "%")
      ),
      name = "3plex score group"
    ) +
    labs(
      title    = paste("Region-Level Correlation (Points):", dataset_name),
      subtitle = paste0(
        "Top ", top_percent * 100, "% 3plex score in yellow; remaining ",
        100 - top_percent * 100, "% in gray (Positives Only)"
      ),
      x = "DNA Shape Value",
      y = "3plex Stability Score"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(
    paste0("RegionPlot_OVERALL_POINTS_", dataset_name, "_Top", top_percent * 100, ".png"),
    p_overall_points, width = 10, height = 8, dpi = 300
  )

  # --- PLOT B: HEXBIN VERSION (unchanged except using Group_global if needed) ---
  p_overall_hex <- ggplot(df_long, aes(x = Shape_Value, y = .data[[score_col]])) +
    stat_bin_hex(aes(fill = ..count..), bins = bins, alpha = 0.9) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
    stat_cor(
      method = "pearson", label.x.npc = "left",
      label.y.npc = "top", color = "black", size = 4
    ) +
    scale_fill_viridis_c(
      option = "C", trans = "log10", name = "Peak count"
    ) +
    facet_wrap(~Feature, scales = "free") +
    labs(
      title    = paste("Region-Level Correlation (Hexbin):", dataset_name),
      subtitle = "Each hexagon summarizes multiple positive genomic peaks",
      x = "DNA Shape Value",
      y = "3plex Stability Score",
      fill = "Peak count"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(
    paste0("RegionPlot_OVERALL_HEX_", dataset_name, ".png"),
    p_overall_hex, width = 10, height = 8, dpi = 300
  )
# --- PLOT B: POINTS COLORED BY DENSITY (Professor's Request) ---
  p_overall_density <- ggplot(df_long, aes(x = Shape_Value, y = .data[[score_col]])) +
    
    # geom_pointdensity calculates the density of neighbors for each point
    geom_pointdensity(adjust = 4, size = 2) + 
    
    # Add the trend line on top
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
    
    # Pearson correlation stats
    stat_cor(
      method = "pearson", 
      label.x.npc = "left", 
      label.y.npc = "top", 
      color = "black", 
      size = 4
    ) +
    
    # Use the 'viridis' color scale (same as your example image)
    # option = "C" is 'plasma' (warm colors), "D" is standard viridis (purple-green-yellow)
    scale_color_viridis_c(
      option = "D",   # Change to "C" if you want the orange/purple look
      name = "Density"
    ) +
    
    facet_wrap(~Feature, scales = "free") +
    
    labs(
      title = paste("Region-Level Correlation (Density):", dataset_name),
      subtitle = "Points colored by local density (lighter color = more overlapping points)",
      x = "DNA Shape Value",
      y = "3plex Stability Score"
    ) +
    theme_bw() +
    theme(legend.position = "right")

  ggsave(
    paste0("RegionPlot_OVERALL_DENSITY_", dataset_name,".png"), 
    p_overall_density, 
    width = 10, 
    height = 8, 
    dpi = 300
  )
  # --- NEW: PER-LncRNA TOP X% GROUP FOR FACETED PLOT ---

  df_pos_lncRNA <- df %>%
    filter(get(label_col) == "pos") %>%
    group_by(lncRNA) %>%
    mutate(
      score_cutoff_lncRNA = quantile(.data[[score_col]], probs = 1 - top_percent, na.rm = TRUE),
      Group_lncRNA = ifelse(.data[[score_col]] >= score_cutoff_lncRNA, "Top", "Rest")
    ) %>%
    ungroup() %>%
    mutate(Group_lncRNA = factor(Group_lncRNA, levels = c("Rest", "Top")))

  df_long_lncRNA <- df_pos_lncRNA %>%
    pivot_longer(
      cols = c(HelT, MGW, ProT, Roll),
      names_to = "Feature",
      values_to = "Shape_Value"
    )

  # --- PLOT C: STRATIFICATION BY GENE WITH PER-LncRNA TOP X% COLOURING ---
  p_faceted <- ggplot(df_long_lncRNA,
                      aes(x = Shape_Value,
                          y = .data[[score_col]],
                          color = Group_lncRNA)) +
    geom_point(alpha = 0.4, size = 0.5) +
    geom_smooth(method = "lm", color = "black", se = FALSE, size = 0.5) +
    facet_grid(lncRNA ~ Feature, scales = "free") +
    scale_color_manual(
      values = c("Rest" = "gray80", "Top" = "gold"),
      labels = c(
        Rest = paste0("Bottom ", 100 - top_percent * 100, "%"),
        Top  = paste0("Top ", top_percent * 100, "%")
      ),
      name = "3plex score group (per lncRNA)"
    ) +
    theme_bw() +
    theme(
      strip.text.y = element_text(angle = 0, size = 7),
      axis.text    = element_text(size = 6),
      legend.position = "bottom"
    ) +
    labs(
      title = paste("Per-LncRNA Correlation with Per-LncRNA Top", top_percent * 100, "% Highlighted:", dataset_name),
      x = "DNA Shape Value",
      y = "3plex Stability Score"
    )

  ggsave(
    paste0("RegionPlot_FACETED_", dataset_name, "_Top", top_percent * 100, ".png"),
    p_faceted, width = 12, height = 30, limitsize = FALSE, dpi = 300
  )
}
run_region_analysis("Random_Negatives", files[["Random_Negatives"]])
