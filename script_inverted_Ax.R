library(tidyverse)
library(ggpubr)
library(data.table)
library(hexbin)  

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

  
  # 3. PLOTTING: Shape (X) vs Score (Y)
  p <- ggplot(df_long, aes(x = Shape_Value, y = mean_score)) +
    
    geom_point(aes(color = Type, shape = Type), size = 3, alpha = 0.7) +
    
    # Add linear regression line to visualize the trend
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE, alpha = 0.5) +
    
    # Add Pearson correlation stats (R and p-value) to answer the reviewer's specific question
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    
    # 'scales = free' allows each shape feature (HelT vs Roll) to have its own X-axis range
    facet_wrap(~Feature, scales = "free") +
    
    scale_color_manual(values = c("Negative" = "blue", "Positive" = "red")) +
    
    labs(
      title = paste("Correlation Analysis:", dataset_name),
      subtitle = "Mean DNA Shape vs. Mean Triplex Stability per lncRNA",
      x = "Normalized Shape Feature Value",  
      y = "Mean 3plex Stability Score"        
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "#f0f0f0"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )

  output_filename <- paste0("Scatterplot_", dataset_name, ".png")
  ggsave(output_filename, plot = p, width = 10, height = 8, dpi = 300)
  message(paste("Saved plot to:", output_filename))
}

# Execute the analysis for all 3 control datasets 
for (name in names(files)) {
  run_analysis(name, files[[name]])
}


###################################################
# --- PART 2: REGION-LEVEL ANALYSIS (Individual Peaks) ---
run_region_analysis <- function(dataset_name, file_path, bins = 25) {
  message(paste("--- Processing Region Plot for:", dataset_name, "---"))
  
  if (!file.exists(file_path)) return(NULL)
  
  df <- fread(file_path)

  # FILTER: Select ONLY Positive peaks
  df_pos <- df %>%
    filter(get(label_col) == "pos")

  # Compute global 75th percentile cutoff for 3plex score among positives
  score_cutoff <- quantile(df_pos[[score_col]], probs = 0.75, na.rm = TRUE)

  # Label Top 25% vs Bottom 75%
  df_pos <- df_pos %>%
    mutate(
      Group = ifelse(.data[[score_col]] >= score_cutoff, "Top25", "Bottom75"),
      Group = factor(Group, levels = c("Bottom75", "Top25"))
    )

  df_long <- df_pos %>%
    pivot_longer(
      cols = c(HelT, MGW, ProT, Roll),
      names_to = "Feature",
      values_to = "Shape_Value"
    )

  # --- PLOT A: SCATTER WITH ALL POINTS, COLORED BY TOP 25% vs REST ---
  p_overall_points <- ggplot(df_long, aes(x = Shape_Value, y = .data[[score_col]], color = Group)) +
    geom_point(alpha = 0.1, size = 1) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed") +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", color = "black") +
    facet_wrap(~Feature, scales = "free") +
    scale_color_manual(
      values = c("Bottom75" = "gray80", "Top25" = "gold"),
      name   = "3plex score group"
    ) +
    labs(
      title = paste("Region-Level Correlation (Points):", dataset_name),
      subtitle = "Top 25% 3plex score in yellow; remaining 75% in gray (Positives Only)",
      x = "DNA Shape Value",
      y = "3plex Stability Score"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  ggsave(
    paste0("RegionPlot_OVERALL_POINTS_", dataset_name, ".png"),
    p_overall_points, width = 10, height = 8, dpi = 300
  )

  # --- PLOT B: HEXBIN VERSION ---
  p_overall_hex <- ggplot(df_long, aes(x = Shape_Value, y = .data[[score_col]])) +
    stat_bin_hex(aes(fill = ..count..), binwidth = c(0.1, 0.1), color = "white", alpha = 0.9) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
    stat_cor(
      method = "pearson",
      label.x.npc = "left",
      label.y.npc = "top",
      color = "black",
      size = 4
    ) +
    scale_fill_viridis_c(
      option = "C",
      trans = "log10",
      name = "Peak count"
    ) +
    facet_wrap(~Feature, scales = "free") +
    labs(
      title = paste("Region-Level Correlation (Hexbin):", dataset_name),
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

  # --- PLOT C: STRATIFICATION BY GENE ---
  p_faceted <- ggplot(df_long, aes(x = Shape_Value, y = .data[[score_col]])) +
    geom_point(alpha = 0.2, size = 0.5, color = "red") +
    geom_smooth(method = "lm", color = "black", se = FALSE, size = 0.5) +
    facet_grid(lncRNA ~ Feature, scales = "free") +
    theme_bw() +
    theme(
      strip.text.y = element_text(angle = 0, size = 7),
      axis.text = element_text(size = 6)
    ) +
    labs(title = paste("Per-LncRNA Correlation:", dataset_name))

  ggsave(
    paste0("RegionPlot_FACETED_", dataset_name, ".png"),
    p_faceted, width = 12, height = 30, limitsize = FALSE
  )
}

# Run Region Analysis ONLY on one file
run_region_analysis("Random_Negatives", files[["Random_Negatives"]])
