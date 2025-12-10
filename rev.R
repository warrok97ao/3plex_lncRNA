library(tidyverse)
library(ggpubr)
library(data.table)
library(ggrepel) # Included to label the specific lncRNAs

# --- CONFIGURATION ---
files <- list(
  "Random_Negatives" = "Random_Negatives/ALL_shape.3plex_stability.matrix",
  "cCRE_Balanced"    = "cCRE_Balanced/ALL_shape.3plex_stability.matrix",
  "Biosample_Specific" = "Biosample_Specific/ALL_shape.3plex_stability.matrix"
)

score_col <- "Stability_best" 
label_col <- "pos_neg"

run_reviewer_response_plot <- function(dataset_name, file_path) {
  message(paste("--- Processing Reviewer Plot for:", dataset_name, "---"))
  
  if (!file.exists(file_path)) return(NULL)
  df <- fread(file_path)

  # 1. AGGREGATE PER lncRNA & LABEL
  # We calculate the mean shape and score for Pos and Neg separately
  df_summary <- df %>%
    group_by(lncRNA, label = get(label_col)) %>%
    summarise(
      mean_score = mean(get(score_col), na.rm = TRUE),
      mean_HelT  = mean(HelT, na.rm = TRUE),
      mean_MGW   = mean(MGW, na.rm = TRUE),
      mean_ProT  = mean(ProT, na.rm = TRUE),
      mean_Roll  = mean(Roll, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(label = ifelse(label == 1 | label == "pos", "pos", "neg"))

  # 2. PIVOT TO WIDE FORMAT
  # We need "Pos" and "Neg" values in the same row to calculate the difference
  df_wide <- df_summary %>%
    pivot_wider(
      id_cols = lncRNA,
      names_from = label,
      values_from = c(mean_score, mean_HelT, mean_MGW, mean_ProT, mean_Roll)
    ) %>%
    mutate(
      # METRIC X: Hoogsteen Capability (Mean Score of POSITIVE sites)
      Hoogsteen_Ability = mean_score_pos,
      
      # METRIC Y: Strength of Difference (Pos - Neg)
      # We keep the sign to see if the direction is consistent
      Diff_HelT = mean_HelT_pos - mean_HelT_neg,
      Diff_MGW  = mean_MGW_pos  - mean_MGW_neg,
      Diff_ProT = mean_ProT_pos - mean_ProT_neg,
      Diff_Roll = mean_Roll_pos - mean_Roll_neg
    )

  # 3. RESHAPE FOR PLOTTING
  # Convert the "Diff" columns back to long format for faceting
  df_plot <- df_wide %>%
    pivot_longer(
      cols = starts_with("Diff_"),
      names_to = "Feature",
      values_to = "Shape_Difference"
    ) %>%
    mutate(Feature = str_remove(Feature, "Diff_"))

  # 4. PLOTTING (X and Y inverted to match other plots)
  p <- ggplot(df_plot, aes(x = Hoogsteen_Ability, y = Shape_Difference)) +
    
    # Add zero line (No difference between Pos and Neg)
    geom_hline(yintercept = 0, linetype = "solid", color = "gray50", linewidth = 0.5) +
    
    # Points representing each lncRNA
    geom_point(size = 3, alpha = 0.7, color = "#2166AC") +
    
    # Linear Regression Line to show the trend
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE, linewidth = 0.8) +
    
    # Add Labels for lncRNAs (Crucial for answering "which lncRNA?")
    geom_text_repel(
      aes(label = lncRNA), 
      size = 2.5, 
      max.overlaps = 15,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "gray50",
      segment.size = 0.3
    ) +
    
    # Correlation Stats
    stat_cor(
      method = "pearson", 
      label.x.npc = 1,     
      label.y.npc = 1,     
      hjust = 1,           
      vjust = 0.3,           
      size = 3.5,
      label.sep = "\n",
      geom = "label",
      fill = "white",
      color = "black",
      alpha = 0.95,
      label.padding = unit(0.5, "lines"),
      label.r = unit(0.15, "lines")
    ) +
    
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    
    facet_wrap(~Feature, scales = "free_y") +
    
    labs(
      title = paste("Reviewer Q2 Analysis:", dataset_name),
      subtitle = "Correlation between Hoogsteen Capability and Strength of Shape Preference",
      x = "Predicted Hoogsteen Ability\n(Mean 3plex Stability Score of Positives)",
      y = "Strength of Shape Preference\n(Mean Positive - Mean Negative)"
    ) +
    
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "#f0f0f0"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      panel.grid.minor = element_blank()
    )

  # Save as PDF for publication
  output_filename_pdf <- paste0("Reviewer_Q2_Correlation_", dataset_name, ".pdf")
  ggsave(output_filename_pdf, plot = p, width = 12, height = 9, device = "pdf")
  message(paste("Saved PDF plot to:", output_filename_pdf))
  
  # Optional PNG backup
  output_filename_png <- paste0("Reviewer_Q2_Correlation_", dataset_name, ".png")
  ggsave(output_filename_png, plot = p, width = 12, height = 9, dpi = 600)
  message(paste("Saved PNG plot to:", output_filename_png))
}

# Run for all datasets
for (name in names(files)) {
  run_reviewer_response_plot(name, files[[name]])
}


library(tidyverse)
library(data.table)
library(ggpubr)

# --- CONFIGURATION ---
# Use the dataset you find most representative (e.g., Random or cCRE)
file_path <- "Random_Negatives/ALL_shape.3plex_stability.matrix" 
dataset_name <- "Random_Negatives"

run_binned_analysis <- function(dataset_name, file_path) {
  if (!file.exists(file_path)) return(NULL)
  
  df <- fread(file_path)
  
  # 1. FILTER FOR POSITIVES ONLY
  # We are asking: "Among triplexes, does better score mean specific shape?"
  df_pos <- df %>% 
    filter(pos_neg == "pos" | pos_neg == 1)
  
  # 2. CREATE STABILITY BINS (QUANTILES)
  # We divide peaks into 4 groups based on Stability_best score
  df_binned <- df_pos %>%
    mutate(
      Score_Bin = ntile(Stability_best, 4), # Divide into 4 Quartiles
      Bin_Label = case_when(
        Score_Bin == 1 ~ "Q1 (Low Stability)",
        Score_Bin == 2 ~ "Q2 (Med Stability)",
        Score_Bin == 3 ~ "Q3 (High Stability)",
        Score_Bin == 4 ~ "Q4 (Best Stability)"
      )
    )
  
  # 3. PIVOT LONG FOR PLOTTING
  df_long <- df_binned %>%
    pivot_longer(
      cols = c(HelT, MGW, ProT, Roll),
      names_to = "Feature",
      values_to = "Shape_Value"
    )

  # 4. PLOT: BOXPLOT OF SHAPE BY STABILITY BIN
  p <- ggplot(df_long, aes(x = Bin_Label, y = Shape_Value, fill = Bin_Label)) +
    
    # Boxplot shows the distribution and consistency (tightness)
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    
    # Add trend line through the MEDIANS of the boxes
    stat_summary(fun = median, geom = "line", aes(group = 1), color = "black", size = 1, linetype = "dashed") +
    stat_summary(fun = median, geom = "point", color = "black", size = 2) +
    
    stat_compare_means(label = "p.signif", method = "kruskal.test", label.y.npc = 0.95) +
    
    facet_wrap(~Feature, scales = "free_y") +
    
    scale_fill_brewer(palette = "OrRd") + # Sequential Red palette
    
    labs(
      title = paste("Dose-Response Analysis:", dataset_name),
      subtitle = "Does DNA Shape change as Triplex Stability increases?",
      x = "Triplex Stability Quartile (Hoogsteen Capability)",
      y = "DNA Shape Value"
    ) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  ggsave(paste0("Binned_Analysis_", dataset_name, ".pdf"), p, width = 10, height = 8)
  message(paste("Saved Binned Analysis for", dataset_name))
}

run_binned_analysis(dataset_name, file_path)

