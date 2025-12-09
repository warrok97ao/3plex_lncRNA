library(tidyverse)
library(ggpubr)
library(data.table) 

files <- list(
  "Random_Negatives" = "/home/aleone/Random_Negatives/ALL_shape.3plex_stability.matrix",
  "cCRE_Balanced"    = "/home/aleone/cCRE_Balanced/ALL_shape.3plex_stability.matrix",
  "Biosample_Specific" = "/home/aleone/Biosample_Specific/ALL_shape.3plex_stability.matrix"
)

shape_features <- c("HelT", "MGW", "ProT", "Roll")
score_col <- "Stability_best" 
label_col <- "pos_neg"

run_analysis <- function(dataset_name, file_path) {
  message(paste("--- Processing:", dataset_name, "---"))
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
 
  df <- fread(file_path)
  
    print(paste("Columns in file:", paste(colnames(df), collapse = ", ")))

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
    mutate(Feature = str_remove(Feature, "mean_")) # Clean up names

  p <- ggplot(df_long, aes(x = mean_score, y = Shape_Value)) +
    geom_point(aes(color = Type, shape = Type), size = 3, alpha = 0.7) +
    
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE, alpha = 0.5) +
    
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    
    facet_wrap(~Feature, scales = "free_y") +
    
    scale_color_manual(values = c("Negative" = "blue", "Positive" = "red")) +
    
    labs(
      title = paste("Correlation Analysis:", dataset_name),
      subtitle = "Mean Triplex Stability vs. Mean DNA Shape per lncRNA",
      x = "Mean 3plex Stability Score",
      y = "Normalized Shape Feature Value"
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


for (name in names(files)) {
  run_analysis(name, files[[name]])
}