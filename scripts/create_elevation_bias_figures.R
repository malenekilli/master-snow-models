# Create elevation bias analysis figures
# Visualizes the relationship between elevation differences and model errors

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "scripts/create_elevation_bias_figures.R"
root <- normalizePath(file.path(dirname(script_path), ".."))
processed <- file.path(root, "data", "processed")

input_csv <- file.path(processed, "elevation_analysis.csv")
output_fig1 <- file.path(root, "figures_thesis", "Figure_Elevation_Bias_Scatter.png")
output_fig2 <- file.path(root, "figures_thesis", "Figure_Elevation_Bias_Boxplot.png")
output_fig3 <- file.path(root, "figures_thesis", "Figure_Elevation_Comparison_Table.png")

# Load data
elev_data <- read.csv(input_csv, stringsAsFactors = FALSE, fileEncoding = "UTF-8")

# Create figures directory if it doesn't exist
if (!dir.exists(file.path(root, "figures_thesis"))) {
  dir.create(file.path(root, "figures_thesis"), showWarnings = FALSE, recursive = TRUE)
}

# Define colors and styles
model_colors <- c("FSM2" = "#3182bd", "seNorge" = "#e6550d", "SNOWPACK" = "#31a354")
high_elev_color <- "#d62728"
low_elev_color <- "#1f77b4"

# Extract metrics
n_stations <- nrow(elev_data)

# ============================================================================
# Figure 1: Scatter plots - |Δz| vs RMSE for each model
# ============================================================================
png(output_fig1, width = 1200, height = 400, res = 100, pointsize = 12)

par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 2), cex = 0.95)

# High elevation threshold
threshold <- 100

models <- c("FSM2", "seNorge", "SNOWPACK")

for (i in 1:length(models)) {
  model <- models[i]
  rmse_col <- paste0(model, "_RMSE")
  
  # Separate by threshold
  high_idx <- elev_data$Abs_Delta_Z > threshold
  low_idx <- elev_data$Abs_Delta_Z <= threshold
  
  # Create scatter plot
  plot(elev_data$Abs_Delta_Z, elev_data[[rmse_col]],
       xlab = "|Δz| (m)",
       ylab = "RMSE (cm)",
       main = model,
       cex.main = 1.2,
       pch = 21,
       bg = ifelse(high_idx, high_elev_color, low_elev_color),
       col = "black",
       cex = 1.5,
       xlim = c(0, 800),
       ylim = c(0.5, 2.5))
  
  # Add threshold line
  abline(v = threshold, lty = 2, col = "#999999", lwd = 1.5)
  
  # Add correlation text
  corr_val <- cor(elev_data$Abs_Delta_Z, elev_data[[rmse_col]], use = "complete.obs")
  usr <- par("usr")
  text(usr[2] - 0.02 * (usr[2] - usr[1]), usr[4] - 0.08 * (usr[4] - usr[3]), 
       paste0("r = ", sprintf("%.3f", corr_val)),
       cex = 1.0,
       adj = c(1, 1),
       bg = "white")
  
  # Add legend (only on first plot)
  if (i == 1) {
    legend("topleft",
           legend = c("|Δz| ≤ 100m", "|Δz| > 100m"),
           pch = 21,
           pt.bg = c(low_elev_color, high_elev_color),
           col = "black",
           cex = 0.9,
           bty = "o",
           bg = "white")
  }
  
  grid(lty = 3, col = "#cccccc")
}

dev.off()
cat("Saved:", output_fig1, "\n")

# ============================================================================
# Figure 2: Box plots - RMSE comparison for high vs low elevation bias
# ============================================================================
png(output_fig2, width = 1200, height = 450, res = 100, pointsize = 12)

par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 2), cex = 0.95)

for (i in 1:length(models)) {
  model <- models[i]
  rmse_col <- paste0(model, "_RMSE")
  
  # Separate by threshold
  high_idx <- elev_data$Abs_Delta_Z > threshold
  low_idx <- elev_data$Abs_Delta_Z <= threshold
  
  rmse_high <- elev_data[[rmse_col]][high_idx]
  rmse_low <- elev_data[[rmse_col]][low_idx]
  
  # Calculate statistics for legend
  mean_high <- mean(rmse_high, na.rm = TRUE)
  mean_low <- mean(rmse_low, na.rm = TRUE)
  
  # Create box plot
  boxplot(list(rmse_low, rmse_high),
          names = c(sprintf("|Δz| ≤ 100m\n(n=%d)", length(rmse_low)),
                    sprintf("|Δz| > 100m\n(n=%d)", length(rmse_high))),
          col = c(low_elev_color, high_elev_color),
          border = "#333333",
          outline = TRUE,
          ylab = "RMSE (cm)",
          main = model,
          cex.main = 1.2,
          ylim = c(0.5, 2.5),
          frame.plot = TRUE)
  
  # Perform t-test (only if both groups have samples)
  if (length(rmse_high) > 0 && length(rmse_low) > 0 && length(rmse_high) > 1 && length(rmse_low) > 1) {
    t_test <- t.test(rmse_high, rmse_low)
    p_val <- t_test$p.value
  } else {
    p_val <- NA
  }
  
  # Add text box with statistics
  p_text <- if (is.na(p_val)) "p-value: N/A" else sprintf("p-value: %.3f", p_val)
  text_str <- sprintf("Mean low: %.3f\nMean high: %.3f\n%s",
                      mean_low, mean_high, p_text)
  
  # Place text in upper right corner using plot coordinates
  usr <- par("usr")
  text(usr[2] - 0.02 * (usr[2] - usr[1]), usr[4] - 0.1 * (usr[4] - usr[3]),
       text_str,
       cex = 0.85,
       adj = c(1, 1),
       bg = "white")
  
  grid(lty = 3, col = "#cccccc")
}

dev.off()
cat("Saved:", output_fig2, "\n")

# ============================================================================
# Figure 3: Summary table as image
# ============================================================================
png(output_fig3, width = 1000, height = 600, res = 100, pointsize = 11)

# Create summary statistics
summary_stats <- data.frame(
  Model = models,
  Low_Δz_RMSE = c(
    mean(elev_data$FSM2_RMSE[elev_data$Abs_Delta_Z <= threshold], na.rm = TRUE),
    mean(elev_data$seNorge_RMSE[elev_data$Abs_Delta_Z <= threshold], na.rm = TRUE),
    mean(elev_data$SNOWPACK_RMSE[elev_data$Abs_Delta_Z <= threshold], na.rm = TRUE)
  ),
  High_Δz_RMSE = c(
    mean(elev_data$FSM2_RMSE[elev_data$Abs_Delta_Z > threshold], na.rm = TRUE),
    mean(elev_data$seNorge_RMSE[elev_data$Abs_Delta_Z > threshold], na.rm = TRUE),
    mean(elev_data$SNOWPACK_RMSE[elev_data$Abs_Delta_Z > threshold], na.rm = TRUE)
  ),
  Difference = NA,
  Percent_Increase = NA,
  Correlation_r = NA,
  P_value = NA
)

# Calculate derived columns
for (i in 1:3) {
  model <- models[i]
  rmse_col <- paste0(model, "_RMSE")
  
  low_rmse <- mean(elev_data[[rmse_col]][elev_data$Abs_Delta_Z <= threshold], na.rm = TRUE)
  high_rmse <- mean(elev_data[[rmse_col]][elev_data$Abs_Delta_Z > threshold], na.rm = TRUE)
  
  summary_stats$Difference[i] <- high_rmse - low_rmse
  summary_stats$Percent_Increase[i] <- (summary_stats$Difference[i] / low_rmse * 100)
  
  corr <- cor(elev_data$Abs_Delta_Z, elev_data[[rmse_col]], use = "complete.obs")
  summary_stats$Correlation_r[i] <- corr
  
  corr_test <- cor.test(elev_data$Abs_Delta_Z, elev_data[[rmse_col]])
  summary_stats$P_value[i] <- corr_test$p.value
}

# Format table for display
table_display <- data.frame(
  Model = summary_stats$Model,
  `Low Δz\nRMSE (cm)` = sprintf("%.3f", summary_stats$Low_Δz_RMSE),
  `High Δz\nRMSE (cm)` = sprintf("%.3f", summary_stats$High_Δz_RMSE),
  `Difference\n(cm)` = sprintf("%.3f", summary_stats$Difference),
  `Increase\n(%)` = sprintf("%.1f", summary_stats$Percent_Increase),
  `Correlation\nr` = sprintf("%.3f", summary_stats$Correlation_r),
  `P-value` = sprintf("%.3f", summary_stats$P_value),
  check.names = FALSE
)

# Create table plot
par(mar = c(1, 1, 4, 1))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

title("Influence of Elevation Difference on Model RMSE", cex.main = 1.3, font.main = 2)

# Create table
x_positions <- c(0.05, 0.2, 0.35, 0.50, 0.65, 0.78, 0.92)
y_start <- 0.85

# Header
header_text <- c("Model", "Low Δz\nRMSE (cm)", "High Δz\nRMSE (cm)", "Difference\n(cm)", 
                 "Increase\n(%)", "Correlation\nr", "P-value")

for (j in 1:ncol(table_display)) {
  text(x_positions[j], y_start, header_text[j], 
       cex = 1.0, font = 2, adj = c(0.5, 0.5))
}

# Draw header line
lines(c(0.02, 0.98), c(y_start - 0.05, y_start - 0.05), lwd = 2)

# Data rows
row_height <- 0.12
colors_bg <- c("#e8f4f8", "#fff8e8", "#e8f8e8")

for (i in 1:nrow(table_display)) {
  y_pos <- y_start - 0.05 - i * row_height
  
  # Background color
  rect(0.02, y_pos - row_height/2 + 0.02, 0.98, y_pos + row_height/2 - 0.02,
       col = colors_bg[i], border = NA)
  
  # Row data
  row_data <- c(
    table_display[i, 1],
    table_display[i, 2],
    table_display[i, 3],
    table_display[i, 4],
    table_display[i, 5],
    table_display[i, 6],
    table_display[i, 7]
  )
  
  for (j in 1:length(row_data)) {
    text(x_positions[j], y_pos, as.character(row_data[j]),
         cex = 0.95, adj = c(0.5, 0.5))
  }
}

# Add footer note
footer_text <- sprintf("Threshold: |Δz| = %dm | N = %d stations | Analysis period: 2024-2025",
                       threshold, nrow(elev_data))
text(0.5, 0.05, footer_text, cex = 0.85, adj = c(0.5, 0.5), font = 3)

dev.off()
cat("Saved:", output_fig3, "\n")

cat("\n=== ELEVATION BIAS ANALYSIS SUMMARY ===\n\n")
cat("Figure 1: Scatter plots (|Δz| vs RMSE)\n")
cat("- Shows relationship between elevation difference and model error\n")
cat("- Red points: |Δz| > 100m (high elevation difference)\n")
cat("- Blue points: |Δz| ≤ 100m (low elevation difference)\n\n")

cat("Figure 2: Box plots (RMSE distribution)\n")
cat("- Compares RMSE distributions for high vs low elevation difference\n")
cat("- Shows mean values and t-test p-values\n\n")

cat("Figure 3: Summary table\n")
cat("- Statistical comparison across all three models\n")
cat("- Shows correlation and p-values\n\n")

# Print summary statistics to console
print(summary_stats, row.names = FALSE)
