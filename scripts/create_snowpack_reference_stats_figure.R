args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "scripts/create_snowpack_reference_stats_figure.R"
root <- normalizePath(file.path(dirname(script_path), ".."))
processed <- file.path(root, "data", "processed")
input_csv <- file.path(processed, "model_comparison_metrics.csv")
output_png <- file.path(processed, "Figure_05b_Model_Comparison_Stats.png")
output_pdf <- file.path(processed, "Figure_05b_Model_Comparison_Stats.pdf")

metrics <- read.csv(input_csv, stringsAsFactors = FALSE)
n_virtual_stations <- nrow(metrics)

model_cols <- list(
  "FSM2-SNP" = c(bias = "bias_fsm", rmse = "rmse_fsm", corr = "corr_fsm"),
  "seNorge-SNP" = c(bias = "bias_senorge", rmse = "rmse_senorge", corr = "corr_senorge")
)

summary_table <- data.frame(
  Model = names(model_cols),
  `Mean RMSE` = sapply(model_cols, function(cols) mean(metrics[[cols["rmse"]]], na.rm = TRUE)),
  `Median RMSE` = sapply(model_cols, function(cols) median(metrics[[cols["rmse"]]], na.rm = TRUE)),
  `Mean Bias` = sapply(model_cols, function(cols) mean(metrics[[cols["bias"]]], na.rm = TRUE)),
  `Median Corr` = sapply(model_cols, function(cols) median(metrics[[cols["corr"]]], na.rm = TRUE)),
  N = n_virtual_stations,
  check.names = FALSE
)

draw_metric_boxplot <- function(metric_name, cols, ylab, title, zero_line = FALSE) {
  boxplot(
    metrics[[cols[1]]],
    metrics[[cols[2]]],
    names = names(cols),
    col = c("#9ecae1", "#fdae6b"),
    border = "#333333",
    outline = TRUE,
    ylab = ylab,
    main = title,
    cex.main = 1.05,
    cex.lab = 0.95,
    cex.axis = 0.9,
    frame.plot = TRUE
  )
  if (zero_line) {
    abline(h = 0, lty = 2, col = "#777777", lwd = 1.2)
  }
  grid(nx = NA, ny = NULL, col = "#e6e6e6", lty = 1)
}

draw_table <- function(table_data) {
  plot.new()
  title("Summary statistics", font.main = 2, cex.main = 1.1)
  display <- table_data
  for (col in names(display)) {
    if (is.numeric(display[[col]]) && col != "N") {
      display[[col]] <- sprintf("%.3f", display[[col]])
    }
  }
  display$N <- as.character(display$N)

  rows <- nrow(display) + 1
  cols <- ncol(display)
  x0 <- 0.02
  y0 <- 0.80
  cell_w <- 0.96 / cols
  cell_h <- 0.18

  for (r in seq_len(rows)) {
    for (c in seq_len(cols)) {
      x_left <- x0 + (c - 1) * cell_w
      y_top <- y0 - (r - 1) * cell_h
      rect(
        x_left,
        y_top - cell_h,
        x_left + cell_w,
        y_top,
        col = if (r == 1) "#f2f2f2" else "white",
        border = "#333333",
        lwd = if (r == 1) 1.2 else 0.9
      )
      label <- if (r == 1) names(display)[c] else display[r - 1, c]
      text(
        x_left + cell_w / 2,
        y_top - cell_h / 2,
        label,
        cex = if (r == 1) 0.82 else 0.78,
        font = if (r == 1) 2 else 1
      )
    }
  }
}

draw_scatter <- function() {
  x <- metrics$rmse_fsm
  y <- metrics$rmse_senorge
  lim <- range(c(x, y), na.rm = TRUE)
  pad <- diff(lim) * 0.04
  lim <- c(max(0, lim[1] - pad), lim[2] + pad)

  plot(
    x,
    y,
    pch = 21,
    bg = adjustcolor("#3182bd", alpha.f = 0.42),
    col = adjustcolor("#08306b", alpha.f = 0.55),
    cex = 0.62,
    xlim = lim,
    ylim = lim,
    xlab = "FSM2-SNP RMSE [m]",
    ylab = "seNorge-SNP RMSE [m]",
    main = "Station-wise RMSE comparison",
    cex.main = 1.05,
    cex.lab = 0.95,
    cex.axis = 0.9
  )
  grid(col = "#e6e6e6", lty = 1)
  abline(0, 1, lty = 2, col = "#777777", lwd = 1.4)
  legend(
    "topleft",
    legend = c("Virtual station", "1:1 line"),
    pch = c(21, NA),
    pt.bg = c(adjustcolor("#3182bd", alpha.f = 0.42), NA),
    lty = c(NA, 2),
    col = c(adjustcolor("#08306b", alpha.f = 0.55), "#777777"),
    bty = "n",
    cex = 0.82
  )
}

draw_figure <- function(output_file, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(output_file, width = 4200, height = 2600, res = 300)
  } else {
    pdf(output_file, width = 14, height = 8.7)
  }

  layout(matrix(c(1, 2, 3, 4, 4, 5), nrow = 2, byrow = TRUE), heights = c(1, 1.18))
  par(
    oma = c(0.5, 0.5, 2.4, 0.5),
    mar = c(4.5, 4.5, 3.0, 1.2),
    family = "Helvetica"
  )

  draw_metric_boxplot(
    "bias",
    c("FSM2-SNP" = "bias_fsm", "seNorge-SNP" = "bias_senorge"),
    "Bias relative to SNP [m]",
    "(a) Bias",
    zero_line = TRUE
  )
  draw_metric_boxplot(
    "rmse",
    c("FSM2-SNP" = "rmse_fsm", "seNorge-SNP" = "rmse_senorge"),
    "RMSE relative to SNP [m]",
    "(b) RMSE"
  )
  draw_metric_boxplot(
    "corr",
    c("FSM2-SNP" = "corr_fsm", "seNorge-SNP" = "corr_senorge"),
    "Temporal correlation [-]",
    "(c) Correlation"
  )

  par(mar = c(1.5, 1.0, 3.0, 1.0))
  draw_table(summary_table)
  par(mar = c(4.5, 4.5, 3.0, 1.2))
  draw_scatter()

  mtext(
    sprintf("Model-to-model statistical agreement relative to SNOWPACK (N = %d virtual stations)", n_virtual_stations),
    outer = TRUE,
    side = 3,
    line = 0.6,
    font = 2,
    cex = 1.35
  )
  dev.off()
}

draw_figure(output_png, "png")
draw_figure(output_pdf, "pdf")

message("Saved: ", output_png)
message("Saved: ", output_pdf)
