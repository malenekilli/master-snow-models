args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "scripts/create_observation_reference_stats_figure.R"
root <- normalizePath(file.path(dirname(script_path), ".."))
processed <- file.path(root, "data", "processed")
input_csv <- file.path(processed, "observation_validation_three_models.csv")
output_png <- file.path(processed, "Figure_12b_Observation_Reference_Model_Comparison_Stats.png")
output_pdf <- file.path(processed, "Figure_12b_Observation_Reference_Model_Comparison_Stats.pdf")

metrics <- read.csv(input_csv, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
n_stations <- nrow(metrics)

model_cols <- list(
  "FSM2" = c(bias = "FSM2_Bias", rmse = "FSM2_RMSE", corr = "FSM2_Corr", n = "FSM2_N"),
  "seNorge" = c(bias = "seNorge_Bias", rmse = "seNorge_RMSE", corr = "seNorge_Corr", n = "seNorge_N"),
  "SNOWPACK" = c(bias = "SNOWPACK_Bias", rmse = "SNOWPACK_RMSE", corr = "SNOWPACK_Corr", n = "SNOWPACK_N")
)

model_colors <- c("FSM2" = "#9ecae1", "seNorge" = "#fdae6b", "SNOWPACK" = "#a1d99b")
model_point_colors <- c("FSM2" = "#3182bd", "seNorge" = "#e6550d", "SNOWPACK" = "#31a354")

summary_table <- data.frame(
  Model = names(model_cols),
  `Mean RMSE` = sapply(model_cols, function(cols) mean(metrics[[cols["rmse"]]], na.rm = TRUE)),
  `Median RMSE` = sapply(model_cols, function(cols) median(metrics[[cols["rmse"]]], na.rm = TRUE)),
  `Mean Bias` = sapply(model_cols, function(cols) mean(metrics[[cols["bias"]]], na.rm = TRUE)),
  `Median Corr` = sapply(model_cols, function(cols) median(metrics[[cols["corr"]]], na.rm = TRUE)),
  `Mean days` = sapply(model_cols, function(cols) mean(metrics[[cols["n"]]], na.rm = TRUE)),
  N = n_stations,
  check.names = FALSE
)

draw_metric_boxplot <- function(cols, ylab, title, zero_line = FALSE) {
  model_names <- names(cols)
  boxplot(
    lapply(cols, function(col) metrics[[col]]),
    names = model_names,
    col = model_colors[model_names],
    border = "#333333",
    outline = TRUE,
    ylab = ylab,
    main = title,
    cex.main = 1.05,
    cex.lab = 0.95,
    cex.axis = 0.88,
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
  y0 <- 0.83
  cell_w <- 0.96 / cols
  cell_h <- 0.145

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
        cex = if (r == 1) 0.72 else 0.70,
        font = if (r == 1) 2 else 1
      )
    }
  }
}

draw_station_rmse <- function() {
  rmse_matrix <- cbind(
    FSM2 = metrics$FSM2_RMSE,
    seNorge = metrics$seNorge_RMSE,
    SNOWPACK = metrics$SNOWPACK_RMSE
  )
  station_order <- order(rowMeans(rmse_matrix, na.rm = TRUE))
  station_labels <- metrics$Station[station_order]
  x <- seq_along(station_order)
  ylim <- c(0, max(rmse_matrix, na.rm = TRUE) * 1.08)

  plot(
    x,
    rmse_matrix[station_order, "FSM2"],
    type = "n",
    xaxt = "n",
    ylim = ylim,
    xlab = "Observation station",
    ylab = "RMSE relative to observations [m]",
    main = "Station-wise RMSE against observations",
    cex.main = 1.05,
    cex.lab = 0.95,
    cex.axis = 0.9
  )
  grid(col = "#e6e6e6", lty = 1)
  axis(1, at = x, labels = seq_along(x), cex.axis = 0.85)

  offsets <- c(FSM2 = -0.16, seNorge = 0, SNOWPACK = 0.16)
  for (model in colnames(rmse_matrix)) {
    points(
      x + offsets[model],
      rmse_matrix[station_order, model],
      pch = 21,
      bg = adjustcolor(model_point_colors[model], alpha.f = 0.70),
      col = "#333333",
      cex = 1.05
    )
  }
  legend(
    "topleft",
    legend = colnames(rmse_matrix),
    pch = 21,
    pt.bg = adjustcolor(model_point_colors[colnames(rmse_matrix)], alpha.f = 0.70),
    col = "#333333",
    bty = "n",
    cex = 0.82
  )
  text(
    x,
    par("usr")[3] - 0.06 * diff(par("usr")[3:4]),
    labels = substr(station_labels, 1, 10),
    srt = 45,
    adj = 1,
    xpd = TRUE,
    cex = 0.55
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
    c("FSM2" = "FSM2_Bias", "seNorge" = "seNorge_Bias", "SNOWPACK" = "SNOWPACK_Bias"),
    "Bias relative to observations [m]",
    "(a) Bias",
    zero_line = TRUE
  )
  draw_metric_boxplot(
    c("FSM2" = "FSM2_RMSE", "seNorge" = "seNorge_RMSE", "SNOWPACK" = "SNOWPACK_RMSE"),
    "RMSE relative to observations [m]",
    "(b) RMSE"
  )
  draw_metric_boxplot(
    c("FSM2" = "FSM2_Corr", "seNorge" = "seNorge_Corr", "SNOWPACK" = "SNOWPACK_Corr"),
    "Temporal correlation [-]",
    "(c) Correlation"
  )

  par(mar = c(1.5, 1.0, 3.0, 1.0))
  draw_table(summary_table)
  par(mar = c(7.0, 4.5, 3.0, 1.2))
  draw_station_rmse()

  mtext(
    sprintf("Model performance relative to observations (N = %d stations)", n_stations),
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
