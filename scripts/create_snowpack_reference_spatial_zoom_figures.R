args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "scripts/create_snowpack_reference_spatial_zoom_figures.R"
root <- normalizePath(file.path(dirname(script_path), ".."))
processed <- file.path(root, "data", "processed")
metrics <- read.csv(file.path(processed, "model_comparison_metrics.csv"), stringsAsFactors = FALSE)

jotun_bounds <- c(
  west = min(metrics$lon, na.rm = TRUE) - 0.10,
  east = max(metrics$lon, na.rm = TRUE) + 0.10,
  south = min(metrics$lat, na.rm = TRUE) - 0.08,
  north = max(metrics$lat, na.rm = TRUE) + 0.05
)

models <- list(
  list(label = "FSM2-SNP", rmse = "rmse_fsm", bias = "bias_fsm", corr = "corr_fsm"),
  list(label = "seNorge-SNP", rmse = "rmse_senorge", bias = "bias_senorge", corr = "corr_senorge")
)

make_palette <- function(type, n = 256) {
  if (type == "rmse") {
    colorRampPalette(c("#f7fbff", "#deebf7", "#9ecae1", "#3182bd", "#08519c"))(n)
  } else if (type == "bias") {
    colorRampPalette(c("#2166ac", "#67a9cf", "#f7f7f7", "#ef8a62", "#b2182b"))(n)
  } else {
    colorRampPalette(c("#fff7bc", "#fec44f", "#31a354", "#006837"))(n)
  }
}

value_to_col <- function(values, limits, palette) {
  idx <- round((values - limits[1]) / diff(limits) * (length(palette) - 1)) + 1
  idx <- pmin(pmax(idx, 1), length(palette))
  palette[idx]
}

draw_colorbar <- function(limits, palette, label) {
  usr <- par("usr")
  x0 <- usr[2] - 0.070 * diff(usr[1:2])
  x1 <- usr[2] - 0.045 * diff(usr[1:2])
  y_bottom <- usr[3] + 0.17 * diff(usr[3:4])
  y_top <- usr[4] - 0.17 * diff(usr[3:4])
  y <- seq(y_bottom, y_top, length.out = length(palette) + 1)
  for (i in seq_along(palette)) {
    rect(x0, y[i], x1, y[i + 1], col = palette[i], border = NA)
  }
  rect(x0, y_bottom, x1, y_top, border = "#333333")
  ticks <- pretty(limits, n = 5)
  ticks <- ticks[ticks >= limits[1] & ticks <= limits[2]]
  tick_y <- y_bottom + (ticks - limits[1]) / diff(limits) * (y_top - y_bottom)
  segments(x0, tick_y, x0 - 0.010 * diff(usr[1:2]), tick_y)
  text(x0 - 0.018 * diff(usr[1:2]), tick_y, labels = ticks, cex = 0.66, adj = 1)
  text(x1 + 0.035 * diff(usr[1:2]), mean(c(y_bottom, y_top)), labels = label, cex = 0.72, srt = 90)
}

draw_panel <- function(data, value_col, title, limits, palette) {
  plot(
    data$lon,
    data$lat,
    type = "n",
    xlim = jotun_bounds[c("west", "east")],
    ylim = jotun_bounds[c("south", "north")],
    xlab = "",
    ylab = "Latitude [deg N]",
    main = title,
    cex.main = 1.05,
    cex.lab = 0.9,
    cex.axis = 0.82,
    asp = 1 / cos(mean(jotun_bounds[c("south", "north")]) * pi / 180)
  )
  grid(col = "#e5e5e5", lty = 1)
  cols <- value_to_col(data[[value_col]], limits, palette)
  points(
    data$lon,
    data$lat,
    pch = 21,
    bg = adjustcolor(cols, alpha.f = 0.78),
    col = adjustcolor("#333333", alpha.f = 0.18),
    cex = 0.72,
    lwd = 0.25
  )
}

make_figure <- function(metric, cols, title, legend_label, output_base, limits = NULL) {
  palette <- make_palette(metric)
  values <- unlist(lapply(cols, function(col) metrics[[col]]))
  values <- values[is.finite(values)]
  if (is.null(limits)) {
    if (metric == "bias") {
      vmax <- max(abs(values), na.rm = TRUE)
      limits <- c(-vmax, vmax)
    } else {
      limits <- range(values, na.rm = TRUE)
    }
  }

  png(file.path(processed, paste0(output_base, ".png")), width = 3600, height = 1800, res = 300)
  layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))
  par(oma = c(0.4, 0.3, 2.2, 0.3), mar = c(3.6, 4.2, 2.8, 1.0), family = "Helvetica")
  draw_panel(metrics, cols[1], names(cols)[1], limits, palette)
  draw_colorbar(limits, palette, legend_label)
  draw_panel(metrics, cols[2], names(cols)[2], limits, palette)
  draw_colorbar(limits, palette, legend_label)
  mtext(title, side = 3, outer = TRUE, line = 0.6, font = 2, cex = 1.25)
  dev.off()

  pdf(file.path(processed, paste0(output_base, ".pdf")), width = 12, height = 6)
  layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))
  par(oma = c(0.4, 0.3, 2.2, 0.3), mar = c(3.6, 4.2, 2.8, 1.0), family = "Helvetica")
  draw_panel(metrics, cols[1], names(cols)[1], limits, palette)
  draw_colorbar(limits, palette, legend_label)
  draw_panel(metrics, cols[2], names(cols)[2], limits, palette)
  draw_colorbar(limits, palette, legend_label)
  mtext(title, side = 3, outer = TRUE, line = 0.6, font = 2, cex = 1.25)
  dev.off()

  message("Saved: ", file.path(processed, paste0(output_base, ".png")))
  message("Saved: ", file.path(processed, paste0(output_base, ".pdf")))
}

make_figure(
  metric = "rmse",
  cols = c("FSM2-SNP" = "rmse_fsm", "seNorge-SNP" = "rmse_senorge"),
  title = sprintf("Spatial RMSE relative to SNOWPACK (N = %d virtual stations)", nrow(metrics)),
  legend_label = "RMSE relative to SNOWPACK [m]",
  output_base = "Figure_05c_Spatial_RMSE_Jotun_Zoomed",
  limits = c(0, max(metrics$rmse_fsm, metrics$rmse_senorge, na.rm = TRUE))
)

make_figure(
  metric = "bias",
  cols = c("FSM2-SNP" = "bias_fsm", "seNorge-SNP" = "bias_senorge"),
  title = sprintf("Spatial bias relative to SNOWPACK (N = %d virtual stations)", nrow(metrics)),
  legend_label = "Bias relative to SNOWPACK [m]",
  output_base = "Figure_06b_Spatial_Bias_Jotun_Zoomed"
)

make_figure(
  metric = "corr",
  cols = c("FSM2-SNP" = "corr_fsm", "seNorge-SNP" = "corr_senorge"),
  title = sprintf("Spatial temporal correlation relative to SNOWPACK (N = %d virtual stations)", nrow(metrics)),
  legend_label = "Temporal correlation [-]",
  output_base = "Figure_07b_Spatial_Correlation_Jotun_Zoomed",
  limits = c(0, 1)
)
