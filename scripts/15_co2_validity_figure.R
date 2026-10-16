#!/usr/bin/env Rscript
# =============================================================================
# 15_co2_validity_figure.R
#
# Produces a 4x2 panel figure showing raw concentration timeseries for 4
# example chamber closures, illustrating CO2 as a measurement validity
# diagnostic (the Validity section of Results; Fig. 5).
#
# Layout: 4 columns (measurements) x 2 rows (CO2 top, CH4 bottom)
#
# Reads: data/from_hf_stem_study/figure_closures/fig5_co2_tracer.csv (the closures' raw readings,
#        written by 14_figure_closures.R; no access to the raw archive)
# =============================================================================

library(ggplot2)
library(patchwork)
library(dplyr)
source("scripts/style.R")

# =============================================================================
# 1. Define measurement metadata
# =============================================================================

measurements <- data.frame(
  panel_label = c("(a) High signal", "(b) Below detection",
                  "(c) No signal", "(d) Ambiguous"),
  UniqueID    = c("45509_414", "45489_153", "45580_439", "45610_414"),
  date_str    = c("2024-08-05", "2024-07-16", "2024-10-15", "2024-11-14"),
  comp_start  = c("11:28:30", "12:04:13", "13:02:55", "11:31:00"),
  comp_end    = c("11:31:30", "12:07:13", "13:05:55", "11:34:00"),
  valid_co2   = c(TRUE, TRUE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

# Enforce a consistent 180-second window for all panels
MAX_ELAPSED_S <- 180

# Column titles; the verdict is also carried by a soft background tint per column
panel_exprs <- list(
  "**Valid**",
  "**Valid**, CH<sub>4</sub> below MDF",
  "**Invalid**",
  "**Ambiguous**"
)

# Readings of the four closures inside their computation windows (180 s), extracted
# from the raw LGR1 records by 14_figure_closures.R
fig5_raw <- read.csv(file.path("data", "from_hf_stem_study", "figure_closures", "fig5_co2_tracer.csv"),
                     stringsAsFactors = FALSE,
                     colClasses = c(time_s = "numeric", CO2d_ppm = "numeric", CH4d_ppm = "numeric"))

# =============================================================================
# 2. Subset each closure to its computation window
# =============================================================================

read_window <- function(id, date_str, comp_start) {
  t_start <- as.POSIXct(paste(as.Date(date_str), comp_start),
                        format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  w <- fig5_raw[fig5_raw$UniqueID == id, ]
  if (nrow(w) == 0) stop("No readings for ", id, " in data/from_hf_stem_study/figure_closures/fig5_co2_tracer.csv")
  out <- data.frame(elapsed_s = w$time_s - as.numeric(t_start),   # seconds from window start
                    CO2_ppm   = w$CO2d_ppm,
                    CH4_ppb   = w$CH4d_ppm * 1000)
  out[out$elapsed_s >= 0 & out$elapsed_s <= MAX_ELAPSED_S, ]
}

# =============================================================================
# 3. Read all 4 measurements
# =============================================================================

all_data <- list()

for (i in seq_len(nrow(measurements))) {
  cat("Reading measurement:", measurements$panel_label[i], "\n")

  df <- read_window(
    id         = measurements$UniqueID[i],
    date_str   = measurements$date_str[i],
    comp_start = measurements$comp_start[i]
  )

  df$panel       <- measurements$panel_label[i]
  df$valid_co2   <- measurements$valid_co2[i]
  df$panel_order <- i

  all_data[[i]] <- df
}

combined <- bind_rows(all_data)

# Set panel as an ordered factor for consistent facet ordering
combined$panel <- factor(combined$panel,
                         levels = measurements$panel_label,
                         ordered = TRUE)

cat("\nData summary:\n")
combined %>%
  group_by(panel) %>%
  summarise(
    n = n(),
    CO2_range = paste0(round(min(CO2_ppm), 1), " - ", round(max(CO2_ppm), 1)),
    CH4_range = paste0(round(min(CH4_ppb), 1), " - ", round(max(CH4_ppb), 1)),
    .groups = "drop"
  ) %>%
  print()

# =============================================================================
# 4. Compute R-squared and linear fit stats for annotations
# =============================================================================

compute_stats <- function(df, gas_col) {
  fit <- lm(df[[gas_col]] ~ df$elapsed_s)
  r2  <- summary(fit)$r.squared
  slope <- coef(fit)[2]
  return(list(r2 = r2, slope = slope))
}

# Store stats in a named list for easy lookup
stats_list <- list()

for (i in seq_len(nrow(measurements))) {
  sub <- combined[combined$panel == measurements$panel_label[i], ]

  co2_stats <- compute_stats(sub, "CO2_ppm")
  ch4_stats <- compute_stats(sub, "CH4_ppb")

  stats_list[[paste0(i, "_CO2")]] <- co2_stats
  stats_list[[paste0(i, "_CH4")]] <- ch4_stats

  cat(sprintf("Panel %s: CO2 R2=%.4f, CH4 R2=%.4f\n",
              measurements$panel_label[i],
              co2_stats$r2, ch4_stats$r2))
}

# =============================================================================
# 5. Build the 4x2 panel figure with ggplot2 + patchwork
# =============================================================================

# Define validity background colors
bg_by_panel <- c("#ebf3ed", "#ebf3ed", "#f8ebe4", "#f1f1ef")  # valid, valid, invalid, ambiguous

# Helper function to create individual panel
make_panel <- function(panel_idx, gas, y_col, y_lab, show_y_axis,
                       show_x_axis, show_title) {

  panel_name <- measurements$panel_label[panel_idx]
  sub <- combined[combined$panel == panel_name, ]
  bg_color <- bg_by_panel[panel_idx]

  # Get stats
  key <- paste0(panel_idx, "_", gas)
  st <- stats_list[[key]]
  r2_text <- bquote(italic(R)^2 == .(formatC(st$r2, format = "f", digits = 3)))

  # Compute axis ranges for annotation placement
  x_range <- range(sub$elapsed_s)
  y_range <- range(sub[[y_col]])
  y_span  <- diff(y_range)

  # Pad y range slightly so annotation does not clip
  if (y_span < 1e-6) y_span <- abs(y_range[1]) * 0.1 + 1

  # Place R2 label in top-left corner
  ann_x <- x_range[1] + 0.05 * diff(x_range)
  ann_y <- y_range[2] - 0.05 * y_span

  p <- ggplot(sub, aes(x = elapsed_s, y = .data[[y_col]])) +
    # Background rectangle for validity tinting
    annotate("rect",
             xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = bg_color) +
    # Raw data points
    geom_point(size = 0.4, alpha = 0.7, color = INK2) +
    # Linear regression line
    geom_smooth(method = "lm", formula = y ~ x,
                se = FALSE, color = INK, linewidth = 0.6) +
    # R-squared annotation (parsed expression; default thin label border, as published)
    annotate("label", x = ann_x, y = ann_y, fill = scales::alpha("white", 0.8),
             label = deparse(r2_text), parse = TRUE,
             hjust = 0, vjust = 1, size = 2.5,
             color = INK) +
    theme_paper(base_size = 8) +
    theme(plot.margin = margin(2, 4, 2, 2))

  # Y-axis label
  if (show_y_axis) {
    p <- p + labs(y = y_lab)
  } else {
    p <- p + labs(y = NULL) +
      theme(axis.text.y = element_text(),
            axis.ticks.y = element_line())
  }

  # X-axis label
  if (show_x_axis) {
    p <- p + labs(x = "Elapsed time (s)")
  } else {
    p <- p + labs(x = NULL)
  }

  # Panel title (column header; markdown via ggtext)
  if (show_title) {
    p <- p + ggtitle(panel_exprs[[panel_idx]]) +
      theme(plot.title = ggtext::element_markdown(size = 8, hjust = 0))
  }

  return(p)
}

# Build all 8 panels
panels <- list()

for (i in seq_len(nrow(measurements))) {
  show_y <- (i == 1)

  # CO2 panel (top row)
  panels[[paste0("co2_", i)]] <- make_panel(
    panel_idx   = i,
    gas         = "CO2",
    y_col       = "CO2_ppm",
    y_lab       = expression(CO[2] ~ "(ppm)"),
    show_y_axis = show_y,
    show_x_axis = FALSE,
    show_title  = TRUE
  )

  # CH4 panel (bottom row)
  panels[[paste0("ch4_", i)]] <- make_panel(
    panel_idx   = i,
    gas         = "CH4",
    y_col       = "CH4_ppb",
    y_lab       = expression(CH[4] ~ "(ppb)"),
    show_y_axis = show_y,
    show_x_axis = TRUE,
    show_title  = FALSE
  )
}

# Assemble with patchwork: top row (CO2), then bottom row (CH4)
top_row <- panels$co2_1 + panels$co2_2 + panels$co2_3 + panels$co2_4 +
  plot_layout(nrow = 1)

bottom_row <- panels$ch4_1 + panels$ch4_2 + panels$ch4_3 + panels$ch4_4 +
  plot_layout(nrow = 1)

full_fig <- (top_row / bottom_row) + plot_annotation(tag_levels = list(c("a", "b", "c", "d", "", "", "", ""))) &
  theme(plot.tag = element_text(face = "bold", size = 9))

# =============================================================================
# 6. Save figure
# =============================================================================

fig_dir <- Sys.getenv("FIG_DIR", "figures")   # as in 00_setup.R
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

out_png <- file.path(fig_dir, "co2_validity_examples.png")
out_pdf <- file.path(fig_dir, "co2_validity_examples.pdf")

ggsave(out_png, full_fig, width = 7.5, height = 3.9, dpi = 300)
cat("\nSaved PNG:", out_png, "\n")

ggsave(out_pdf, full_fig, width = 7.5, height = 3.9)
cat("Saved PDF:", out_pdf, "\n")

cat("\nDone.\n")
