#!/usr/bin/env Rscript
# =============================================================================
# 15_co2_validity_figure.R
#
# Produces a 4x2 panel figure showing raw concentration timeseries for 4
# example chamber closures, illustrating CO2 as a measurement validity
# diagnostic (manuscript section 4.2).
#
# Layout: 4 columns (measurements) x 2 rows (CO2 top, CH4 bottom)
# =============================================================================

library(ggplot2)
library(patchwork)
library(dplyr)

# =============================================================================
# 1. Define measurement metadata
# =============================================================================

measurements <- data.frame(
  panel_label = c("(a) High signal", "(b) Below detection",
                  "(c) No signal", "(d) Ambiguous"),
  UniqueID    = c("45509_414", "45468_880", "45580_439", "45610_414"),
  date_str    = c("2024-08-05", "2024-06-25", "2024-10-15", "2024-11-14"),
  comp_start  = c("11:28:30", "10:56:01", "13:02:55", "11:31:00"),
  comp_end    = c("11:31:30", "10:59:01", "13:05:55", "11:34:00"),
  valid_co2   = c(TRUE, TRUE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

# Enforce a consistent 180-second window for all panels
MAX_ELAPSED_S <- 180

# For panel titles, use plotmath expressions for subscripts
panel_exprs <- list(
  expression(bold("(a) High signal")),
  expression(bold("(b) Below detection")),
  expression(bold("(c) No signal")),
  expression(bold("(d) Ambiguous"))
)

# Base path for LGR data
lgr_base <- file.path(getwd(), "data", "stem_raw", "LGR", "LGR1")

# =============================================================================
# 2. Function to read LGR file(s) and subset to computation window
# =============================================================================

read_lgr_window <- function(date_str, comp_start, comp_end) {
  dir_path <- file.path(lgr_base, date_str)

  # Build computation window POSIXct boundaries
  date_for_window <- as.Date(date_str)
  t_start <- as.POSIXct(paste(date_for_window, comp_start),
                         format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  t_end   <- as.POSIXct(paste(date_for_window, comp_end),
                         format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

  # Find all f*.txt files (unzipped) in the directory
  txt_files <- list.files(dir_path, pattern = "^micro_.*_f[0-9]+\\.txt$",
                          full.names = TRUE)

  # Also check for zipped f*.txt.zip files
  zip_files <- list.files(dir_path, pattern = "^micro_.*_f[0-9]+\\.txt\\.zip$",
                          full.names = TRUE)

  all_subsets <- list()

  # Read unzipped files
  for (fpath in txt_files) {
    cat("  Reading:", basename(fpath), "\n")
    raw <- read.csv(fpath, skip = 1, stringsAsFactors = FALSE,
                    check.names = FALSE)
    names(raw) <- trimws(names(raw))

    raw$Time_parsed <- as.POSIXct(trimws(raw$Time),
                                  format = "%m/%d/%Y %H:%M:%OS",
                                  tz = "UTC")

    subset_data <- raw[!is.na(raw$Time_parsed) &
                         raw$Time_parsed >= t_start &
                         raw$Time_parsed <= t_end, ]

    if (nrow(subset_data) > 0) {
      all_subsets[[length(all_subsets) + 1]] <- subset_data
    }
  }

  # Read zipped files if we haven't found data yet
  if (length(all_subsets) == 0 && length(zip_files) > 0) {
    for (zpath in zip_files) {
      cat("  Unzipping and reading:", basename(zpath), "\n")
      tmp_dir <- tempdir()
      unzipped <- unzip(zpath, exdir = tmp_dir)

      for (fpath in unzipped) {
        raw <- read.csv(fpath, skip = 1, stringsAsFactors = FALSE,
                        check.names = FALSE)
        names(raw) <- trimws(names(raw))

        raw$Time_parsed <- as.POSIXct(trimws(raw$Time),
                                      format = "%m/%d/%Y %H:%M:%OS",
                                      tz = "UTC")

        subset_data <- raw[!is.na(raw$Time_parsed) &
                             raw$Time_parsed >= t_start &
                             raw$Time_parsed <= t_end, ]

        if (nrow(subset_data) > 0) {
          all_subsets[[length(all_subsets) + 1]] <- subset_data
        }

        # Clean up temp file
        file.remove(fpath)
      }
    }
  }

  if (length(all_subsets) == 0) {
    stop("No data found in computation window for date ", date_str,
         " (", comp_start, " - ", comp_end, ")")
  }

  # Combine all subsets
  combined_raw <- do.call(rbind, all_subsets)

  # Compute elapsed time in seconds from start
  combined_raw$elapsed_s <- as.numeric(
    difftime(combined_raw$Time_parsed, t_start, units = "secs")
  )

  # Extract CO2 and CH4 (dry mole fractions)
  combined_raw$CO2_ppm <- as.numeric(trimws(combined_raw[["[CO2]d_ppm"]]))
  combined_raw$CH4_ppb <- as.numeric(trimws(combined_raw[["[CH4]d_ppm"]])) * 1000

  # Enforce 180-second window cap
  combined_raw <- combined_raw[combined_raw$elapsed_s >= 0 &
                                 combined_raw$elapsed_s <= MAX_ELAPSED_S, ]

  return(combined_raw[, c("elapsed_s", "CO2_ppm", "CH4_ppb")])
}

# =============================================================================
# 3. Read all 4 measurements
# =============================================================================

all_data <- list()

for (i in seq_len(nrow(measurements))) {
  cat("Reading measurement:", measurements$panel_label[i], "\n")

  df <- read_lgr_window(
    date_str   = measurements$date_str[i],
    comp_start = measurements$comp_start[i],
    comp_end   = measurements$comp_end[i]
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
bg_valid   <- "#E8F5E9"  # light green
bg_invalid <- "#FFEBEE"  # light red

# Helper function to create individual panel
make_panel <- function(panel_idx, gas, y_col, y_lab, show_y_axis,
                       show_x_axis, show_title) {

  panel_name <- measurements$panel_label[panel_idx]
  sub <- combined[combined$panel == panel_name, ]
  is_valid <- measurements$valid_co2[panel_idx]
  bg_color <- ifelse(is_valid, bg_valid, bg_invalid)

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
             fill = bg_color, alpha = 0.5) +
    # Raw data points
    geom_point(size = 0.8, alpha = 0.5, color = "grey30") +
    # Linear regression line
    geom_smooth(method = "lm", formula = y ~ x,
                se = FALSE, color = "#1565C0", linewidth = 0.9) +
    # R-squared annotation using plotmath
    annotate("text", x = ann_x, y = ann_y,
             label = deparse(r2_text), parse = TRUE,
             hjust = 0, vjust = 1, size = 3.2,
             color = "grey20") +
    theme_bw(base_size = 11) +
    theme(
      plot.margin = margin(3, 6, 3, 3),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
    )

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

  # Panel title (column header) using plotmath expression
  if (show_title) {
    p <- p + ggtitle(panel_exprs[[panel_idx]]) +
      theme(plot.title = element_text(size = 10, hjust = 0.5))
  }

  return(p)
}

# Build all 8 panels
panels <- list()

for (i in seq_len(nrow(measurements))) {
  show_y <- (i == 1)
  show_title <- TRUE

  # CO2 panel (top row)
  panels[[paste0("co2_", i)]] <- make_panel(
    panel_idx   = i,
    gas         = "CO2",
    y_col       = "CO2_ppm",
    y_lab       = expression(CO[2] ~ "(ppm)"),
    show_y_axis = show_y,
    show_x_axis = FALSE,
    show_title  = show_title
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

full_fig <- top_row / bottom_row

# =============================================================================
# 6. Save figure
# =============================================================================

fig_dir <- file.path(getwd(), "figures")
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

out_png <- file.path(fig_dir, "co2_validity_examples.png")
out_pdf <- file.path(fig_dir, "co2_validity_examples.pdf")

ggsave(out_png, full_fig, width = 12, height = 6, dpi = 300)
cat("\nSaved PNG:", out_png, "\n")

ggsave(out_pdf, full_fig, width = 12, height = 6)
cat("Saved PDF:", out_pdf, "\n")

cat("\nDone.\n")
