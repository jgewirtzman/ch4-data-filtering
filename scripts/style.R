## ============================================================================
## style.R — One visual language for every figure in the paper and SI
##
## Sourced by 00_setup.R (and by 15_co2_validity_figure.R directly). Colour
## follows meaning, and every role keeps its colour in every figure:
##
##   Detection / filter criterion: campaign-sigma MDF = slate blue, R^2 =
##   terracotta, regression (slope p-value or interval) = sage teal,
##   datasheet MDF and "set to zero" = mid grey, no filter / keep all = near-black.
##   Analyzer: noisier = near-black (solid lines, filled marks), precise = grey
##   (dashed lines, open marks).
##   Magnitude (e.g. near-zero -> large flux): one blue ramp, pale to dark.
##
## Three muted hues plus greys and near-black, on theme_classic().
## ============================================================================

INK    <- "#1f1f1f"   # primary text and the "keep all" reference (near-black)
INK2   <- "#606060"   # secondary text
MUTED  <- "#9a9a9a"   # reference lines, annotations
AXIS   <- "#d2d2d2"   # baselines
GREY   <- "#c8c8c6"   # below-MDF / not detected / kept points
GREY_D <- "#a3a3a0"   # datasheet / set to zero
TINT_OK  <- "#e9f2ee" # pale teal: region where only the regression detects a flux (Fig 4)
TINT_NEU <- "#f3f3f2" # neutral tint: dormant season / shaded region

# Role colours
COL <- list(
  mdf      = "#3d6f8e",   # slate blue
  r2       = "#c2683f",   # terracotta
  reg      = "#5b9a85",   # sage teal
  datasheet = GREY_D,
  keep     = INK,
  zero     = GREY_D,
  A        = INK,         # noisier analyzer: black (solid lines / filled marks)
  B        = "#8c8c88"    # precise analyzer: grey (dashed lines / open marks)
)

PAL_ANALYZER <- c("Noisier analyzer" = COL$A, "Precise analyzer" = COL$B)
# Ordered ramp (light -> dark) for magnitude
RAMP_BLUE <- c("#9cc0dc", "#4f8dbd", "#173f5f")   # e.g. near-zero -> large flux

PAPER_FONT <- "Helvetica"

# Base theme: theme_classic() scaled up by 1.35, legend above the panels and
# left-aligned, panel titles at the panel, bold tags, no strip box, no minor grid.
theme_paper <- function(base_size = 9, base_family = PAPER_FONT) {
  theme_classic(base_size = base_size * 1.35, base_family = base_family) +
    theme(legend.position = "top", legend.justification = "left", legend.location = "plot",
          legend.margin = margin(0, 0, 0, 0), legend.key.size = unit(10, "pt"),
          plot.title = element_text(size = rel(1.0)), plot.title.position = "panel",
          plot.tag = element_text(face = "bold", size = rel(1.15)),
          strip.background = element_blank(),
          panel.grid.minor = element_blank())
}
theme_set(theme_paper())
update_geom_defaults("point", list(colour = INK, size = 1))
update_geom_defaults("line",  list(colour = INK, linewidth = 0.5))
update_geom_defaults("text",  list(colour = INK2, family = PAPER_FONT, size = 2.6))
update_geom_defaults("label", list(colour = INK2, family = PAPER_FONT, size = 2.6))
update_geom_defaults("hline", list(colour = MUTED, linewidth = 0.3))
update_geom_defaults("vline", list(colour = MUTED, linewidth = 0.35))
update_geom_defaults("abline", list(colour = MUTED, linewidth = 0.35))

# Shared axis label for CH4 flux
LAB_FLUX <- expression(CH[4]~flux~(nmol~m^{-2}~s^{-1}))

# Figure saving: every figure of the paper goes through here. Titles and
# subtitles of single-plot figures are dropped (the caption carries them);
# PDFs go through cairo so that sigma and R-squared glyphs render.
ggsave <- function(filename, plot = ggplot2::last_plot(), device = NULL,
                   keep_titles = FALSE, bg = "white", ...) {
  if (!keep_titles && inherits(plot, "ggplot") && !inherits(plot, "patchwork"))
    plot <- plot + labs(title = NULL, subtitle = NULL)
  if (is.null(device) && grepl("\\.pdf$", filename, ignore.case = TRUE))
    device <- grDevices::cairo_pdf
  if (is.null(device) && grepl("\\.png$", filename, ignore.case = TRUE))
    device <- ragg::agg_png
  ggplot2::ggsave(filename, plot = plot, device = device, bg = bg, ...)
}

# Legend placement: above the panels, left-aligned (default), or inside a panel
# where it has empty space. Never directly under the x-axis title.
legend_top <- function() theme(legend.position = "top", legend.justification = "left", legend.location = "plot")
legend_in <- function(x, y) theme(legend.position = "inside", legend.position.inside = c(x, y),
                                  legend.justification = c(x, y),
                                  legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
                                  legend.box.background = element_blank())
