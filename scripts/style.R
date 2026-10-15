## ============================================================================
## style.R — One visual language for every figure in the paper and SI
##
## Sourced by 00_setup.R. Colour follows MEANING and every role keeps its colour
## in every figure:
##
##   Detection / filter criterion (what a closure is judged by, and what is
##   removed by it):  campaign-sigma MDF = orange, R^2 = violet,
##   regression (slope p-value or interval) = aqua, datasheet MDF = grey,
##   no filter / keep all = ink.
##   Sign of a detected flux: emission = red, uptake = blue, below MDF = grey.
##   Analyzer: A = green, B = magenta.
##   Temperature: one warm sequential ramp, pale (cold) to dark (warm).
##
## Blues, greens, greys and black, with one red. Every pair shown together was
## checked for separation under normal, protan, deutan and tritan vision.
## ============================================================================

INK    <- "#1f1f1f"   # primary text and the "keep all" reference (near-black)
INK2   <- "#606060"   # secondary text
MUTED  <- "#9a9a9a"   # reference lines, annotations
GRID   <- "#ececec"   # hairline grid
AXIS   <- "#d2d2d2"   # baselines
GREY   <- "#c8c8c6"   # below-MDF / not detected / kept points
GREY_D <- "#a3a3a0"   # datasheet / set to zero
TINT_OK  <- "#e9f2ee" # pale teal: region where only the regression detects a flux (Fig 4)
TINT_NEU <- "#f3f3f2" # neutral tint: dormant season / shaded region

# A small set of muted hues (slate blue, terracotta, sage teal) plus greys and
# black. Each figure picks what reads best for its own content; the names below
# are just convenient handles.
COL <- list(
  mdf      = "#3d6f8e",   # slate blue
  r2       = "#c2683f",   # terracotta
  reg      = "#5b9a85",   # sage teal
  datasheet = GREY_D,
  keep     = INK,
  zero     = GREY_D,
  valid    = "#7fa98a",   # sage: Fig 5 verdict chip
  invalid  = "#b87a62",   # clay: Fig 5 verdict chip
  emission = "#b87a62",
  uptake   = "#2b6a9a",
  below    = GREY,
  A        = INK,         # noisier analyzer: black (solid lines / filled marks)
  B        = "#8c8c88"    # precise analyzer: grey (dashed lines / open marks)
)

PAL_ANALYZER <- c("Noisier analyzer" = COL$A, "Precise analyzer" = COL$B)
# Ordered ramps (light -> dark) for variants of one criterion, and for magnitude
RAMP_R2   <- c("#eac0aa", "#d8946f", "#c2683f", "#8c4526")
RAMP_X3T  <- c("#b9b9b6", "#8c8c88", "#5f5f5c")
RAMP_BLUE <- c("#9cc0dc", "#4f8dbd", "#173f5f")   # e.g. near-zero -> large flux
R2_LIGHT <- RAMP_R2[1]; R2_DARK <- RAMP_R2[4]

PAPER_FONT <- "Helvetica"

# Base theme: white surface, horizontal hairline grid, no panel box, bold
# left-aligned panel titles, strips as plain bold labels.
# Optional standard ggplot looks for comparison: PAPER_THEME = "bw" or "classic"
theme_std <- function(kind, base_size, base_family) {
  base <- if (kind == "bw") theme_bw(base_size = base_size, base_family = base_family)
          else theme_classic(base_size = base_size, base_family = base_family)
  base + theme(legend.position = "top", legend.justification = "left", legend.location = "plot",
               legend.margin = margin(0, 0, 0, 0), legend.key.size = unit(10, "pt"),
               plot.title = element_text(size = rel(1.0)), plot.title.position = "panel",
               plot.tag = element_text(face = "bold", size = rel(1.15)),
               strip.background = if (kind == "bw") element_rect(fill = "grey92", colour = "grey20") else element_blank(),
               panel.grid.minor = element_blank())
}
theme_paper <- function(base_size = 9, grid = c("y", "xy", "none"),
                        base_family = PAPER_FONT) {
  grid <- match.arg(grid)
  kind <- Sys.getenv("PAPER_THEME", "classic")   # paper default: theme_classic(); "minimal" restores the earlier look
  if (kind %in% c("bw", "classic")) return(theme_std(kind, base_size * 1.35, base_family))
  t <- theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      text             = element_text(colour = INK),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = GRID, linewidth = 0.35),
      axis.line        = element_blank(),
      axis.ticks       = element_blank(),
      axis.ticks.length = unit(0, "pt"),
      axis.text        = element_text(colour = INK2, size = rel(0.9)),
      axis.title       = element_text(colour = INK2, size = rel(0.95)),
      strip.text       = element_text(face = "bold", hjust = 0, size = rel(0.95), colour = INK,
                                      margin = margin(2, 0, 3, 0)),
      legend.title     = element_text(size = rel(0.9), colour = INK2),
      legend.text      = element_text(size = rel(0.9), colour = INK),
      legend.key.size  = unit(10, "pt"),
      legend.position  = "top",
      legend.justification = "left",
      legend.location  = "plot",
      legend.margin    = margin(0, 0, 0, 0),
      plot.title       = element_text(face = "bold", hjust = 0, size = rel(1.0), colour = INK,
                                      margin = margin(0, 0, 4, 0)),
      plot.subtitle    = element_text(colour = INK2, hjust = 0, size = rel(0.9),
                                      margin = margin(0, 0, 4, 0)),
      plot.caption     = element_text(colour = MUTED, hjust = 0, size = rel(0.8)),
      plot.tag         = element_text(face = "bold", size = rel(1.15)),
      plot.title.position = "plot",
      panel.spacing    = unit(10, "pt"),
      plot.margin      = margin(4, 6, 4, 4)
    )
  if (grid == "y")    t <- t + theme(panel.grid.major.x = element_blank())
  if (grid == "none") t <- t + theme(panel.grid.major = element_blank())
  t
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
