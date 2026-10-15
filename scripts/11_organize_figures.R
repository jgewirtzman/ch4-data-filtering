# ============================================================================
# 11_organize_figures.R — Copy the paper's figures into figures/main_text and
#                         figures/SI under their published numbers
#
# Main figures are complete compositions built in R (17, 15, 19); nothing is
# pasted together here. SI numbers follow the order in which
# the Supplement presents them (Figure S1 is the worked example in Sect. S1).
# Run last (see run_all.R).
# ============================================================================

fig_dir  <- "figures"
main_dir <- file.path(fig_dir, "main_text")
si_dir   <- file.path(fig_dir, "SI")
unlink(main_dir, recursive = TRUE); unlink(si_dir, recursive = TRUE)
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(si_dir,   recursive = TRUE, showWarnings = FALSE)

resolve_fig <- function(f) {
  if (file.exists(file.path(fig_dir, f))) return(file.path(fig_dir, f))
  stop("Figure not found: ", f)
}
copy_set <- function(map, dest) {
  for (nm in names(map)) {
    for (ext in c(".png", ".pdf")) {
      src <- tryCatch(resolve_fig(sub("\\.png$", ext, map[[nm]])), error = function(e) NULL)
      if (!is.null(src)) file.copy(src, file.path(dest, sub("\\.png$", ext, nm)), overwrite = TRUE)
    }
    message("  ", nm, "  <-  ", map[[nm]])
  }
}

message("\n=== Main text figures ===")
copy_set(c(
  "fig1_processing_chain.png"       = "main_fig1_processing_chain.png",
  "fig2_sigma_explained.png"        = "sigma_explained.png",
  "fig3_precision_detection.png"    = "main_fig2_precision_detection.png",
  "fig4_regression_vs_mdf.png"      = "main_fig3_regression_vs_mdf.png",
  "fig5_co2_tracer.png"             = "co2_validity_examples.png",
  "fig6_r2_seasonal.png"            = "main_fig5_r2_seasonal.png",
  "fig7_treatment_budget.png"       = "main_fig6_treatment_budget.png",
  "fig8_transformation.png"         = "main_fig9_transformation.png",
  "fig9_generality.png"             = "main_fig8_generality.png"
), main_dir)

message("\n=== SI figures ===")
copy_set(c(
  "figS01_worked_example_mdf.png"          = "worked_example_mdf.png",
  "figS02_sigma_alternatives.png"          = "sigma_alternatives.png",
  "figS03_logging_interval.png"            = "si_logging_interval.png",
  "figS04_known_truth_simulation.png"      = "known_truth_simulation.png",
  "figS05_mdf_vs_duration.png"             = "mdf_vs_duration.png",
  "figS06_transform_demo.png"              = "si_transform_demo.png"
), si_dir)

message("\nMain: ", length(list.files(main_dir, "\\.png$")), " figures; SI: ",
        length(list.files(si_dir, "\\.png$")), " figures")
