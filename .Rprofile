# Project-only library holding the goFlux fork (v0.5.0.9002), so that other projects on
# this machine keep their own goFlux. Install with:
#   remotes::install_github("jgewirtzman/goFlux@v0.5.0.9002",
#                           lib = "~/Library/R/ch4-data-filtering-lib", dependencies = FALSE)
local({
  lib <- path.expand("~/Library/R/ch4-data-filtering-lib")
  if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
})
