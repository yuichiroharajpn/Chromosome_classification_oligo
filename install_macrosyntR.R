#!/usr/bin/env Rscript
# Run after activating the conda environment defined in environment.<OS>.yml.
repository <- "https://cloud.r-project.org"
required_version <- "0.2.19"
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = repository)
}
message("Installing macrosyntR ", required_version)
remotes::install_version(
  "macrosyntR", version = "0.2.19", repos = repository,
  dependencies = NA, upgrade = "never"
)
if (!requireNamespace("macrosyntR", quietly = TRUE)) {
  stop("macrosyntR installation failed; inspect the installation output above.")
}
if (as.character(packageVersion("macrosyntR")) != required_version) {
  stop("Installed macrosyntR is not 0.2.19; inspect installation output.")
}
if (!all(c("sp1_bed", "sp2_bed") %in% names(formals(macrosyntR::load_orthologs)))) {
  stop("Installed macrosyntR lacks the API required by the existing analysis scripts.")
}
message("Installed macrosyntR ", packageVersion("macrosyntR"), " with R ", getRversion())
sessionInfo()
