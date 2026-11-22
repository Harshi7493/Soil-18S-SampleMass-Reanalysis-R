# =============================================================================
# 00_setup.R — install and check the packages used in this project
# Run once. Safe to re-run: it only installs what is missing.
# =============================================================================

cran_pkgs <- c("here", "tidyverse", "vegan", "metafor", "iNEXT", "broom")
bioc_pkgs <- c("dada2", "ShortRead", "Biostrings", "phyloseq", "microbiome")

# Install missing CRAN packages
missing_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing_cran) > 0) install.packages(missing_cran)

# Install missing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
missing_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

# Check that everything loads, and record the versions
all_pkgs <- c(cran_pkgs, bioc_pkgs)
versions <- data.frame(
  package = all_pkgs,
  version = sapply(all_pkgs, function(p) as.character(packageVersion(p)))
)
print(versions, row.names = FALSE)

writeLines(
  c(R.version.string, "", paste(versions$package, versions$version, sep = "\t")),
  "docs/package_versions.txt"
)
message("Setup complete. Versions saved to docs/package_versions.txt")