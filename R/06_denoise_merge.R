# =============================================================================
# 06_denoise_merge.R — dada2 denoising and read merging, one sample at a time
# (low-memory design: each sample is saved to disk, then memory is freed)
# =============================================================================

library(here)
library(tidyverse)
library(dada2)

samples   <- read_tsv(here("data", "metadata", "samples.tsv"), show_col_types = FALSE)
filt_dir  <- here("data", "filtered")
rds_dir   <- here("output", "rds")
merge_dir <- file.path(rds_dir, "merged")
dir.create(merge_dir, recursive = TRUE, showWarnings = FALSE)

errF <- readRDS(file.path(rds_dir, "errF.rds"))
errR <- readRDS(file.path(rds_dir, "errR.rds"))

# Paper samples first, "small" samples last
samples <- samples %>%
  mutate(priority = match(extraction, c("single", "split_pooled", "aliquot_8g", "small"))) %>%
  arrange(priority, sort)

for (i in seq_len(nrow(samples))) {
  id       <- samples$sample_alias[i]
  out_file <- file.path(merge_dir, paste0(id, ".rds"))
  if (file.exists(out_file)) next                          # already done
  
  t0 <- Sys.time()
  fF <- file.path(filt_dir, paste0(samples$run_accession[i], "_1.fastq.gz"))
  fR <- file.path(filt_dir, paste0(samples$run_accession[i], "_2.fastq.gz"))
  
  ddF <- dada(fF, err = errF, multithread = TRUE, verbose = FALSE)
  ddR <- dada(fR, err = errR, multithread = TRUE, verbose = FALSE)
  mm  <- mergePairs(ddF, fF, ddR, fR)
  
  saveRDS(list(merged    = mm,
               denoisedF = sum(getUniques(ddF)),
               denoisedR = sum(getUniques(ddR))), out_file)
  
  message(sprintf("[%d/%d] %s  (%s, %s g): %d merged reads, %.1f min",
                  i, nrow(samples), id, samples$extraction[i], samples$mass_g[i],
                  sum(mm$abundance), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  
  rm(ddF, ddR, mm); invisible(gc())                       # free memory
}

message("All samples done.")