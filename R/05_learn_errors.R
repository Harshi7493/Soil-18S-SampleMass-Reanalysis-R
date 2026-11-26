# =============================================================================
# 05_learn_errors.R — learn the sequencing error model (forward and reverse)
# =============================================================================

library(here)
library(tidyverse)
library(dada2)

samples  <- read_tsv(here("data", "metadata", "samples.tsv"), show_col_types = FALSE)
filt_dir <- here("data", "filtered")
rds_dir  <- here("output", "rds")            # large R objects, not tracked by git
dir.create(rds_dir, showWarnings = FALSE)

filtFs <- file.path(filt_dir, paste0(samples$run_accession, "_1.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(samples$run_accession, "_2.fastq.gz"))

set.seed(100)   # makes the run reproducible

# --- 1. Learn errors (uses ~100 million bases, i.e. the first ~9-10 samples) ------
errF <- learnErrors(filtFs, nbases = 2e7, randomize = TRUE,
                    multithread = TRUE, verbose = TRUE)
errR <- learnErrors(filtRs, nbases = 2e7, randomize = TRUE,
                    multithread = TRUE, verbose = TRUE)

saveRDS(errF, file.path(rds_dir, "errF.rds"))
saveRDS(errR, file.path(rds_dir, "errR.rds"))

# --- 2. Plot the error models ---------------------------------------------------------
pF <- plotErrors(errF, nominalQ = TRUE) + ggtitle("Error model, forward reads")
pR <- plotErrors(errR, nominalQ = TRUE) + ggtitle("Error model, reverse reads")
print(pF); print(pR)
ggsave(here("output", "figures", "05_error_model_forward.png"), pF, width = 9, height = 8)
ggsave(here("output", "figures", "05_error_model_reverse.png"), pR, width = 9, height = 8)