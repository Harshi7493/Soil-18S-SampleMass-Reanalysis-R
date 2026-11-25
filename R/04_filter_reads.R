# =============================================================================
# 04_filter_reads.R — quality filtering with dada2::filterAndTrim
# =============================================================================

library(here)
library(tidyverse)
library(dada2)

samples  <- read_tsv(here("data", "metadata", "samples.tsv"), show_col_types = FALSE)
trim_dir <- here("data", "trimmed")
filt_dir <- here("data", "filtered")
dir.create(filt_dir, showWarnings = FALSE)

cutFs  <- file.path(trim_dir, paste0(samples$run_accession, "_1.fastq.gz"))
cutRs  <- file.path(trim_dir, paste0(samples$run_accession, "_2.fastq.gz"))
filtFs <- file.path(filt_dir, basename(cutFs))
filtRs <- file.path(filt_dir, basename(cutRs))

# --- 1. Filter (multithread must be FALSE on Windows) ------------------------------
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
                     maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 50,
                     rm.phix = TRUE, compress = TRUE, multithread = FALSE)

# --- 2. Summary -------------------------------------------------------------------------
filt_summary <- as_tibble(out, rownames = "file") %>%
  mutate(run_accession = str_remove(file, "_1.fastq.gz")) %>%
  left_join(select(samples, run_accession, sample_alias, extraction, mass_g),
            by = "run_accession") %>%
  mutate(pct_kept = round(100 * reads.out / reads.in, 1)) %>%
  select(run_accession, sample_alias, extraction, mass_g, reads.in, reads.out, pct_kept)

write_tsv(filt_summary, here("output", "tables", "04_filtering_summary.tsv"))

print(summary(filt_summary$pct_kept))
print(summary(filt_summary$reads.out))
filt_summary %>% arrange(reads.out) %>% head(5) %>% print()   # samples with fewest reads

# --- 3. Quality after filtering (same two example samples) --------------------------
ex <- which(samples$sample_alias %in% c("18S-EMB-21", "18S-EMB-60"))
qF <- plotQualityProfile(filtFs[ex]) + ggtitle("Forward reads, filtered")
qR <- plotQualityProfile(filtRs[ex]) + ggtitle("Reverse reads, filtered")
print(qF); print(qR)
ggsave(here("output", "figures", "04_quality_forward_filtered.png"), qF, width = 9, height = 4)
ggsave(here("output", "figures", "04_quality_reverse_filtered.png"), qR, width = 9, height = 4)