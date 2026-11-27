# =============================================================================
# 07_seqtab_chimeras.R — ASV table, chimera removal, read tracking
# =============================================================================

library(here)
library(tidyverse)
library(dada2)

samples   <- read_tsv(here("data", "metadata", "samples.tsv"), show_col_types = FALSE)
rds_dir   <- here("output", "rds")
merge_dir <- file.path(rds_dir, "merged")

# --- 1. Load the per-sample results and build the ASV table ----------------------
res <- map(file.path(merge_dir, paste0(samples$sample_alias, ".rds")), readRDS)
names(res) <- samples$sample_alias

seqtab <- makeSequenceTable(map(res, "merged"))
cat("ASV table (samples x ASVs):", dim(seqtab), "\n")

# --- 2. ASV lengths -----------------------------------------------------------------
asv_len <- nchar(getSequences(seqtab))
print(summary(asv_len))
len_plot <- tibble(length = asv_len, reads = colSums(seqtab)) %>%
  ggplot(aes(length, weight = reads)) +
  geom_histogram(binwidth = 2) +
  labs(x = "ASV length (bases)", y = "Reads", title = "Length of merged V9 sequences") +
  theme_bw()
print(len_plot)
ggsave(here("output", "figures", "07_asv_length.png"), len_plot, width = 7, height = 4)

# --- 3. Remove chimeras ----------------------------------------------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)
cat("ASVs after chimera removal:", ncol(seqtab.nochim), "\n")
cat("Share of reads kept:", round(100 * sum(seqtab.nochim) / sum(seqtab), 1), "%\n")

saveRDS(seqtab,        file.path(rds_dir, "seqtab.rds"))
saveRDS(seqtab.nochim, file.path(rds_dir, "seqtab_nochim.rds"))

# --- 4. Track reads through the pipeline ------------------------------------------------
trim <- read_tsv(here("output", "tables", "03_trimming_summary.tsv"), show_col_types = FALSE)
filt <- read_tsv(here("output", "tables", "04_filtering_summary.tsv"), show_col_types = FALSE)

track <- samples %>%
  select(sample_alias, run_accession, extraction, mass_g, replicate) %>%
  left_join(select(trim, run_accession, raw = pairs_in, trimmed = pairs_out), by = "run_accession") %>%
  left_join(select(filt, run_accession, filtered = reads.out), by = "run_accession") %>%
  mutate(denoisedF = map_dbl(res, "denoisedF")[sample_alias],
         denoisedR = map_dbl(res, "denoisedR")[sample_alias],
         merged    = rowSums(seqtab)[sample_alias],
         nonchim   = rowSums(seqtab.nochim)[sample_alias],
         pct_final = round(100 * nonchim / raw, 1))

write_tsv(track, here("output", "tables", "07_read_tracking.tsv"))
print(summary(track$pct_final))
track %>% arrange(pct_final) %>% select(sample_alias, extraction, raw, nonchim, pct_final) %>% head(5) %>% print()