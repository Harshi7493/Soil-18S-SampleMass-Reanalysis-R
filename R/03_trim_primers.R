# =============================================================================
# 03_trim_primers.R — remove primers from both ends with cutadapt
# Primers: 1391f (forward) and EukBr (reverse), 18S rRNA V9 region
# =============================================================================

library(here)
library(tidyverse)
library(dada2)
library(ShortRead)
library(Biostrings)

raw_dir  <- here("data", "raw")
trim_dir <- here("data", "trimmed")
log_dir  <- here("output", "cutadapt_logs")
dir.create(trim_dir, showWarnings = FALSE)
dir.create(log_dir,  showWarnings = FALSE)

samples <- read_tsv(here("data", "metadata", "samples.tsv"), show_col_types = FALSE)
fnFs  <- file.path(raw_dir,  paste0(samples$run_accession, "_1.fastq.gz"))
fnRs  <- file.path(raw_dir,  paste0(samples$run_accession, "_2.fastq.gz"))
cutFs <- file.path(trim_dir, basename(fnFs))
cutRs <- file.path(trim_dir, basename(fnRs))

FWD <- "GTACACACCGCCCGTC"            # 1391f
REV <- "TGATCCTTCTGCAGGTTCACCTAC"    # EukBr
FWD.RC <- rc(FWD)
REV.RC <- rc(REV)

# Run cutadapt through the Python launcher
cutadapt <- function(args) system2("py", c("-m", "cutadapt", args), stdout = TRUE, stderr = TRUE)
print(cutadapt("--version"))

# --- 1. Trim every sample ---------------------------------------------------------
for (i in seq_along(fnFs)) {
  if (file.exists(cutFs[i]) && file.exists(cutRs[i])) next    # already done
  message(sprintf("[%d/%d] %s", i, length(fnFs), samples$run_accession[i]))
  out <- cutadapt(c(
    "-a", paste0("^", FWD, "...", REV.RC),   # R1: 1391f at start (required) ... EukBr-RC at end
    "-A", paste0("^", REV, "...", FWD.RC),   # R2: EukBr at start (required) ... 1391f-RC at end
    "--discard-untrimmed",                    # drop pairs without the start primer
    "--minimum-length", "50",                 # drop very short reads
    "-o", cutFs[i], "-p", cutRs[i],
    fnFs[i], fnRs[i]
  ))
  writeLines(out, file.path(log_dir, paste0(samples$run_accession[i], ".log")))
}

# --- 2. How many read pairs survived? ------------------------------------------------
get_num <- function(lines, label) {
  x <- grep(label, lines, value = TRUE)[1]
  as.numeric(gsub(",", "", str_match(x, ":\\s+([0-9,]+)")[, 2]))
}

trim_summary <- samples %>%
  select(run_accession, sample_alias, extraction, mass_g) %>%
  mutate(log = map(file.path(log_dir, paste0(run_accession, ".log")), readLines),
         pairs_in  = map_dbl(log, get_num, "Total read pairs processed"),
         pairs_out = map_dbl(log, get_num, "Pairs written"),
         pct_kept  = round(100 * pairs_out / pairs_in, 1)) %>%
  select(-log)

write_tsv(trim_summary, here("output", "tables", "03_trimming_summary.tsv"))
print(summary(trim_summary$pct_kept))

# --- 3. Check: primers should now be gone ---------------------------------------------
allOrients <- function(primer) {
  dna <- DNAString(primer)
  sapply(c(Forward = dna, Complement = Biostrings::complement(dna),
           Reverse = Biostrings::reverse(dna), RevComp = reverseComplement(dna)), toString)
}
primerHits <- function(primer, fn) sum(vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE) > 0)

ex <- which(samples$sample_alias == "18S-EMB-21")
primer_check <- rbind(
  FWD.ForwardReads = sapply(allOrients(FWD), primerHits, fn = cutFs[ex]),
  FWD.ReverseReads = sapply(allOrients(FWD), primerHits, fn = cutRs[ex]),
  REV.ForwardReads = sapply(allOrients(REV), primerHits, fn = cutFs[ex]),
  REV.ReverseReads = sapply(allOrients(REV), primerHits, fn = cutRs[ex])
)
print(primer_check)
print(summary(width(sread(readFastq(cutFs[ex])))))