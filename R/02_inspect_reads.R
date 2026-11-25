# =============================================================================
# 02_inspect_reads.R — read quality and primer positions (before trimming)
# =============================================================================

library(here)
library(tidyverse)
library(dada2)
library(ShortRead)
library(Biostrings)

raw_dir <- here("data", "raw")
samples <- read_tsv(here("data", "metadata", "samples.tsv"), show_col_types = FALSE)

fnFs <- file.path(raw_dir, paste0(samples$run_accession, "_1.fastq.gz"))
fnRs <- file.path(raw_dir, paste0(samples$run_accession, "_2.fastq.gz"))
stopifnot(all(file.exists(fnFs)), all(file.exists(fnRs)))

# --- 1. Quality profiles for two example samples (1 g and 32 g pooled) ----------
ex <- which(samples$sample_alias %in% c("18S-EMB-21", "18S-EMB-60"))

qF <- plotQualityProfile(fnFs[ex]) + ggtitle("Forward reads")
qR <- plotQualityProfile(fnRs[ex]) + ggtitle("Reverse reads")
print(qF); print(qR)
ggsave(here("output", "figures", "02_quality_forward.png"), qF, width = 9, height = 4)
ggsave(here("output", "figures", "02_quality_reverse.png"), qR, width = 9, height = 4)

# --- 2. Read lengths in one sample ----------------------------------------------
reads_F <- readFastq(fnFs[ex[1]])
length(reads_F)                       # number of reads
summary(width(sread(reads_F)))        # read lengths

# --- 3. Where are the primers? --------------------------------------------------
FWD <- "GTACACACCGCCCGTC"            # 1391f
REV <- "TGATCCTTCTGCAGGTTCACCTAC"    # EukBr

# Every orientation a primer can appear in
allOrients <- function(primer) {
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna),
               Reverse = Biostrings::reverse(dna), RevComp = reverseComplement(dna))
  sapply(orients, toString)
}

# Number of reads in a file that contain a sequence
primerHits <- function(primer, fn) {
  sum(vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE) > 0)
}

primer_table <- rbind(
  FWD.ForwardReads = sapply(allOrients(FWD), primerHits, fn = fnFs[ex[1]]),
  FWD.ReverseReads = sapply(allOrients(FWD), primerHits, fn = fnRs[ex[1]]),
  REV.ForwardReads = sapply(allOrients(REV), primerHits, fn = fnFs[ex[1]]),
  REV.ReverseReads = sapply(allOrients(REV), primerHits, fn = fnRs[ex[1]])
)
primer_table