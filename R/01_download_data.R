# =============================================================================
# 01_download_data.R — run table and FASTQ files for PRJNA1359375
# =============================================================================

library(here)
library(tidyverse)

project  <- "PRJNA1359375"
meta_dir <- here("data", "metadata")

# --- 1. Run table from ENA ----------------------------------------------------
fields <- c("run_accession", "sample_accession", "sample_alias", "sample_title",
            "library_name", "read_count", "fastq_ftp", "fastq_md5", "fastq_bytes")

ena_url <- paste0(
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=", project,
  "&result=read_run&fields=", paste(fields, collapse = ","), "&format=tsv"
)

runinfo <- read_tsv(ena_url, show_col_types = FALSE)
write_tsv(runinfo, file.path(meta_dir, "ena_runinfo.tsv"))

# --- Look at it ---------------------------------------------------------------
nrow(runinfo)                                    # expect 40
runinfo %>% select(run_accession, sample_alias, sample_title, library_name) %>% print(n = 40)
sum(runinfo$fastq_bytes %>% str_split(";") %>% unlist() %>% as.numeric()) / 1e9  # total GB

# --- 2. Sample attributes from ENA (sort number, replicate, etc.) --------------
library(xml2)

get_attrs <- function(acc) {
  x <- read_xml(paste0("https://www.ebi.ac.uk/ena/browser/api/xml/", acc))
  tibble(sample_accession = acc,
         tag   = xml_text(xml_find_all(x, "//SAMPLE_ATTRIBUTE/TAG")),
         value = xml_text(xml_find_all(x, "//SAMPLE_ATTRIBUTE/VALUE")))
}

attrs_long <- map_dfr(unique(runinfo$sample_accession), get_attrs)

# Which attribute tags exist, and in how many samples?
count(attrs_long, tag) %>% print(n = Inf)

# One row per sample with the fields we care about
attrs <- attrs_long %>%
  filter(tag %in% c("sort", "rep", "collection_date", "depth")) %>%
  distinct(sample_accession, tag, .keep_all = TRUE) %>%
  pivot_wider(names_from = tag, values_from = value) %>%
  left_join(select(runinfo, run_accession, sample_accession, sample_alias), by = "sample_accession") %>%
  mutate(sort = as.integer(sort)) %>%
  arrange(sort)

write_tsv(attrs, file.path(meta_dir, "biosample_attributes.tsv"))

# Look for the pattern
count(attrs, rep)
count(attrs, collection_date, depth)
attrs %>% select(sort, rep, sample_alias) %>% print(n = 80)


# --- 3. Sample sheet from the authors' ARC files -------------------------------
library(readxl)

# Soil mass and group for each sample (from the authors' study file)
study <- read_excel(file.path(meta_dir, "isa.study.xlsx"), sheet = "study_samples") %>%
  transmute(sample_alias = `Output [Sample Name]`,
            group        = Description,
            mass_g       = as.numeric(`Factor [soil_mass]`))

# SRR run accession for each sample (from the authors' assay file)
assay <- read_excel(file.path(meta_dir, "isa.assay.xlsx"), sheet = "seq_file_id") %>%
  distinct(sample_alias = `Input [Sample Name]`, run_accession = `Comment [SRR accession]`)

# Replicate (s1-s5) from the BioSample attributes we saved earlier
attrs <- read_tsv(file.path(meta_dir, "biosample_attributes.tsv"), show_col_types = FALSE) %>%
  select(sample_alias, replicate = rep)

samples <- study %>%
  left_join(assay, by = "sample_alias") %>%
  left_join(attrs, by = "sample_alias") %>%
  mutate(sort = as.integer(str_extract(sample_alias, "\\d+$")),
         extraction = case_when(
           str_detect(group, "_small$")  ~ "small",         # alternative small extraction
           str_detect(group, "_pooled$") ~ "split_pooled",  # 16/32 g, portions pooled
           mass_g >= 16                  ~ "aliquot_8g",    # single 8 g portion of a 16/32 g sample
           TRUE                          ~ "single"),       # one extraction, 0.25-8 g
         in_paper = extraction %in% c("single", "split_pooled")) %>%
  arrange(sort)

write_tsv(samples, file.path(meta_dir, "samples.tsv"))

# --- Checks ---------------------------------------------------------------------
count(samples, extraction, mass_g, in_paper) %>% print(n = 30)
stopifnot(nrow(samples) == 80,
          sum(samples$in_paper) == 40,
          !anyNA(samples$run_accession),
          !anyNA(samples$replicate))
# The authors' SRR mapping must agree with ENA's
stopifnot(all(samples$run_accession %in% runinfo$run_accession))
message("Sample sheet OK: 80 samples, 40 in the paper.")