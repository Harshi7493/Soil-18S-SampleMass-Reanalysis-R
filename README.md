# Soil-18S-SampleMass-Reanalysis-R

Does the body-size effect of soil sample mass hold when the split-extracted samples are removed?
A reanalysis in R of Wang et al. (2026), *Applied Soil Ecology* 227, 107400.

## Background

Wang et al. extracted DNA from 0.25 to 32 g of soil (8 mass levels × 5 replicates = 40 samples)
and sequenced the 18S rRNA V9 region. Detected diversity rose with soil mass for larger soil
animals (nematodes, arachnids, springtails, insects) but not for fungi and protists.

The DNA kit takes at most 10 g, so the 16 g and 32 g samples were split into 8 g portions,
extracted separately and pooled. The authors note this may add technical variation.

## Question

Using only the single-extraction samples (0.25–8 g, n = 30):

1. Do the group-level slopes of diversity against soil mass keep their sign and significance?
2. Do larger-bodied groups still respond more strongly than small-bodied groups?
3. Does replicate dispersion change sharply at 16 g and 32 g?

## Data

- Raw reads: NCBI SRA BioProject PRJNA1359375 (Illumina MiSeq 2 × 300, 18S V9, primers 1391f / EukBr)
- Taxonomy reference: PR2 v5.0.0

## Reference

Wang L, Leonar AL, Cesarz S, Eisenhauer N, Chatzinotas A, Jurburg SD (2026). Sample size shapes
metabarcoding-driven biodiversity assessments across body sizes in soil.
*Applied Soil Ecology* 227: 107400. https://doi.org/10.1016/j.apsoil.2026.107400