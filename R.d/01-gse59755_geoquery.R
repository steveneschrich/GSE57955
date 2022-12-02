# GEO GSE57955
# Use the `GEOquery` package to quickly retrieve the pre-normalized data.

d <- GEOquery::getGEO("GSE57955", AnnotGPL=TRUE)
d<-d$GSE57955_series_matrix.txt.gz

# The title contains useful information for mapping from their paper (and subsequent paper),
# so we extract the sample information to store.

d$sample_name <- stringr::str_match(d$title, "Penile_Cancer_([0-9A-Za-z]+)_.*$")[,2]

saveRDS(d, file=here::here("data-raw/work/01-gse57955.rds"))
