# 02-raw_gse59755
# We are going to try and read the raw data files (downloaded previously) using
# the agilp Bioconductor package.
#

# Create a working dir for output
unlink(here::here("data-raw/work/02_geo_raw"), recursive=TRUE)
dir.create(here::here("data-raw/work/02_geo_raw"), showWarnings=FALSE)
unlink(here::here("data-raw/work/02_geo_loess"), recursive=TRUE)
dir.create(here::here("data-raw/work/02_geo_loess"), showWarnings=FALSE)
unlink(here::here("data-raw/work/02_geo_id"), recursive = TRUE)
dir.create(here::here("data-raw/work/02_geo_id"), showWarnings=FALSE)

# Load the data using agilp
#
# Per the manual (https://www.bioconductor.org/packages/release/bioc/vignettes/agilp/inst/doc/agilp_manual.pdf)
#
#

# Process first
cli::cli_inform("Running the processing first step")
agilp::AAProcess(
  input = here::here("data-raw/geo_raw"),
  output = here::here("data-raw/work/02_geo_raw/")
)

# Map agilent probe reporters to ensembl annotations. We will use the existing GPL in GEO for this purpose,
# then the agilp package will annotate/remove probes with too high of variability.
cli::cli_inform("Retrieving GPL annotation.")
gpl <- GEOquery::getGEO("GPL6480")
annot <- GEOquery::Table(gpl)
readr::write_tsv(annot, file = here::here("data-raw/work/02-gpl_annotation.txt"))


# Note: ERR of 1 will not discard any multiple probes for the same genes, whereas lower numbers represent
# the coefficient of variation threshold below which probe groups are discarded vs. summarized (using mean
# by default).
#
# NB: This is not strictly needed since we have annotations of the individual probes. Each of the probes
# will be annotated. So for now this is disabled.
#cli::cli_inform("Reannotating GEO data with Ensembl ID's.")
#agilp::IDswop(
#  input = here::here("data-raw/work/02_geo_raw/"),
#  output = here::here("data-raw/work/02_geo_id/"),
#  annotation = here::here("data-raw/work/02-gpl_annotation.txt"),
#  source_ID = "ID",
#  target_ID = "ENSEMBL_ID",
#  ERR = 1
#)

#agilp::filenamex(input=here::here("data-raw/work/02_geo_raw"),output=here::here("data-raw/work/02_geo_raw/"))

# Calculate the baseline normalization factors
cli::cli_inform("Calculating baseline normalization factors.")
agilp::Baseline(NORM="LOG",allfiles=TRUE,r=2,A=2,B=3,input=here::here("data-raw/work/02_geo_raw/"),
         baseout=here::here("data-raw/work/02_geo_raw_baseline.txt"))



# This does a normalization using loess normalization approach, with baselines calculated from last step.
cli::cli_inform("Normalizing data via loess normalization.")
agilp::AALoess(
  input =  here::here("data-raw/work/02_geo_raw/"),
  output=here::here("data-raw/work/02_geo_loess/"),
  baseline = here::here("data-raw/work/02_geo_raw_baseline.txt"),
  LOG=TRUE
)

# Now we should have normalized, individual files. The normalization provides a SSE score which indicates problems
sse <- vroom::vroom(here::here("data-raw/work/02_geo_loess/SSE_2022-12-02.txt"), delim="\t")
sse |>
  dplyr::arrange(dplyr::desc(SSE)) |>
  dplyr::slice_head(n=5)

# This indicates that gRaw_GSM1398516_PE25T.gz_ENSEMBL_ID.txt (SSE = 103393) is possibly an outlier. We can
# at least flag it in the output, and look further when QC is done.

# We load the individual files, annotating as needed,
cli::cli_inform("Loading normalized data.")
arrays <- list.files(here::here("data-raw/work/02_geo_loess"), pattern="*.gz", full.names = TRUE)
array_data <- tibble::tibble(
  channel = ifelse(stringr::str_detect(arrays, "/nr_GSM"), "Cy5","Cy3"),
  reference = ifelse(channel=="Cy5", TRUE, FALSE),
  combined_sample_name = stringr::str_match(arrays, "_(GSM.*).txt.gz")[,2],
  tumor_geo_accession = stringr::str_split_fixed(combined_sample_name, "_", 2)[,1],
  geo_accession = ifelse(reference, paste0(tumor_geo_accession,"_reference"), tumor_geo_accession),
  tumor_sample_name = stringr::str_split_fixed(combined_sample_name, "_", 2)[,2],
  sample_name = ifelse(reference, "", tumor_sample_name),
  exprs = purrr::map2(arrays, geo_accession, function(f, sample_name) {
    vroom::vroom(f, delim="\t", col_names=c("probeid",sample_name), col_types = "cn", skip=1)
  })
)

# Next, we need to extract tumor and reference data separately (because samples are labeled by
# the tumor name).
cli::cli_inform("Creating ExpressionSet.")

expression_data <- array_data |>
  dplyr::pull("exprs") |>
  purrr::map_dfc(~.x[,2]) |>
  dplyr::mutate(probeid = array_data$exprs[[1]]$probeid) |>
  tibble::column_to_rownames("probeid") |>
  as.matrix()


# Subset the gene-level annotation
gpl <- GEOquery::getGEO("GPL6480")
annot <- GEOquery::Table(gpl)

annot <- tibble::enframe(rownames(expression_data), name=NULL, value="probeid") |>
  dplyr::left_join(annot,by=c("probeid"="ID"))



# Now create the Expression object.
x <- Biobase::ExpressionSet(
  assayData = expression_data,
  phenoData = Biobase::AnnotatedDataFrame(
    array_data |>
      dplyr::select(-exprs) |>
      tibble::column_to_rownames("geo_accession")
  ),
  featureData = Biobase::AnnotatedDataFrame(
    tibble::column_to_rownames(annot, "probeid")
  )
)

saveRDS(x, file = here::here("data-raw/work/02-sgse59755.rds"))
