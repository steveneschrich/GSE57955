# 10-annotate_geo
# Metadata
# As noted, Daniel extracted the HPV status by sample name. So we can merge this information in here.

d <- readRDS(here::here("data-raw/work/01-gse57955.rds"))
md <- readxl::read_excel(
  here::here(
    "data-raw/hpv_curation/Penile Cancer Dataset_Brazilian_HPV_query_20210218_updated.xlsx"
  ),
  sheet = "Genomic Profiling Paper"
) |>
  dplyr::mutate(sample_name = stringr::str_remove_all(Sample, "^P|T$"))

# NB: The two sets are not the same:
#> setdiff(d$sample_name, md$sample_name)
#[1] "A20" "A22"
#> setdiff(md$sample_name, d$sample_name)
#[1] "A26"
#stopifnot(nrow(md)==ncol(d))
#
# This is ok, we just left-join with the metadata and what's missing is missing
# from the expression side.
pd <- tibble::rownames_to_column(Biobase::pData(d), "rn") |>
  dplyr::left_join(md, by="sample_name") |>
  dplyr::mutate(HPV = dplyr::case_when(
    HPV=="Pos" ~ "+",
    HPV=="Neg" ~ "-"
  )) |>
  tibble::column_to_rownames("rn")

Biobase::pData(d) <- pd

#d$HPV <- (tibble::column_to_rownames(md, "sample_name"))[d$sample_name,"HPV"]
#d$HPV <- c("-","+")[as.numeric(d$HPV=="Pos")+1]

saveRDS(d, here::here("data-raw/work/10-gse57955_annotated.rds"))


# Do the sample-level GSE59755
d <- readRDS(here::here("data-raw/work/02-sgse59755.rds"))
pd <- tibble::rownames_to_column(Biobase::pData(d), "rn") |>
  dplyr::left_join(md |> dplyr::select(-sample_name), by = c("sample_name"="Sample"))
Biobase::pData(d) <- pd
saveRDS(d, here::here("data-raw/work/10-sgse57955_annotated.rds"))
