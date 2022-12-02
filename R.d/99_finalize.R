# 99-finalize

d <- readRDS(here::here("data-raw/work/10-gse57955_annotated.rds"))
saveRDS(d, file = here::here("data/pscc_gse57955_ratio.rds"))

d <- readRDS(here::here("data-raw/work/10-sgse57955_annotated.rds"))
saveRDS(d, file = here::here("data/pscc_gse57955.rds"))
