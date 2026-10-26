## SPARTA demo: settings for R/SPARTA_null_model.R (override its User Input block)
##   Rscript R/SPARTA_null_model.R demo/config_null_model.R
field_all   <- "demo"
pheno_file  <- function(field) "demo/data/pheno_demo.txt"
covar_file  <- "demo/data/covar_demo.txt"
pcs         <- c("sex", "age", "age2", "agesex", "age2sex", paste0("pc", 1:5))
samid_file  <- "demo/data/sample_ids_demo.txt"
grm_file    <- "demo/data/sGRM_demo.RData"
output_path <- "demo/output/null_model/"
