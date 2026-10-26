###########################################################
#      SPARTA demo: gene-centric SPARTA tests on the simulated genes
###########################################################
# Prerequisite (Step 3.2 on the demo data):
#   Rscript R/SPARTA_null_model.R demo/config_null_model.R
# Then, from the repository root:
#   Rscript demo/run_demo.R
# Runs the SPARTA test (STAAR_sp_2norm_SPA, the function SPARTA_coding.R calls
# for every variant mask) on the 20 simulated genes and writes
#   demo/output/SPARTA_demo_results.csv   (one row per gene, STAAR-style columns)
###########################################################
t0 <- Sys.time()
suppressPackageStartupMessages({library(STAAR); library(Matrix); library(data.table)})

## SPARTA test functions = the function blocks at the top of R/SPARTA_coding.R
src <- parse("R/SPARTA_coding.R")
first_library <- which(sapply(src, function(e) is.call(e) && identical(e[[1]], as.name("library"))))[1]
for (e in src[seq_len(first_library - 1)]) if (is.call(e) && identical(e[[1]], as.name("if"))) eval(e)

## null model from Step 3.2 (obj_nullmodel, resid_ratio, pi1) and the simulated genes
load("demo/output/null_model/obj.STAAR_500K_demo_norm_pseudo_rf100_v5.rda")
load("demo/data/genotype_demo.RData")
ids <- as.character(obj_nullmodel$id_include)
cat(sprintf("Null model: %d samples, pi1 = %.3f, resid_ratio = %.3f\n", length(ids), pi1, resid_ratio))

res <- list()
for (g in names(genotype_list)) {
  G   <- genotype_list[[g]][ids, , drop = FALSE]
  MAF <- Matrix::colMeans(G) / 2
  log <- capture.output(
    pv <- STAAR_sp_2norm_SPA(G, MAF, obj_nullmodel, annotation_list[[g]],
                             rare_maf_cutoff = 0.01, rv_num_cutoff = 2,
                             pi1 = pi1, resid_ratio = resid_ratio))
  res[[g]] <- data.frame(gene = g, num_variant = pv$num_variant, cMAC = pv$cMAC,
                         SPA = any(grepl("use SPA", log)),
                         pv$results_STAAR_S_1_25, pv$results_STAAR_S_1_1,
                         pv$results_STAAR_B_1_25, pv$results_STAAR_B_1_1,
                         pv$results_STAAR_A_1_25, pv$results_STAAR_A_1_1,
                         "ACAT-O" = pv$results_ACAT_O, "STAAR-O" = pv$results_STAAR_O,
                         check.names = FALSE)
}
res <- do.call(rbind, res)
res <- res[order(res[["STAAR-O"]]), ]
dir.create("demo/output", recursive = TRUE, showWarnings = FALSE)
write.csv(res, "demo/output/SPARTA_demo_results.csv", row.names = FALSE)

show <- res[, c("gene", "num_variant", "cMAC", "STAAR-O", "Burden(1,25)", "SKAT(1,25)", "ACAT-V(1,25)", "SPA")]
show[, 4:7] <- lapply(show[, 4:7], signif, 3)
print(show, row.names = FALSE)
cat(sprintf("Results written to demo/output/SPARTA_demo_results.csv  (%.1f s)\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
