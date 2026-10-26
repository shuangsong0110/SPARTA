###########################################################
#      SPARTA demo: simulate a small dataset
###########################################################
# Generates the files in demo/data/. Everything is simulated (sample ids are
# demo00001 ... demo05000). Run from the repository root:
#   Rscript demo/simulate_demo_data.R
#
# Design
#   N = 5,000 samples; covariates age, sex, pc1-pc5
#   250 sibling pairs (kinship 0.5)                -> sGRM_demo.RData (sparse GRM `sGRM`)
#   20 genes x 100 rare variants (MAF 1/2N - 1%)   -> genotype_demo.RData
#     genes 1-2 causal (gene01: same-direction effects, gene02: mixed signs), genes 3-20 null;
#     two simulated annotation scores (PHRED-like) per variant
#   s_pheno: predicted phenotype for ALL samples (correlation ~0.85 with the trait)
#   pheno  : observed in a random 30% of samples, NA otherwise  -> pheno_demo.txt (Step 2 format)
###########################################################
suppressPackageStartupMessages({library(Matrix); library(data.table)})
set.seed(2026)

N <- 5000; n_gene <- 20; m <- 100
pi_obs   <- 0.30     # fraction of samples with an observed phenotype
cor_pred <- 0.85     # correlation between predicted and true phenotype
n_pair   <- 250      # sibling pairs
beta1 <- 0.45; beta2 <- 0.5; frac_causal <- 0.4

ids <- sprintf("demo%05d", seq_len(N))
dir.create("demo/data", recursive = TRUE, showWarnings = FALSE)

## ---- covariates ----
covar <- data.table(IID = ids, age = round(runif(N, 40, 70)), sex = rbinom(N, 1, 0.5))
pcs <- matrix(round(rnorm(N * 5), 3), N, 5, dimnames = list(NULL, paste0("pc", 1:5)))
covar <- cbind(covar, as.data.table(pcs))

## ---- sparse GRM: identity + sibling pairs ----
idx <- sample(N, 2 * n_pair); i1 <- idx[1:n_pair]; i2 <- idx[-(1:n_pair)]
sGRM <- sparseMatrix(i = c(1:N, i1, i2), j = c(1:N, i2, i1),
                     x = c(rep(1, N), rep(0.5, 2 * n_pair)),
                     dims = c(N, N), dimnames = list(ids, ids))
save(sGRM, file = "demo/data/sGRM_demo.RData")

## ---- genotypes and annotations ----
genotype_list <- annotation_list <- list()
y_gene <- numeric(N)
for (g in seq_len(n_gene)) {
  gname <- sprintf("gene%02d", g)
  maf <- 10^runif(m, log10(1 / (2 * N)), log10(0.01))
  G <- Matrix(rbinom(N * m, 2, rep(maf, each = N)), nrow = N, ncol = m, sparse = TRUE,
              dimnames = list(ids, sprintf("%s_v%03d", gname, seq_len(m))))
  genotype_list[[gname]]   <- as(G, "dgCMatrix")
  annotation_list[[gname]] <- data.frame(Annotation1 = round(rexp(m, 1 / 5), 2),
                                         Annotation2 = round(runif(m, 0, 30), 2))
  if (g <= 2) {
    causal <- sample(m, round(frac_causal * m))
    beta <- numeric(m)
    beta[causal] <- if (g == 1) beta1 else beta2 * sample(c(-1, 1), length(causal), replace = TRUE)
    y_gene <- y_gene + as.numeric(G %*% beta)
  }
}
save(genotype_list, annotation_list, file = "demo/data/genotype_demo.RData")

## ---- phenotypes ----
fam <- rnorm(n_pair, 0, 0.6)                       # shared sibling-pair effect
y <- 0.02 * (covar$age - 55) + 0.3 * covar$sex + as.numeric(pcs %*% rep(0.1, 5)) + y_gene + rnorm(N)
y[i1] <- y[i1] + fam; y[i2] <- y[i2] + fam
s_pheno <- y + rnorm(N, 0, sd(y) * sqrt(1 / cor_pred^2 - 1))   # cor(s_pheno, y) ~ cor_pred
R <- rbinom(N, 1, pi_obs)
pheno <- data.table(ID = ids, pheno = round(ifelse(R == 1, y, NA), 4), s_pheno = round(s_pheno, 4))

fwrite(pheno, "demo/data/pheno_demo.txt", sep = " ", na = "NA", quote = FALSE)
fwrite(covar, "demo/data/covar_demo.txt", sep = " ")
fwrite(data.table(ids), "demo/data/sample_ids_demo.txt", col.names = FALSE)
cat(sprintf("N = %d, observed = %d (%.1f%%), cor(s_pheno, y) = %.3f\n", N, sum(R), 100 * mean(R), cor(s_pheno, y)))
