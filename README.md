# SPARTA
A powerful and computationally efficient framework for rare variant association testing that leverages machine-learning-predicted phenotypes to boost power in large-scale biobank studies.

SPARTA combines an **observed** phenotype (measured in only a subset of samples) with an **AI/ML-predicted** phenotype (available for *all* samples) into a debiased *pseudo phenotype*, fits a STAAR null model on it, and performs gene-centric rare variant tests with a two-normal-mixture saddlepoint approximation (SPA) that keeps the type-I error calibrated. This repository uses local file paths; for the UK Biobank Research Analysis Platform (DNANexus) version, see [`SPARTA-RAP`](https://github.com/shuangsong0110/SPARTA-RAP).


## Table of contents

* [Installation](#package-step-1-installation)
* [Prepare the Predicted Phenotype](#dna-step-2-prepare-the-predicted-phenotype)
* [Fit the Null Model](#scroll-step-3-fit-the-null-model)
* [Gene-centric Rare Variant Tests](#rocket-step-4-gene-centric-rare-variant-tests)
* [Summary and Figures](#framed_picture-step-5-summary-and-figures)


## :package: Step 1: Installation

SPARTA is built on [STAAR](https://github.com/xihaoli/STAAR) and [STAARpipeline](https://github.com/xihaoli/STAARpipeline).

**Recommended (Docker).** All required R packages are pre-installed in the STAARpipeline Docker image, pulled from Docker Hub:
```
docker pull zilinli/staarpipeline:0.9.9
```

**Manual install.** Alternatively, install the dependencies in R:
```
## CRAN
install.packages(c("data.table","Matrix","dplyr","devtools"))

## Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("gdsfmt","SeqArray","SeqVarTools"))

## GitHub (STAAR ecosystem)
devtools::install_github("xihaoli/STAAR")
devtools::install_github("xihaoli/STAARpipeline")
devtools::install_github("xihaoli/STAARpipelineSummary")
```

The scripts load the following at run time:
```
library(data.table)
library(STAAR)
library(SeqArray)
library(SeqVarTools)
library(gdsfmt)
library(STAARpipeline)
library(STAARpipelineSummary)
library(STAARpipelinePheWAS)
```
Note: `genes_info` (the genome-wide gene list used to partition jobs) ships with `STAARpipeline`; all packages above are pre-installed in the Docker image.


## :dna: Step 2: Prepare the predicted phenotype

The only SPARTA-specific input is a per-phenotype table of **observed** and **predicted** phenotype values. The predicted phenotype can be trained with any generative-AI / machine-learning model of your choice (e.g. a random forest, gradient boosting, or a deep network) using whatever features are available for the full cohort.

**Format:** a headed, whitespace-delimited `data.frame` with three columns:

| ID | pheno | s_pheno |
|----|-------|---------|
| sample_1 | 1.24 | 0.68 |
| sample_2 | -0.53 | -0.11 |
| sample_3 | NA | 0.42 |
| ... | ... | ... |

*(The rows above are synthetic and for illustration only — not real sample IDs or phenotype data.)*

- `ID`: sample id.
- `pheno`: the **observed** phenotype; set to `NA` for samples in which it was not measured.
- `s_pheno`: the **predicted** phenotype, provided for **all** samples (observed and unobserved).

SPARTA internally defines the observation indicator `R = 1` if `pheno` is observed and `R = 0` otherwise, and builds the pseudo phenotype

```
pheno_pseudo = s_pheno + (1 / mean(R)) * R * (pheno - s_pheno)
```


## :scroll: Step 3: Fit the null model

### Step 3.1 Sparse Genetic Relatedness Matrix (GRM) — *optional*

If the study contains related samples, model the relatedness with a sparse GRM. The R package [**FastSparseGRM**](https://github.com/rounakdey/FastSparseGRM) efficiently computes genetic principal components and an ancestry-adjusted sparse GRM; the PCs can be used as fixed-effect covariates and the sparse GRM as the random effect in the null model. See the [FastSparseGRM](https://github.com/rounakdey/FastSparseGRM) repository and the [STAARpipeline-Tutorial](https://github.com/xihaoli/STAARpipeline-Tutorial#generate-sparse-genetic-relatedness-matrix-grm) for details. The output is an `.RData` file containing a sparse matrix `sGRM`.

If the samples are **unrelated**, this step can be omitted (set `grm_file <- NULL` in Step 3.2).

### Step 3.2 SPARTA null model

#### Script: <a href="R/SPARTA_null_model.R">**SPARTA_null_model.R**</a>
Builds the pseudo phenotype from the observed + predicted phenotypes (Step 2), fits a STAAR null model adjusting for covariates (and, optionally, the sparse GRM), and stores the two SPARTA scaling factors used by the SPA correction.

**Input:** the predicted-phenotype file (Step 2), a covariate file (keyed by sample id), the sequenced-sample id list, and (optionally) the sparse GRM.<br>
**Output:** one `.rda` per phenotype (`obj.STAAR_500K_<field>_norm_pseudo_rf100_v5.rda`) storing:
- `obj_nullmodel` — the fitted STAAR null model,
- `resid_ratio` — the ratio of residual variance in observed (`R=1`) vs unobserved (`R=0`) samples,
- `pi1` — the fraction of samples with an observed phenotype.

Edit the `User Input` block at the top of the script, then run:
```
Rscript R/SPARTA_null_model.R
```


## :rocket: Step 4: Gene-centric rare variant tests

Both scripts partition the genome into **195 jobs** (blocks of 100 genes per chromosome). The array index `arrayid` (1..195) selects one block; results are written per job. Edit the `User Input` block (aGDS directory, null-model directory, annotation catalog, output path) before running.

### Step 4.1 Gene-centric coding regions

#### Script: <a href="R/SPARTA_coding.R">**SPARTA_coding.R**</a>
Performs the SPARTA gene-centric coding analysis. Five functional categories aggregate the coding rare variants of each protein-coding gene: (1) putative loss of function (stop gain, stop loss, splice), (2) missense, (3) disruptive missense, (4) putative loss of function + disruptive missense, and (5) synonymous.
```
Rscript R/SPARTA_coding.R ${arrayid}      # arrayid = 1..195
```
#### Output: one `Coding_<kk>.Rdata` file per gene block, written to `output_path`.

### Step 4.2 Gene-centric non-coding regions

#### Script: <a href="R/SPARTA_noncoding.R">**SPARTA_noncoding.R**</a>
Performs the SPARTA gene-centric noncoding analysis. Seven regulatory categories aggregate the noncoding rare variants: (1) promoter RVs overlaid with CAGE sites, (2) promoter RVs overlaid with DHS sites, (3) enhancer RVs overlaid with CAGE sites, (4) enhancer RVs overlaid with DHS sites, (5) untranslated region (UTR) RVs, (6) upstream RVs, and (7) downstream RVs.
```
Rscript R/SPARTA_noncoding.R ${arrayid}   # arrayid = 1..195
```
#### Output: one `Noncoding_<kk>.Rdata` file per gene block, written to `output_path`.

Note: Hardy-Weinberg equilibrium filtering is assumed to be handled upstream in variant QC (`QC_label`), so no separate HWE step is required. The noncoding scripts read a precomputed per-gene noncoding variant-id file (`Noncoding_<kk>.Rdata`, one variant set per regulatory category) as a prerequisite input (path set in the `User Input` block).


## :framed_picture: Step 5: Summary and figures

#### Script: <a href="R/summary_coding.R">**summary_coding.R**</a> and <a href="R/summary_noncoding.R">**summary_noncoding.R**</a>
Summarize the gene-centric results from Step 4. Each script sources its summary function (`R/Gene_Centric_Coding_Results_Summary_anno.R` / `R/Gene_Centric_Noncoding_Results_Summary_anno.R`), aggregates the per-job result files for each phenotype, and produces genome-wide tables and figures.

**Input:** paths to the Step 4 per-job result files and an output path (set in the `User Input` block).<br>
**Output:** a Manhattan plot, a Q-Q plot, and CSVs of the significant coding / noncoding masks.
```
Rscript R/summary_coding.R
Rscript R/summary_noncoding.R
```


## Data Availability
The whole-genome functional annotation data and the precomputed annotation principal components are available at the [Functional Annotation of Variant - Online Resource (FAVOR)](https://favor.genohub.org) site and [FAVOR Essential Database](https://doi.org/10.7910/DVN/1VGTJI).


## Citation
If you use SPARTA, please cite:

Shuang Song, Ruoyu Wang, Haoyu Yang, Max Elliott, Xihao Li, Jingnan Du, Randy Lee Buckner, Xihong Lin. *Augmented synthetic phenotype analysis enhances rare variant association discovery in biobank-scale whole genome sequencing studies.*


## :busts_in_silhouette: Maintainer
Please contact Shuang Song (shuangsong@hsph.harvard.edu) if there are any problems or questions.


## References
[STAAR](https://github.com/xihaoli/STAAR/)
[STAARpipeline](https://github.com/xihaoli/STAARpipeline/)
[STAARpipelineSummary](https://github.com/xihaoli/STAARpipelineSummary/)
[FastSparseGRM](https://github.com/rounakdey/FastSparseGRM/)
