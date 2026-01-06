############################################################
# Lee et al. 2004 — Proliferation signature
# TCGA-LIHC overall survival analysis
#
# Biological context:
# Lee et al. (2004) identified a proliferation-driven
# subclass of hepatocellular carcinoma associated with
# aggressive tumor behavior and poor prognosis.
#
# This script computes a proliferation score using
# canonical cell-cycle genes and evaluates its association
# with overall survival in TCGA-LIHC.
#
# Data sources:
# - TCGA-LIHC RNA-seq (TPM)
# - TCGA curated survival data
#
# Output:
# - Kaplan–Meier survival plot
############################################################

library(survival)
library(survminer)
library(readr)
library(readxl)
library(org.Hs.eg.db)
library(AnnotationDbi)

# -----------------------------
# 1. Paths (EDIT THESE)
# -----------------------------
base_dir <- "PATH_TO_PROJECT_ROOT"

tcga_expr_file <- file.path(
  base_dir,
  "data/raw/TCGA_LIHC/TCGA-LIHC.star_tpm.tsv"
)

tcga_surv_file <- file.path(
  base_dir,
  "data/raw/TCGA_LIHC/phenotype_curated_survival.xlsx"
)

out_png <- file.path(
  base_dir,
  "figures/lee2004/Lee2004_TCGA_LIHC_proliferation_survival.png"
)

# -----------------------------
# 2. Load survival data
# -----------------------------
surv_data <- read_excel(tcga_surv_file)

surv_data$sample  <- substr(surv_data$sample, 1, 12)
surv_data$OS      <- as.numeric(surv_data$OS)
surv_data$OS.time <- as.numeric(surv_data$OS.time)

surv_data <- surv_data[
  complete.cases(surv_data[, c("sample", "OS", "OS.time")]),
  c("sample", "OS", "OS.time")
]

# -----------------------------
# 3. Load TCGA expression (TPM)
# -----------------------------
expr_data <- read_tsv(tcga_expr_file, show_col_types = FALSE)

expr_mat <- as.matrix(expr_data[, -1])
rownames(expr_mat) <- expr_data[[1]]
mode(expr_mat) <- "numeric"

expr_mat <- log2(expr_mat + 1)
colnames(expr_mat) <- substr(colnames(expr_mat), 1, 12)

# -----------------------------
# 4. Ensembl → gene symbol mapping
# -----------------------------
ensembl_ids <- sub("\\..*$", "", rownames(expr_mat))

symbol_map <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

keep <- !is.na(symbol_map)
expr_mat <- expr_mat[keep, ]
rownames(expr_mat) <- symbol_map[keep]

expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]

# -----------------------------
# 5. Lee 2004 proliferation genes
# -----------------------------
lee_prolif_genes <- c(
  "MKI67","PCNA","TOP2A","AURKA","AURKB",
  "CDC20","CDK1","CDK2","CCNB1","CCNA2",
  "MCM2","MCM3","MCM4","MCM5","MCM6","MCM7",
  "BUB1","BUB1B","UBE2C","FOXM1"
)

lee_prolif_genes <- intersect(lee_prolif_genes, rownames(expr_mat))

if (length(lee_prolif_genes) < 10) {
  stop("Proliferation gene mapping failed")
}

# -----------------------------
# 6. Proliferation score
# -----------------------------
prolif_score <- colMeans(
  expr_mat[lee_prolif_genes, , drop = FALSE],
  na.rm = TRUE
)

prolif_score <- scale(prolif_score)[,1]

# -----------------------------
# 7. Define tumor groups
# -----------------------------
lee_group <- ifelse(
  prolif_score >= median(prolif_score, na.rm = TRUE),
  "Proliferation-high tumors",
  "Proliferation-low tumors"
)

names(lee_group) <- names(prolif_score)

# -----------------------------
# 8. Merge with survival
# -----------------------------
surv_data <- surv_data[surv_data$sample %in% names(lee_group), ]
surv_data$Lee_Group <- lee_group[surv_data$sample]

surv_data$Lee_Group <- factor(
  surv_data$Lee_Group,
  levels = c(
    "Proliferation-high tumors",
    "Proliferation-low tumors"
  )
)

# -----------------------------
# 9. Kaplan–Meier survival
# -----------------------------
surv_obj <- Surv(
  time  = surv_data$OS.time,
  event = surv_data$OS
)

fit <- survfit(surv_obj ~ Lee_Group, data = surv_data)

km_plot <- ggsurvplot(
  fit,
  data = surv_data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  legend.title = "Tumor biology (Lee 2004)",
  palette = c("#D55E00", "#0072B2"),
  risk.table.height = 0.30,
  ggtheme = theme_bw()
)

# -----------------------------
# 10. Save figure
# -----------------------------
ggsave(
  filename = out_png,
  plot = arrange_ggsurvplots(list(km_plot), print = FALSE),
  width  = 14,
  height = 10,
  dpi = 300
)

############################################################
# END
############################################################
