############################################################
# Hoshida et al. 2009 — S1 subclass
# TCGA-LIHC overall survival analysis
#
# Biological context:
# Hoshida et al. (2009) defined molecular subclasses of
# hepatocellular carcinoma. The S1 subclass is associated
# with inflammatory signaling and poorer prognosis.
#
# This script computes an S1 score per tumor and evaluates
# its association with overall survival in TCGA-LIHC.
#
# Data sources:
# - TCGA-LIHC RNA-seq (TPM)
# - TCGA curated survival data
#
# Output:
# - Kaplan–Meier survival plot (S1-like vs non–S1-like)
############################################################

library(dplyr)
library(readr)
library(readxl)
library(survival)
library(survminer)
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

hoshida_gene_file <- file.path(
  base_dir,
  "data/processed/hoshida2009_s1_genes.xlsx"
)

out_png <- file.path(
  base_dir,
  "figures/hoshida2009/Hoshida2009_TCGA_LIHC_S1_survival.png"
)

# -----------------------------
# 2. Load Hoshida S1 genes
# -----------------------------
hoshida <- read_excel(hoshida_gene_file)

stopifnot(all(c("Subtype", "Symbol") %in% colnames(hoshida)))

s1_genes <- unique(hoshida$Symbol[hoshida$Subtype == "S1"])
s1_genes <- s1_genes[!is.na(s1_genes)]

if (length(s1_genes) < 30) {
  stop("Too few S1 genes loaded")
}

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
expr_mat <- expr_mat[keep, , drop = FALSE]
rownames(expr_mat) <- symbol_map[keep]

expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]

# -----------------------------
# 5. Keep S1 genes present in TCGA
# -----------------------------
s1_genes <- intersect(s1_genes, rownames(expr_mat))

if (length(s1_genes) < 20) {
  stop("Too few S1 genes overlap with TCGA")
}

# -----------------------------
# 6. Compute S1 score per tumor
# -----------------------------
s1_score <- colMeans(
  expr_mat[s1_genes, , drop = FALSE],
  na.rm = TRUE
)

# -----------------------------
# 7. Load survival data
# -----------------------------
surv_data <- read_excel(tcga_surv_file)

surv_data$sample  <- substr(surv_data$sample, 1, 12)
surv_data$OS      <- as.numeric(surv_data$OS)
surv_data$OS.time <- as.numeric(surv_data$OS.time)

surv_data <- surv_data %>%
  dplyr::select(sample, OS, OS.time) %>%
  filter(!is.na(OS), !is.na(OS.time))

# Match patients
surv_data <- surv_data %>%
  filter(sample %in% names(s1_score)) %>%
  mutate(s1_score = s1_score[sample])

# -----------------------------
# 8. Define S1-like vs non–S1-like
# -----------------------------
cutoff <- median(surv_data$s1_score, na.rm = TRUE)

surv_data <- surv_data %>%
  mutate(
    group = ifelse(
      s1_score >= cutoff,
      "S1-like tumors",
      "Non–S1-like tumors"
    )
  )

# -----------------------------
# 9. Kaplan–Meier survival
# -----------------------------
surv_obj <- Surv(
  time  = surv_data$OS.time,
  event = surv_data$OS
)

fit <- survfit(surv_obj ~ group, data = surv_data)

km_plot <- ggsurvplot(
  fit,
  data = surv_data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  legend.title = "",
  legend.labs = c("S1-like tumors", "Non–S1-like tumors"),
  xlab = "Time (days)",
  ylab = "Overall survival probability",
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
