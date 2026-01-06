############################################################
# Boyault et al. 2007 — Molecular subclasses
# TCGA-LIHC overall survival analysis
#
# Biological context:
# Boyault et al. (2007) defined molecular subclasses of
# hepatocellular carcinoma (G1–G6) with distinct biology
# and clinical behavior. G1 tumors are associated with
# high proliferation and poorer prognosis.
#
# This script assigns TCGA-LIHC tumors to Boyault subclasses
# based on dominant expression scores and evaluates
# overall survival differences.
#
# Data sources:
# - TCGA-LIHC RNA-seq (TPM)
# - TCGA curated survival data
# - Boyault subclass gene lists
#
# Output:
# - Kaplan–Meier survival plot by Boyault subclass
############################################################

library(survival)
library(survminer)
library(readr)
library(readxl)
library(dplyr)
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

boyault_gene_dir <- file.path(
  base_dir,
  "data/processed/boyault2007"
)

out_png <- file.path(
  base_dir,
  "figures/boyault2007/Boyault2007_TCGA_LIHC_subclass_survival.png"
)

# -----------------------------
# 2. Load TCGA survival data
# -----------------------------
surv_data <- read_excel(tcga_surv_file)

surv_data$sample  <- substr(surv_data$sample, 1, 12)
surv_data$OS      <- as.numeric(surv_data$OS)
surv_data$OS.time <- as.numeric(surv_data$OS.time)

surv_data <- surv_data %>%
  dplyr::select(sample, OS, OS.time) %>%
  filter(!is.na(OS), !is.na(OS.time))

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

# Standardize symbols (important for Boyault lists)
rownames(expr_mat) <- toupper(trimws(rownames(expr_mat)))

# -----------------------------
# 5. Load Boyault gene lists
# -----------------------------
boyault_files <- list(
  G1 = "G1.xlsx",
  G2 = "G2.xlsx",
  G3 = "G3.xlsx",
  G5 = "G5.xlsx",
  G6 = "G6.xlsx"
)

boyault_genes <- lapply(boyault_files, function(f) {
  genes <- read_excel(file.path(boyault_gene_dir, f))[[1]]
  genes <- toupper(trimws(genes))
  unique(na.omit(genes))
})

# Keep only genes present in TCGA
boyault_genes <- lapply(boyault_genes, intersect, y = rownames(expr_mat))

# Diagnose overlap
overlap_counts <- sapply(boyault_genes, length)
print(overlap_counts)

# Keep subclasses with sufficient overlap
boyault_genes <- boyault_genes[overlap_counts >= 5]

if (length(boyault_genes) < 2) {
  stop("Too few Boyault subclasses overlap with TCGA data")
}

# -----------------------------
# 6. Score TCGA tumors
# -----------------------------
score_matrix <- sapply(boyault_genes, function(genes) {
  colMeans(expr_mat[genes, , drop = FALSE], na.rm = TRUE)
})

# -----------------------------
# 7. Assign dominant Boyault subclass
# -----------------------------
dominant_class <- apply(score_matrix, 1, function(x) {
  names(which.max(x))
})

surv_data <- surv_data %>%
  filter(sample %in% names(dominant_class)) %>%
  mutate(Boyault_Class = dominant_class[sample])

surv_data$Boyault_Class <- factor(
  surv_data$Boyault_Class,
  levels = names(boyault_genes)
)

# -----------------------------
# 8. Kaplan–Meier survival
# -----------------------------
surv_obj <- Surv(
  time  = surv_data$OS.time,
  event = surv_data$OS
)

fit <- survfit(surv_obj ~ Boyault_Class, data = surv_data)

km_plot <- ggsurvplot(
  fit,
  data = surv_data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  legend.title = "Boyault subclass",
  risk.table.height = 0.30,
  ggtheme = theme_bw() +
    theme(
      legend.position = "top",
      plot.margin = margin(10, 20, 10, 40)
    )
)

# -----------------------------
# 9. Save figure
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
