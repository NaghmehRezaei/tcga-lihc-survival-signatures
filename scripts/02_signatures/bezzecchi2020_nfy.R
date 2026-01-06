############################################################
# Bezzecchi et al. 2020 — NF-Y–driven tumor biology
# TCGA-LIHC overall survival analysis
#
# Biological context:
# Bezzecchi et al. (2020) identified a proliferative
# hepatocellular carcinoma subtype driven by NF-Y activity,
# distinct from metabolically active tumors.
#
# This script computes an NF-Y–weighted expression score
# and evaluates its association with overall survival in
# TCGA-LIHC.
#
# Data sources:
# - TCGA-LIHC RNA-seq (TPM)
# - TCGA curated survival data
# - Bezzecchi DEG signature (UP/DOWN)
#
# Output:
# - Kaplan–Meier survival plot
############################################################

library(survival)
library(survminer)
library(readxl)
library(readr)
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

bezzecchi_deg_file <- file.path(
  base_dir,
  "data/processed/bezzecchi2020/Table_S3_DEG.xlsx"
)

out_png <- file.path(
  base_dir,
  "figures/bezzecchi2020/Bezzecchi2020_TCGA_LIHC_NFY_survival.png"
)

# -----------------------------
# 2. Load survival data
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

# -----------------------------
# 5. Load Bezzecchi DEG signature
# -----------------------------
deg <- read_excel(bezzecchi_deg_file)

bezz_up   <- unique(na.omit(deg$UP))
bezz_down <- unique(na.omit(deg$DOWN))

bezz_up   <- intersect(bezz_up, rownames(expr_mat))
bezz_down <- intersect(bezz_down, rownames(expr_mat))

if (length(bezz_up) < 20 || length(bezz_down) < 20) {
  stop("Insufficient Bezzecchi DEG overlap with TCGA")
}

# -----------------------------
# 6. Compute NF-Y weighted score
# -----------------------------
sig_genes <- c(bezz_up, bezz_down)
sig_mat <- expr_mat[sig_genes, , drop = FALSE]

gene_weights <- c(
  setNames(rep( 1, length(bezz_up)),   bezz_up),
  setNames(rep(-1, length(bezz_down)), bezz_down)
)

weighted_expr <- sweep(
  sig_mat,
  1,
  gene_weights[rownames(sig_mat)],
  `*`
)

nfy_score <- colMeans(weighted_expr, na.rm = TRUE)

# -----------------------------
# 7. Merge with survival
# -----------------------------
surv_data <- surv_data %>%
  filter(sample %in% names(nfy_score)) %>%
  mutate(score = nfy_score[sample])

# -----------------------------
# 8. Define proliferative vs metabolic tumors
# -----------------------------
q_hi <- quantile(surv_data$score, 0.75, na.rm = TRUE)
q_lo <- quantile(surv_data$score, 0.25, na.rm = TRUE)

surv_data <- surv_data %>%
  filter(score >= q_hi | score <= q_lo) %>%
  mutate(
    group = ifelse(
      score >= q_hi,
      "Proliferative-like tumors (NF-Y–high)",
      "Metabolic-like tumors (NF-Y–low)"
    )
  )

# -----------------------------
# 9. Kaplan–Meier survival
# -----------------------------
surv_obj <- Surv(surv_data$OS.time, surv_data$OS)
fit <- survfit(surv_obj ~ group, data = surv_data)

km_plot <- ggsurvplot(
  fit,
  data = surv_data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  legend.title = "",
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
