############################################################
# Chiang-style tumor biology
# CTNNB1-like vs Proliferation-like tumors
# TCGA-LIHC overall survival analysis
#
# Biological context:
# Chiang et al. described liver tumor biology driven by
# either CTNNB1 (beta-catenin) signaling or high
# proliferation programs. These tumor types show
# distinct clinical outcomes.
#
# This script compares CTNNB1-driven versus
# proliferation-driven tumors using MSigDB Hallmark
# gene sets and evaluates overall survival in TCGA-LIHC.
#
# Data sources:
# - TCGA-LIHC RNA-seq (TPM)
# - TCGA curated survival data
# - MSigDB Hallmark gene sets
#
# Output:
# - Kaplan–Meier survival plot
############################################################

library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(survival)
library(survminer)
library(readr)
library(readxl)
library(dplyr)

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
  "figures/chiang/Chiang_TCGA_LIHC_CTNNB1_vs_Proliferation.png"
)

# -----------------------------
# 2. Load TCGA expression
# -----------------------------
expr_data <- read_tsv(tcga_expr_file, show_col_types = FALSE)

expr_mat <- as.matrix(expr_data[, -1])
rownames(expr_mat) <- expr_data[[1]]
mode(expr_mat) <- "numeric"

expr_mat <- log2(expr_mat + 1)
colnames(expr_mat) <- substr(colnames(expr_mat), 1, 12)

# -----------------------------
# 3. Ensembl → gene symbol mapping
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
# 4. Load TCGA survival data
# -----------------------------
surv_data <- read_excel(tcga_surv_file)

surv_data$sample  <- substr(surv_data$sample, 1, 12)
surv_data$OS      <- as.numeric(surv_data$OS)
surv_data$OS.time <- as.numeric(surv_data$OS.time)

surv_data <- surv_data %>%
  dplyr::select(sample, OS, OS.time) %>%
  filter(!is.na(OS), !is.na(OS.time))

# -----------------------------
# 5. Load MSigDB Hallmark programs
# -----------------------------
msig <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

ctnnb1_genes <- msig %>%
  filter(gs_name == "HALLMARK_WNT_BETA_CATENIN_SIGNALING") %>%
  pull(gene_symbol) %>%
  intersect(rownames(expr_mat))

prolif_genes <- msig %>%
  filter(gs_name %in% c("HALLMARK_E2F_TARGETS",
                        "HALLMARK_G2M_CHECKPOINT")) %>%
  pull(gene_symbol) %>%
  intersect(rownames(expr_mat))

if (length(ctnnb1_genes) < 10 || length(prolif_genes) < 30) {
  stop("Insufficient gene overlap for Chiang-style programs")
}

# -----------------------------
# 6. Score tumors
# -----------------------------
ctnnb1_score <- colMeans(expr_mat[ctnnb1_genes, , drop = FALSE])
prolif_score <- colMeans(expr_mat[prolif_genes, , drop = FALSE])

score_df <- data.frame(
  sample = names(ctnnb1_score),
  CTNNB1 = ctnnb1_score,
  Proliferation = prolif_score
)

score_df$Chiang_Class <- ifelse(
  score_df$Proliferation > score_df$CTNNB1,
  "Proliferation-like tumors",
  "CTNNB1-like tumors"
)

# -----------------------------
# 7. Merge with survival
# -----------------------------
surv_data <- surv_data %>%
  filter(sample %in% score_df$sample) %>%
  left_join(score_df, by = "sample")

# -----------------------------
# 8. Kaplan–Meier survival
# -----------------------------
surv_obj <- Surv(surv_data$OS.time, surv_data$OS)
fit <- survfit(surv_obj ~ Chiang_Class, data = surv_data)

km_plot <- ggsurvplot(
  fit,
  data = surv_data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  legend.title = "Tumor biology",
  legend.labs = c("CTNNB1-like", "Proliferation-like"),
  palette = c("#0072B2", "#D55E00"),
  risk.table.height = 0.30,
  ggtheme = theme_bw()
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
