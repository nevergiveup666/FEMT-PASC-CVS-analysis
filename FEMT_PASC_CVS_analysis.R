# ====================================================
# Network Pharmacology Analysis of Five-Element Music
# Therapy (FEMT) for PASC-CVS
# ====================================================
# Runtime environment
#   R version: 4.2.2
#   Key packages: limma 3.50.3, clusterProfiler 4.4.4, pROC 1.18.0,
#                 caret 6.0-94, igraph 1.3.5, org.Hs.eg.db 3.16
#
# Input files (place in working directory):
#   GSE171110_Data_RNAseq_norm_counts_geo.txt
#   GSE157103_genes.tpm.tsv
#   GSE171110_series_matrix.txt   (used to derive sample groups)
#
# Outputs:
#   ./results/                  figures (PDF/PNG) and RData
#   ./results/10_tables/        supplementary tables (CSV)
#
# This script is provided for reproducibility of the manuscript.
# All outputs are hypothesis-generating; no causal inference is made.
# ====================================================

rm(list = ls())
gc()

# ---- Packages ----
library(tidyverse)
library(limma)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(igraph)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
library(pROC)
library(randomForest)
library(caret)
library(glmnet)
library(ggvenn)
library(gridExtra)
library(grid)

# ---- Paths ----
setwd("...")  # edit before upload
output_dir <- "..."
if (!dir.exists(output_dir)) dir.create(output_dir)
if (!dir.exists(file.path(output_dir, "...")))
  dir.create(file.path(output_dir, "..."))

# ---- Helpers ----
print_line <- function(char = "=", length = 50) {
  cat(paste(rep(char, length), collapse = ""), "\n")
}
print_title <- function(title) {
  print_line("=", 60); cat(title, "\n"); print_line("=", 60)
}

# ====================================================
# 1. FEMT target database
# ====================================================
create_enhanced_FEMT_database <- function() {
  cat("\n========================================\n")
  cat("Creating comprehensive FEMT database\n")
  cat("========================================\n")
  
  literature_evidence <- data.frame(
    Target_Class = c(
      "NLRP3 inflammasome pathway",
      "BDNF/FNDC5 pathway",
      "Endorphin/mu opioid receptor system",
      "Cardiovascular regulation",
      "Inflammatory & oxidative stress markers",
      "Calcium signaling pathway",
      "Immune modulation targets"
    ),
    Specific_Genes = c(
      "NLRP3, IL1B, PYCARD, CASP1",
      "BDNF, NTRK2, FNDC5, CREB1",
      "POMC, OPRM1, PENK, PDYN",
      "NOS3, ACE, AGTR1, EDN1, ADRB1, ADRB2, CHRNA7",
      "NOS2, SOD1, CAT, GPX1, IL6, TNF, IL10, IL4",
      "CAMK2A, CAMK2B, CAMK2D, CAMK2G",
      "IL2RB, TUBA1B, CD4, CD8A, FOXP3, CD19"
    ),
    Primary_Reference = c(
      "Cheng et al. 2020", "Hung et al. 2021", "Dan et al. 2024",
      "Gao et al. 2025", "Wang et al. 2025", "Lin et al. 2018",
      "Wu et al. 2025"
    ),
    Evidence_Level = c(
      "Strong (animal)", "Strong (animal)", "Moderate (review)",
      "Moderate (review)", "Moderate (clinical trial)",
      "Moderate (conference)", "Moderate (bioinformatics)"
    ),
    stringsAsFactors = FALSE
  )
  
  femt_targets <- list(
    Zhi = list(
      targets = c(
        "NOS3", "ACE", "AGTR1", "EDN1", "NPPA", "NPPB",
        "ADRB1", "ADRB2", "CHRM2", "CHRM3", "CHRNA7",
        "HTR1A", "HTR2A", "HTR2C", "DRD1", "DRD2",
        "NOS2", "SOD1", "SOD2", "CAT", "GPX1", "GPX4",
        "IL6", "TNF", "IL1B"
      ),
      functions = c("Cardiovascular regulation", "Autonomic nervous system",
                    "Mood-neural coupling", "Oxidative stress response",
                    "Inflammatory response"),
      references = c("Gao et al. 2025", "Wang et al. 2025", "Dan et al. 2024")
    ),
    Yu = list(
      targets = c(
        "NR3C1", "FKBP5", "CRH", "CRHR1", "AVP", "OXT",
        "POMC", "OPRM1", "PENK", "PDYN",
        "BDNF", "NTRK2", "FNDC5", "CREB1",
        "AQP2", "SLC12A3", "KCNJ1"
      ),
      functions = c("Stress response (HPA axis)", "Endogenous opioid system",
                    "Neuroprotection", "Fluid balance"),
      references = c("Dan et al. 2024", "Hung et al. 2021")
    ),
    Jue = list(
      targets = c(
        "MAOA", "MAOB", "COMT", "SLC6A4",
        "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G",
        "CALM1", "CALM2", "CALM3",
        "CYP3A4", "CYP2D6", "GSTP1", "UGT1A1", "CYP1A1", "CYP1A2",
        "GABRA1", "GABRA2", "GABRB3", "GABRG2"
      ),
      functions = c("Neurotransmitter metabolism", "Calcium signaling",
                    "Detoxification", "Inhibitory neurotransmission"),
      references = c("Dan et al. 2024", "Lin et al. 2018")
    ),
    Shang = list(
      targets = c(
        "NLRP3", "IL1B", "PYCARD", "CASP1",
        "IL6", "TNF", "TGFB1", "IL10", "IL4", "IL13", "IFNG",
        "ACE2", "TMPRSS2", "CTSL", "CTSB",
        "IL2RB", "TUBA1B", "CD4", "CD8A", "FOXP3", "CD19",
        "SFTPC", "SFTPA1", "SFTPB", "SFTPD"
      ),
      functions = c("NLRP3 inflammasome pathway", "Inflammatory response",
                    "Viral entry mechanisms", "Immune modulation",
                    "Surfactant production"),
      references = c("Cheng et al. 2020", "Wang et al. 2025", "Wu et al. 2025")
    ),
    Gong = list(
      targets = c(
        "BDNF", "NTRK2", "FNDC5", "CREB1",
        "SNAP25", "SYN1", "SYN2", "DLG4",
        "SHANK1", "SHANK2", "SHANK3",
        "NLGN1", "NLGN2", "NLGN3",
        "FOS", "JUN", "EGR1", "EGR2", "EGR3", "ARC",
        "DCX", "MAP2", "NES", "SOX2", "TUBB3",
        "INS", "IGF1", "LEP", "ADIPOQ", "GHRL", "GCG"
      ),
      functions = c("BDNF signaling pathway", "Synaptic plasticity",
                    "Immediate early gene response", "Neurogenesis",
                    "Metabolic regulation"),
      references = c("Hung et al. 2021", "Wang et al. 2025")
    )
  )
  
  femt_df <- tibble()
  for (tone in names(femt_targets)) {
    organ <- switch(tone, Zhi = "Heart", Yu = "Kidney", Jue = "Liver",
                    Shang = "Lung", Gong = "Spleen")
    element <- switch(tone, Zhi = "Fire", Yu = "Water", Jue = "Wood",
                      Shang = "Metal", Gong = "Earth")
    targets_vec <- unique(femt_targets[[tone]]$targets)
    functions_vec <- rep(paste(femt_targets[[tone]]$functions, collapse = "; "),
                         length.out = length(targets_vec))
    refs_vec <- rep(paste(femt_targets[[tone]]$references, collapse = "; "),
                    length.out = length(targets_vec))
    femt_df <- bind_rows(femt_df, tibble(
      Tone = tone, Element = element, Organ = organ,
      Target = targets_vec, Function = functions_vec, References = refs_vec
    ))
  }
  femt_df <- femt_df %>% distinct(Target, .keep_all = TRUE)
  
  suppressMessages({
    entrez_ids <- mapIds(org.Hs.eg.db, keys = femt_df$Target,
                         column = "ENTREZID", keytype = "SYMBOL",
                         multiVals = "first")
    femt_df$EntrezID <- as.character(entrez_ids)
    uniprot_ids <- mapIds(org.Hs.eg.db, keys = femt_df$Target,
                          column = "UNIPROT", keytype = "SYMBOL",
                          multiVals = "first")
    femt_df$UniProtID <- uniprot_ids
  })
  femt_df <- femt_df %>% filter(!is.na(EntrezID))
  
  strong_evidence_genes <- c("NLRP3", "IL1B", "BDNF", "FNDC5", "NOS2", "IL6",
                             "TNF", "POMC", "OPRM1", "NR3C1", "FKBP5", "CRH",
                             "FOS", "JUN")
  moderate_evidence_genes <- c("NOS3", "ACE", "AGTR1", "EDN1", "SOD1", "CAT",
                               "GPX1", "CAMK2A", "CAMK2B", "COMT", "MAOA",
                               "MAOB", "SLC6A4", "ADRB1", "ADRB2", "CHRNA7",
                               "HTR1A", "HTR2A", "DRD1", "DRD2", "SNAP25",
                               "SYN1", "DLG4", "EGR1", "ARC", "IL10", "IL4",
                               "IL2RB", "TUBA1B")
  femt_df <- femt_df %>%
    mutate(Evidence_Level = case_when(
      Target %in% strong_evidence_genes ~ "Strong (animal/multi-omics)",
      Target %in% moderate_evidence_genes ~ "Moderate",
      TRUE ~ "Preliminary/Theoretical"
    ))
  
  pathway_mapping <- list(
    "NLRP3|IL1B|PYCARD|CASP1" = "NLRP3 inflammasome",
    "BDNF|NTRK2|FNDC5|CREB1" = "BDNF signaling",
    "POMC|OPRM1|PENK|PDYN" = "Opioid system",
    "NOS2|SOD|CAT|GPX" = "Oxidative stress",
    "IL6|TNF|TGFB1|IL10|IL4|IL13|IFNG" = "Cytokine signaling",
    "SNAP25|SYN|DLG4|SHANK|NLGN" = "Synaptic plasticity",
    "CD4|CD8A|FOXP3|CD19|IL2RB|TUBA1B" = "Immune modulation",
    "ACE|AGTR1|EDN1|NPP" = "Cardiovascular",
    "CAMK2|CALM" = "Calcium signaling",
    "HTR|DRD|SLC6A4|MAO|COMT" = "Neurotransmitter",
    "ADRB|CHRM|CHRNA" = "Autonomic nervous",
    "NR3C1|FKBP5|CRH|AVP|OXT" = "HPA axis",
    "DCX|MAP2|NES|SOX2|TUBB3" = "Neurogenesis",
    "FOS|JUN|EGR|ARC" = "Immediate early genes"
  )
  femt_df$Pathway <- "Other"
  for (pattern in names(pathway_mapping)) {
    idx <- grepl(pattern, femt_df$Target, ignore.case = TRUE)
    femt_df$Pathway[idx] <- pathway_mapping[[pattern]]
  }
  
  write.csv(femt_df,
            file.path(output_dir, "10_tables", "FEMT_database_complete.csv"),
            row.names = FALSE)
  write.csv(literature_evidence,
            file.path(output_dir, "10_tables",
                      "Literature_evidence_summary.csv"),
            row.names = FALSE)
  
  cat("\n=== FEMT Database Summary ===\n")
  cat(sprintf("Total unique targets: %d\n", nrow(femt_df)))
  cat("\nDistribution by Tone:\n"); print(table(femt_df$Tone))
  cat("\nEvidence Level Distribution:\n"); print(table(femt_df$Evidence_Level))
  return(femt_df)
}

# ====================================================
# 2. Dataset processing
# ====================================================
# Derive sample groups from GSE171110_series_matrix.txt.
# The expression matrix columns follow the series matrix sample order.
parse_series_matrix_groups <- function(series_file, n_samples) {
  if (!file.exists(series_file)) return(NULL)
  sm_lines <- readLines(series_file, n = 300, warn = FALSE)
  sm_lines <- sm_lines[grepl("^!Sample_", sm_lines)]
  
  get_field <- function(key) {
    ln <- sm_lines[grepl(paste0("^", key), sm_lines)]
    if (length(ln) == 0) return(NULL)
    vals <- strsplit(ln[1], "\t")[[1]][-1]
    gsub('^"|"$', '', vals)
  }
  
  titles <- get_field("!Sample_title")
  chars  <- get_field("!Sample_characteristics_ch1")
  
  src <- NULL
  if (!is.null(chars) && any(grepl("group|disease|condition", chars,
                                   ignore.case = TRUE))) {
    src <- chars
  } else if (!is.null(titles)) {
    src <- titles
  }
  if (is.null(src)) return(NULL)
  
  cat("  - Series matrix sample titles (first 5):\n")
  if (!is.null(titles)) print(head(titles, 5))
  
  grp_raw <- tolower(src)
  grp <- ifelse(grepl("healthy|control|^hv$|^hc$", grp_raw), "Control",
                ifelse(grepl("severe|covid|patient|case", grp_raw), "Severe", NA))
  
  if (length(grp) == n_samples) {
    cat("  - Group vector aligned by series-matrix sample order\n")
    return(grp)
  }
  warning("Series matrix sample count (", length(grp),
          ") differs from expression matrix columns (", n_samples, ")")
  NULL
}

process_GSE171110 <- function() {
  cat("\n========================================\n")
  cat("Processing main dataset GSE171110\n")
  cat("========================================\n")
  
  gse_file <- "GSE171110_Data_RNAseq_norm_counts_geo.txt"
  if (!file.exists(gse_file)) stop("Main dataset file not found: ", gse_file)
  
  gse_data <- read.delim(gse_file, check.names = FALSE)
  cat("Data dimensions:", dim(gse_data), "\n")
  
  genes <- gse_data[, 1]
  expr_matrix <- as.matrix(gse_data[, -1])
  expr_matrix <- apply(expr_matrix, 2, as.numeric)
  rownames(expr_matrix) <- genes
  cat("Expression matrix dimensions:", dim(expr_matrix), "\n")
  
  # ---- Group assignment ----
  group <- parse_series_matrix_groups(
    "GSE171110_series_matrix.txt", ncol(expr_matrix))
  
  if (is.null(group)) {
    warning("Falling back to name-based grouping. Verify manually.")
    group <- rep(NA_character_, ncol(expr_matrix))
    for (i in seq_along(colnames(expr_matrix))) {
      nm <- colnames(expr_matrix)[i]
      if (grepl("^V|healthy|control|^HC", nm, ignore.case = TRUE)) {
        group[i] <- "Control"
      } else if (grepl("^C|severe|covid", nm, ignore.case = TRUE)) {
        group[i] <- "Severe"
      } else if (grepl("^[0-9]+$", nm)) {
        group[i] <- "Severe"
      }
    }
  }
  
  cat("\nGroup distribution (before filtering):\n")
  print(table(group, useNA = "ifany"))
  
  keep <- !is.na(group)
  expr_matrix <- expr_matrix[, keep, drop = FALSE]
  group <- factor(group[keep], levels = c("Control", "Severe"))
  
  if (sum(group == "Control") == 0 || sum(group == "Severe") == 0)
    stop("Group assignment failed: one group has zero samples.")
  
  cat(sprintf("\nFinal samples: Control = %d, Severe = %d\n",
              sum(group == "Control"), sum(group == "Severe")))
  cat(sprintf("Genes: %d\n", nrow(expr_matrix)))
  
  if (any(is.na(expr_matrix))) {
    cat("Removing genes with NA values...\n")
    na_genes <- apply(expr_matrix, 1, function(x) any(is.na(x)))
    expr_matrix <- expr_matrix[!na_genes, ]
    cat(sprintf("Genes after removing NA: %d\n", nrow(expr_matrix)))
  }
  
  return(list(expr_data = expr_matrix, genes = rownames(expr_matrix),
              group = group))
}

process_GSE157103 <- function() {
  cat("\n========================================\n")
  cat("Processing validation dataset GSE157103\n")
  cat("========================================\n")
  
  gse_file <- "GSE157103_genes.tpm.tsv"
  if (!file.exists(gse_file)) stop("Validation file not found: ", gse_file)
  
  gse_data <- read.delim(gse_file, check.names = FALSE)
  genes <- if ("#symbol" %in% colnames(gse_data)) gse_data$`#symbol` else gse_data[, 1]
  expr_matrix <- as.matrix(gse_data[, -1])
  rownames(expr_matrix) <- genes
  
  sample_names <- colnames(expr_matrix)
  group <- ifelse(grepl("^C", sample_names), "COVID",
                  ifelse(grepl("^NC", sample_names), "Control", NA))
  keep <- !is.na(group)
  expr_matrix <- expr_matrix[, keep, drop = FALSE]
  group <- factor(group[keep], levels = c("Control", "COVID"))
  
  cat(sprintf("Samples: Control = %d, COVID = %d\n",
              sum(group == "Control"), sum(group == "COVID")))
  cat(sprintf("Genes: %d\n", nrow(expr_matrix)))
  
  keep_genes <- rowMeans(expr_matrix) > 0.1
  expr_matrix <- expr_matrix[keep_genes, ]
  cat(sprintf("Genes after filtering: %d\n", nrow(expr_matrix)))
  
  return(list(expr_data = expr_matrix, genes = rownames(expr_matrix),
              group = group))
}

# ====================================================
# 3. Differential expression
# ====================================================
perform_differential_expression <- function(expr_data, group) {
  cat("Performing differential expression analysis...\n")
  
  group <- droplevels(group)
  if (nlevels(group) != 2)
    stop("Expected exactly 2 group levels; got: ",
         paste(levels(group), collapse = ", "))
  tab <- table(group)
  cat("  - Group sizes:\n"); print(tab)
  if (any(tab == 0))
    stop("One group has zero samples: ",
         paste(names(tab)[tab == 0], collapse = ", "))
  if (any(is.na(group))) stop("NA values present in group vector")
  
  group_levels <- levels(group)
  control_group <- group_levels[1]
  case_group <- group_levels[2]
  cat(sprintf("  - Comparing %s vs %s\n", case_group, control_group))
  
  use_idx <- c(which(group == control_group), which(group == case_group))
  expr_subset <- expr_data[, use_idx, drop = FALSE]
  group_subset <- factor(group[use_idx],
                         levels = c(control_group, case_group))
  
  design <- model.matrix(~ 0 + group_subset)
  colnames(design) <- levels(group_subset)
  contrast.matrix <- makeContrasts(
    contrasts = paste0(case_group, " - ", control_group), levels = design)
  
  fit <- lmFit(expr_subset, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  degs <- topTable(fit2, number = Inf, adjust.method = "BH")
  degs$Gene <- rownames(degs)
  degs$DEG_strict <- ifelse(abs(degs$logFC) > 1.5 & degs$adj.P.Val < 0.01,
                            ifelse(degs$logFC > 0, "Up", "Down"), "No")
  degs$DEG_moderate <- ifelse(abs(degs$logFC) > 1 & degs$adj.P.Val < 0.05,
                              ifelse(degs$logFC > 0, "Up", "Down"), "No")
  degs$DEG_lenient <- ifelse(abs(degs$logFC) > 0.5 & degs$adj.P.Val < 0.1,
                             ifelse(degs$logFC > 0, "Up", "Down"), "No")
  
  write.csv(degs, file.path(output_dir, "DEG_results.csv"), row.names = FALSE)
  cat(sprintf("  Strict (|FC|>1.5, FDR<0.01): %d up, %d down\n",
              sum(degs$DEG_strict == "Up"), sum(degs$DEG_strict == "Down")))
  cat(sprintf("  Moderate (|FC|>1, FDR<0.05): %d up, %d down\n",
              sum(degs$DEG_moderate == "Up"), sum(degs$DEG_moderate == "Down")))
  cat(sprintf("  Lenient (|FC|>0.5, FDR<0.1): %d up, %d down\n",
              sum(degs$DEG_lenient == "Up"), sum(degs$DEG_lenient == "Down")))
  return(degs)
}

# ====================================================
# 4. Pathway enrichment
# ====================================================
perform_pathway_enrichment <- function(gene_list, background_genes = NULL) {
  cat("Performing pathway enrichment analysis...\n")
  if (length(gene_list) < 3) {
    warning("Too few genes for meaningful enrichment analysis")
    return(NULL)
  }
  suppressWarnings({
    entrez_list <- tryCatch(
      bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID",
           OrgDb = org.Hs.eg.db),
      error = function(e) { cat("  - Gene ID conversion error:",
                                e$message, "\n"); NULL })
  })
  if (is.null(entrez_list) || nrow(entrez_list) < 3) {
    cat("  - Insufficient genes converted to Entrez ID\n"); return(NULL)
  }
  cat(sprintf("  - Converted %d/%d genes to Entrez ID\n",
              nrow(entrez_list), length(gene_list)))
  
  pathway_results <- list()
  pathway_results$KEGG <- tryCatch({
    if (!is.null(background_genes)) {
      suppressWarnings({
        entrez_background <- bitr(background_genes, fromType = "SYMBOL",
                                  toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      })
      if (!is.null(entrez_background) && nrow(entrez_background) > 0) {
        enrichKEGG(gene = entrez_list$ENTREZID, organism = "hsa",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                   universe = entrez_background$ENTREZID)
      } else {
        enrichKEGG(gene = entrez_list$ENTREZID, organism = "hsa",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2)
      }
    } else {
      enrichKEGG(gene = entrez_list$ENTREZID, organism = "hsa",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    }
  }, error = function(e) { cat("  - KEGG enrichment failed:", e$message, "\n"); NULL })
  
  pathway_results$GO_BP <- tryCatch(
    enrichGO(gene = entrez_list$ENTREZID, OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
             pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
    error = function(e) { cat("  - GO enrichment failed:", e$message, "\n"); NULL })
  
  if (!is.null(pathway_results$KEGG) && nrow(pathway_results$KEGG) > 0)
    cat(sprintf("  - %d significant KEGG pathways\n", nrow(pathway_results$KEGG)))
  else cat("  - No significant KEGG pathways found\n")
  
  if (!is.null(pathway_results$GO_BP) && nrow(pathway_results$GO_BP) > 0)
    cat(sprintf("  - %d significant GO BP terms\n", nrow(pathway_results$GO_BP)))
  else cat("  - No significant GO terms found\n")
  return(pathway_results)
}

# ====================================================
# 5. ROC analysis
# ====================================================
perform_ROC_analysis <- function(expr_data, group, important_genes,
                                 dataset_name = "Main") {
  cat(sprintf("\nPerforming ROC analysis for %s dataset...\n", dataset_name))
  common_genes <- intersect(important_genes, rownames(expr_data))
  if (length(common_genes) < 2) {
    cat(sprintf("  - Only %d genes found, insufficient for multi-gene model\n",
                length(common_genes)))
    if (length(common_genes) == 0) return(NULL)
  }
  cat(sprintf("  - Using %d genes for ROC analysis\n", length(common_genes)))
  
  expr_subset <- expr_data[common_genes, , drop = FALSE]
  plot_data <- as.data.frame(t(expr_subset))
  plot_data$Group <- as.numeric(group) - 1
  if (length(unique(plot_data$Group)) < 2) {
    cat("  - Error: only one group present\n"); return(NULL)
  }
  
  results <- list()
  single_roc_list <- list()
  for (gene in common_genes) {
    if (sd(plot_data[[gene]]) > 0) {
      roc_obj <- tryCatch(roc(plot_data$Group, plot_data[[gene]], quiet = TRUE),
                          error = function(e) NULL)
      if (!is.null(roc_obj)) single_roc_list[[gene]] <- roc_obj
    } else {
      cat(sprintf("  - %s has zero variance, skipping\n", gene))
    }
  }
  
  if (length(common_genes) >= 3 && nrow(plot_data) >= 8) {
    set.seed(123)
    plot_data$Group_factor <- factor(
      ifelse(plot_data$Group == 1, "COVID", "Control"),
      levels = c("Control", "COVID"))
    min_class_n <- min(table(plot_data$Group_factor))
    cv_folds <- max(3, min(5, min_class_n))
    if (cv_folds >= 3) {
      train_control <- trainControl(method = "cv", number = cv_folds,
                                    classProbs = TRUE,
                                    summaryFunction = twoClassSummary)
      model_formula <- as.formula(paste("Group_factor ~",
                                        paste(common_genes, collapse = " + ")))
      model <- tryCatch(
        train(model_formula, data = plot_data, method = "glm",
              family = "binomial", trControl = train_control,
              metric = "ROC", na.action = na.omit),
        error = function(e) { cat("  - Logistic regression failed:",
                                  e$message, "\n"); NULL })
      if (!is.null(model)) {
        plot_data$pred_prob <- predict(model, plot_data, type = "prob")$COVID
        multi_roc <- tryCatch(roc(plot_data$Group, plot_data$pred_prob,
                                  quiet = TRUE), error = function(e) NULL)
        if (!is.null(multi_roc)) {
          results$multi_roc <- multi_roc
          results$model <- model
          results$cv_roc <- max(model$results$ROC, na.rm = TRUE)
          cat(sprintf("  - Cross-validated AUC: %.3f\n", results$cv_roc))
        }
      }
    }
  }
  
  if (length(single_roc_list) > 0) {
    gene_auc <- data.frame(
      Gene = names(single_roc_list),
      AUC = sapply(single_roc_list, function(x) as.numeric(auc(x))),
      stringsAsFactors = FALSE) %>% arrange(desc(AUC))
    results$gene_auc <- gene_auc
    results$single_roc <- single_roc_list
    cat("  - Top genes by AUC:\n")
    for (i in 1:min(3, nrow(gene_auc)))
      cat(sprintf("    %s: AUC = %.3f\n", gene_auc$Gene[i], gene_auc$AUC[i]))
  }
  return(results)
}

plot_ROC_curves <- function(roc_results, dataset_name = "Main") {
  if (is.null(roc_results) || length(roc_results$single_roc) == 0) {
    cat(sprintf("  - No ROC results to plot for %s\n", dataset_name)); return(NULL)
  }
  plot_data_list <- list()
  if (!is.null(roc_results$multi_roc)) {
    plot_data_list[[1]] <- data.frame(
      FPR = 1 - roc_results$multi_roc$specificities,
      TPR = roc_results$multi_roc$sensitivities,
      Model = "Multi-gene model",
      AUC = round(auc(roc_results$multi_roc), 3))
  }
  if (!is.null(roc_results$gene_auc) && nrow(roc_results$gene_auc) > 0) {
    top_n <- min(6, nrow(roc_results$gene_auc))
    for (gene in head(roc_results$gene_auc$Gene, top_n)) {
      if (!is.null(roc_results$single_roc[[gene]])) {
        roc_obj <- roc_results$single_roc[[gene]]
        plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
          FPR = 1 - roc_obj$specificities,
          TPR = roc_obj$sensitivities,
          Model = gene, AUC = round(auc(roc_obj), 3))
      }
    }
  }
  plot_df <- do.call(rbind, plot_data_list)
  plot_df$Label <- paste0(plot_df$Model, " (AUC = ",
                          sprintf("%.3f", plot_df$AUC), ")")
  n_models <- length(unique(plot_df$Model))
  colors <- brewer.pal(min(n_models, 8), "Set1")
  if (n_models > 8) colors <- c(colors, rainbow(n_models - 8))
  model_names <- unique(plot_df$Model)
  color_map <- setNames(colors[1:length(model_names)], model_names)
  if ("Multi-gene model" %in% model_names)
    color_map["Multi-gene model"] <- "#E41A1C"
  plot_df$LineSize <- ifelse(plot_df$Model == "Multi-gene model", 1.5, 1)
  
  p <- ggplot(plot_df, aes(x = FPR, y = TPR, color = Model, size = LineSize)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "grey50", size = 0.8) +
    scale_color_manual(values = color_map,
                       labels = setNames(unique(plot_df$Label),
                                         unique(plot_df$Model)), name = "") +
    scale_size_continuous(range = c(0.8, 1.5), guide = "none") +
    labs(title = ifelse(dataset_name == "Main", "A", "B"),
         x = "1 - Specificity (False Positive Rate)",
         y = "Sensitivity (True Positive Rate)") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 16),
          legend.position = c(0.65, 0.35),
          legend.background = element_rect(fill = "white", color = "grey80",
                                           size = 0.5),
          legend.text = element_text(size = 9),
          legend.key.size = unit(0.5, "cm")) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 1)))
  
  suffix <- ifelse(dataset_name == "Main", "4A", "4B")
  ggsave(file.path(output_dir, paste0("Figure_", suffix, ".pdf")),
         p, width = 9, height = 8)
  ggsave(file.path(output_dir, paste0("Figure_", suffix, ".png")),
         p, width = 9, height = 8, dpi = 300)
  cat(sprintf("  - Saved Figure %s\n", suffix))
  return(p)
}

plot_all_roc_curves <- function(roc_results, dataset_name = "Main",
                                output_suffix = "1") {
  cat(sprintf("\n--- Individual ROC curves (%s) ---\n", dataset_name))
  if (is.null(roc_results) || is.null(roc_results$single_roc)) {
    cat("  - No ROC results available\n"); return(NULL)
  }
  genes <- names(roc_results$single_roc)
  n_col <- 4; n_row <- ceiling(length(genes) / n_col)
  plot_list <- list()
  for (i in seq_along(genes)) {
    gene <- genes[i]
    roc_obj <- roc_results$single_roc[[gene]]
    auc_val <- round(auc(roc_obj), 3)
    rank <- which(roc_results$gene_auc$Gene == gene)
    df <- data.frame(FPR = 1 - roc_obj$specificities,
                     TPR = roc_obj$sensitivities)
    p <- ggplot(df, aes(x = FPR, y = TPR)) +
      geom_line(color = ifelse(rank <= 3, "#E41A1C", "#377EB8"), size = 1) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                  color = "grey50", size = 0.5) +
      labs(title = paste(gene, "(AUC =", sprintf("%.3f", auc_val), ")"),
           x = "1 - Specificity", y = "Sensitivity") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
            axis.title = element_text(size = 6),
            axis.text = element_text(size = 5),
            plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))
    plot_list[[i]] <- p
  }
  combined_plot <- wrap_plots(plot_list, ncol = n_col) +
    plot_annotation(
      title = paste("Supplementary Figure S", output_suffix, ": ",
                    dataset_name, " Dataset - Individual ROC Curves", sep = ""),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = 12)))
  height <- max(6, n_row * 2.5)
  ggsave(file.path(output_dir,
                   paste0("Supplementary_Figure_S", output_suffix, ".pdf")),
         combined_plot, width = 12, height = height)
  ggsave(file.path(output_dir,
                   paste0("Supplementary_Figure_S", output_suffix, ".png")),
         combined_plot, width = 12, height = height, dpi = 300)
  cat(sprintf("  - Saved Supplementary Figure S%s (%d genes)\n",
              output_suffix, length(genes)))
  return(combined_plot)
}

generate_auc_table <- function(roc_main, roc_valid, output_file) {
  cat("\n--- Generating AUC values table ---\n")
  main_auc <- roc_main$gene_auc; colnames(main_auc) <- c("Gene", "AUC_Main")
  valid_auc <- roc_valid$gene_auc
  colnames(valid_auc) <- c("Gene", "AUC_Validation")
  auc_table <- merge(main_auc, valid_auc, by = "Gene", all = TRUE)
  auc_table$CI_Main_Lower <- NA; auc_table$CI_Main_Upper <- NA
  auc_table$CI_Valid_Lower <- NA; auc_table$CI_Valid_Upper <- NA
  for (gene in roc_main$gene_auc$Gene) {
    if (!is.null(roc_main$single_roc[[gene]])) {
      ci_val <- ci.auc(roc_main$single_roc[[gene]])
      auc_table$CI_Main_Lower[auc_table$Gene == gene] <- round(ci_val[1], 3)
      auc_table$CI_Main_Upper[auc_table$Gene == gene] <- round(ci_val[3], 3)
    }
  }
  for (gene in roc_valid$gene_auc$Gene) {
    if (!is.null(roc_valid$single_roc[[gene]])) {
      ci_val <- ci.auc(roc_valid$single_roc[[gene]])
      auc_table$CI_Valid_Lower[auc_table$Gene == gene] <- round(ci_val[1], 3)
      auc_table$CI_Valid_Upper[auc_table$Gene == gene] <- round(ci_val[3], 3)
    }
  }
  auc_table$AUC_Main_CI <- sprintf("%.3f (%.3f-%.3f)", auc_table$AUC_Main,
                                   auc_table$CI_Main_Lower,
                                   auc_table$CI_Main_Upper)
  auc_table$AUC_Validation_CI <- sprintf("%.3f (%.3f-%.3f)",
                                         auc_table$AUC_Validation,
                                         auc_table$CI_Valid_Lower,
                                         auc_table$CI_Valid_Upper)
  auc_table <- auc_table[order(-auc_table$AUC_Main), ]
  output_table <- auc_table[, c("Gene", "AUC_Main_CI", "AUC_Validation_CI")]
  colnames(output_table) <- c("Gene", "AUC (Main Dataset)",
                              "AUC (Validation Dataset)")
  write.csv(output_table,
            file.path(output_dir, "10_tables", output_file), row.names = FALSE)
  cat(sprintf("  - Saved AUC table: %s\n", output_file))
  return(output_table)
}

# ====================================================
# 6. Network construction
# ====================================================
construct_network <- function(femt_db, overlap_genes, pathway_results,
                              deg_results) {
  cat("\nConstructing FEMT-PASC-CVS network...\n")
  cat("  - Input overlap_genes:", paste(overlap_genes, collapse = ", "), "\n")
  
  tones <- c("Zhi", "Yu", "Jue", "Shang", "Gong")
  tone_colors <- c(Zhi = "#E41A1C", Yu = "#377EB8", Jue = "#4DAF4A",
                   Shang = "#984EA3", Gong = "#FF7F00")
  tone_nodes <- data.frame(id = tones, label = tones, node_type = "Tone",
                           color = tone_colors[tones], size = 20,
                           stringsAsFactors = FALSE)
  
  # Map genes to tones dynamically from the FEMT database
  gene_tone_mapping <- setNames(femt_db$Tone, femt_db$Target)
  
  gene_entrez <- data.frame()
  for (gene in overlap_genes) {
    entrez <- tryCatch(
      mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID",
             keytype = "SYMBOL"),
      error = function(e) NA)
    gene_entrez <- rbind(gene_entrez, data.frame(
      symbol = gene, entrez = as.character(entrez), stringsAsFactors = FALSE))
  }
  
  gene_nodes <- data.frame()
  for (gene in overlap_genes) {
    tone <- ifelse(gene %in% names(gene_tone_mapping),
                   gene_tone_mapping[[gene]], "Unknown")
    logFC_val <- if (gene %in% deg_results$Gene)
      deg_results$logFC[deg_results$Gene == gene][1] else 0
    entrez_val <- gene_entrez$entrez[gene_entrez$symbol == gene][1]
    gene_nodes <- rbind(gene_nodes, data.frame(
      id = gene, label = gene, node_type = "Gene", tone = tone,
      logFC = logFC_val, entrez = entrez_val,
      size = 12 + 3 * abs(logFC_val), stringsAsFactors = FALSE))
  }
  cat("\n    Final gene nodes:\n")
  print(gene_nodes[, c("id", "tone", "logFC", "entrez")])
  
  pathway_nodes <- data.frame()
  if (!is.null(pathway_results$KEGG) &&
      !is.null(pathway_results$KEGG@result) &&
      nrow(pathway_results$KEGG@result) > 0) {
    top_pathways <- pathway_results$KEGG@result %>%
      arrange(pvalue) %>% head(5)
    for (i in 1:nrow(top_pathways)) {
      full_name <- top_pathways$Description[i]
      short_name <- case_when(
        grepl("T cell receptor", full_name) ~ "T cell receptor signaling",
        grepl("Th1 and Th2", full_name) ~ "Th1/Th2 differentiation",
        grepl("Th17", full_name) ~ "Th17 differentiation",
        grepl("Amphetamine", full_name) ~ "Amphetamine addiction",
        grepl("Tryptophan", full_name) ~ "Tryptophan metabolism",
        TRUE ~ substr(full_name, 1, 25))
      pathway_nodes <- rbind(pathway_nodes, data.frame(
        id = paste0("Pathway_", i), label = short_name,
        full_name = full_name, node_type = "Pathway", size = 12,
        stringsAsFactors = FALSE))
    }
  }
  
  nodes_list <- list(tone = tone_nodes[, c("id", "label", "node_type", "size")])
  if (nrow(gene_nodes) > 0)
    nodes_list$gene <- gene_nodes[, c("id", "label", "node_type", "size")]
  if (nrow(pathway_nodes) > 0)
    nodes_list$pathway <- pathway_nodes[, c("id", "label", "node_type", "size")]
  nodes <- do.call(rbind, nodes_list); rownames(nodes) <- NULL
  cat(sprintf("\n  Total nodes: %d\n", nrow(nodes)))
  print(table(nodes$node_type))
  
  edges <- data.frame()
  if (nrow(gene_nodes) > 0) {
    for (i in 1:nrow(gene_nodes)) {
      if (gene_nodes$tone[i] %in% tone_nodes$id)
        edges <- rbind(edges, data.frame(from = gene_nodes$tone[i],
                                         to = gene_nodes$id[i],
                                         stringsAsFactors = FALSE))
    }
  }
  if (nrow(gene_nodes) > 0 && nrow(pathway_nodes) > 0 &&
      !is.null(pathway_results$KEGG)) {
    kegg_result <- pathway_results$KEGG@result
    for (i in 1:nrow(pathway_nodes)) {
      if (i <= nrow(kegg_result)) {
        gene_ids <- kegg_result$geneID[i]
        if (!is.na(gene_ids) && gene_ids != "") {
          pathway_genes_entrez <- strsplit(gene_ids, "/")[[1]]
          for (j in 1:nrow(gene_nodes)) {
            if (!is.na(gene_nodes$entrez[j]) &&
                gene_nodes$entrez[j] %in% pathway_genes_entrez)
              edges <- rbind(edges, data.frame(
                from = gene_nodes$id[j], to = pathway_nodes$id[i],
                stringsAsFactors = FALSE))
          }
        }
      }
    }
  }
  if (nrow(edges) > 0) {
    edges$type <- ifelse(edges$from %in% tone_nodes$id,
                         "Tone-Gene", "Gene-Pathway")
    print(table(edges$type))
  }
  list(nodes = nodes, edges = edges, gene_details = gene_nodes)
}

generate_degree_centrality_table_manual <- function(
    network_result, femt_db,
    output_file = "Supplementary_Table_3_Degree_centrality.csv") {
  cat("\n========================================\n")
  cat("Generating Supplementary Table 3: Degree Centrality\n")
  cat("========================================\n")
  if (is.null(network_result) || nrow(network_result$edges) == 0 ||
      nrow(network_result$nodes) == 0) {
    cat("  - No network data available\n"); return(NULL)
  }
  nodes <- network_result$nodes
  edges <- network_result$edges
  degree_raw <- rep(0, nrow(nodes)); names(degree_raw) <- nodes$id
  for (i in 1:nrow(edges)) {
    if (edges$from[i] %in% names(degree_raw))
      degree_raw[edges$from[i]] <- degree_raw[edges$from[i]] + 1
    if (edges$to[i] %in% names(degree_raw))
      degree_raw[edges$to[i]] <- degree_raw[edges$to[i]] + 1
  }
  degree_raw <- as.numeric(degree_raw); names(degree_raw) <- nodes$id
  degree_norm <- degree_raw / (nrow(nodes) - 1)
  
  centrality_table <- data.frame(
    Node_ID = nodes$id,
    Node_Label = ifelse(is.na(nodes$label), nodes$id, nodes$label),
    Node_Type = nodes$node_type,
    Degree = degree_raw,
    Degree_Normalized = round(degree_norm, 4),
    stringsAsFactors = FALSE)
  centrality_table$Tone <- NA
  centrality_table$Organ <- NA
  centrality_table$Element <- NA
  for (i in 1:nrow(centrality_table)) {
    if (centrality_table$Node_Type[i] == "Gene") {
      db_entry <- femt_db[femt_db$Target == centrality_table$Node_Label[i], ]
      if (nrow(db_entry) > 0) {
        centrality_table$Tone[i] <- db_entry$Tone[1]
        centrality_table$Organ[i] <- db_entry$Organ[1]
        centrality_table$Element[i] <- db_entry$Element[1]
      }
    }
  }
  centrality_table <- centrality_table[order(-centrality_table$Degree), ]
  centrality_table$Rank <- 1:nrow(centrality_table)
  
  tone_nodes <- centrality_table[centrality_table$Node_Type == "Tone", ]
  gene_nodes <- centrality_table[centrality_table$Node_Type == "Gene", ]
  pathway_nodes <- centrality_table[centrality_table$Node_Type == "Pathway", ]
  
  if (nrow(tone_nodes) > 0) {
    tone_table <- tone_nodes[, c("Rank", "Node_Label", "Degree",
                                 "Degree_Normalized")]
    colnames(tone_table) <- c("Rank", "Tone", "Degree", "Degree_Normalized")
    write.csv(tone_table,
              file.path(output_dir, "10_tables", "Tone_centrality.csv"),
              row.names = FALSE)
  }
  if (nrow(gene_nodes) > 0) {
    gene_table <- gene_nodes[, c("Rank", "Node_Label", "Tone", "Organ",
                                 "Element", "Degree", "Degree_Normalized")]
    colnames(gene_table) <- c("Rank", "Gene", "Tone", "Organ", "Element",
                              "Degree", "Degree_Normalized")
    write.csv(gene_table,
              file.path(output_dir, "10_tables", "Gene_centrality.csv"),
              row.names = FALSE)
  }
  if (nrow(pathway_nodes) > 0) {
    pathway_table <- pathway_nodes[, c("Rank", "Node_Label", "Degree",
                                       "Degree_Normalized")]
    colnames(pathway_table) <- c("Rank", "Pathway", "Degree",
                                 "Degree_Normalized")
    write.csv(pathway_table,
              file.path(output_dir, "10_tables", "Pathway_centrality.csv"),
              row.names = FALSE)
  }
  write.csv(centrality_table,
            file.path(output_dir, "10_tables", output_file), row.names = FALSE)
  
  cat("\n  - Network summary:\n")
  cat(sprintf("    Total nodes: %d\n", nrow(nodes)))
  cat(sprintf("    Total edges: %d\n", nrow(edges)))
  cat(sprintf("    Network density: %.4f\n",
              2 * nrow(edges) / (nrow(nodes) * (nrow(nodes) - 1))))
  return(centrality_table)
}

plot_network <- function(network_result, filename = "FEMT_network") {
  if (is.null(network_result) || nrow(network_result$edges) == 0 ||
      nrow(network_result$nodes) == 0) {
    cat("  - No network to plot\n"); return(NULL)
  }
  cat("  - Creating Figure 5: FEMT-PASC-CVS interaction network...\n")
  g <- graph_from_data_frame(d = network_result$edges,
                             vertices = network_result$nodes,
                             directed = FALSE)
  node_colors <- c(Tone = "#E41A1C", Gene = "#377EB8", Pathway = "#4DAF4A")
  V(g)$color <- node_colors[V(g)$node_type]
  V(g)$color[is.na(V(g)$color)] <- "grey50"
  V(g)$size <- if ("size" %in% vertex_attr_names(g)) V(g)$size else 15
  V(g)$label <- if ("label" %in% vertex_attr_names(g)) V(g)$label else V(g)$name
  set.seed(123); l <- layout_with_fr(g)
  
  pdf(file.path(output_dir, "Figure_5.pdf"), width = 12, height = 8)
  plot(g, layout = l, vertex.color = V(g)$color, vertex.size = V(g)$size,
       vertex.label = V(g)$label, vertex.label.cex = 0.8,
       vertex.label.color = "black", vertex.label.family = "sans",
       edge.color = "grey80", edge.width = 1.5,
       main = "Figure 5. FEMT-PASC-CVS interaction network")
  legend("topright",
         legend = c("Tone (5 elements)", "Gene", "Pathway"),
         fill = c("#E41A1C", "#377EB8", "#4DAF4A"),
         title = "Node Type", cex = 0.9, bty = "n")
  dev.off()
  
  png(file.path(output_dir, "Figure_5.png"), width = 12, height = 8,
      units = "in", res = 300)
  plot(g, layout = l, vertex.color = V(g)$color, vertex.size = V(g)$size,
       vertex.label = V(g)$label, vertex.label.cex = 0.8,
       vertex.label.color = "black", vertex.label.family = "sans",
       edge.color = "grey80", edge.width = 1.5,
       main = "Figure 5. FEMT-PASC-CVS interaction network")
  legend("topright",
         legend = c("Tone (5 elements)", "Gene", "Pathway"),
         fill = c("#E41A1C", "#377EB8", "#4DAF4A"),
         title = "Node Type", cex = 0.9, bty = "n")
  dev.off()
  cat(sprintf("  - Network has %d nodes and %d edges\n", vcount(g), ecount(g)))
  return(g)
}

# ====================================================
# 7. Basic visualizations (Figures 1-3)
# ====================================================
generate_visualizations <- function(femt_db, deg_results, overlap_genes,
                                    pathway_results, geo_data = NULL) {
  cat("Generating visualizations for manuscript figures...\n")
  select <- dplyr::select
  
  # ---- Figure 1 ----
  cat("\n--- Figure 1 ---\n")
  p1A_data <- femt_db %>% group_by(Tone, Organ) %>%
    summarise(n_targets = n(), .groups = "drop")
  p1A <- ggplot(p1A_data, aes(x = Organ, y = n_targets, fill = Tone)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = n_targets), position = position_dodge(width = 0.9),
              vjust = -0.5, size = 3) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "A", x = "Organ System", y = "Number of Targets") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", size = 16),
          legend.position = "right")
  
  p1B_data <- femt_db %>% group_by(Evidence_Level) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(Percentage = round(Count / sum(Count) * 100, 1))
  evidence_levels <- c("Strong (animal/multi-omics)", "Moderate",
                       "Preliminary/Theoretical")
  evidence_colors <- c("#E41A1C", "#377EB8", "#4DAF4A")
  p1B <- ggplot(p1B_data, aes(x = "", y = Count,
                              fill = factor(Evidence_Level,
                                            levels = evidence_levels))) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    geom_text(aes(label = paste0(Count, " (", Percentage, "%)")),
              position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(values = evidence_colors, name = "Evidence Level") +
    labs(title = "B") + theme_void() +
    theme(plot.title = element_text(face = "bold", size = 16),
          legend.position = "bottom", legend.direction = "horizontal",
          legend.box = "horizontal", legend.text = element_text(size = 9),
          legend.title = element_text(size = 10, face = "bold"),
          legend.spacing.x = unit(0.5, "cm")) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  
  p1C_data <- femt_db %>% group_by(Pathway) %>%
    summarise(Count = n(), .groups = "drop") %>%
    arrange(desc(Count)) %>% head(10)
  p1C <- ggplot(p1C_data, aes(x = reorder(Pathway, Count), y = Count,
                              fill = Pathway)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), hjust = -0.1, size = 3) +
    scale_fill_brewer(palette = "Paired") + coord_flip() +
    labs(title = "C", x = "", y = "Number of Targets") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 16))
  
  figure1 <- (p1A | p1B) / p1C +
    plot_annotation(title = "Figure 1. FEMT target database characteristics",
                    theme = theme(plot.title = element_text(
                      hjust = 0.5, face = "bold", size = 14))) +
    plot_layout(heights = c(1, 0.8))
  ggsave(file.path(output_dir, "Figure_1.pdf"), figure1, width = 16, height = 12)
  ggsave(file.path(output_dir, "Figure_1.png"), figure1, width = 16,
         height = 12, dpi = 300)
  cat("  - Saved Figure 1\n")
  
  # ---- Figure 2 ----
  cat("\n--- Figure 2 ---\n")
  volcano_data <- deg_results
  volcano_data$Significance <- factor(
    ifelse(volcano_data$DEG_moderate != "No",
           ifelse(volcano_data$logFC > 0, "Up", "Down"),
           "Not Significant"),
    levels = c("Up", "Down", "Not Significant"))
  volcano_data$FEMT_Overlap <- ifelse(volcano_data$Gene %in% overlap_genes,
                                      "Yes", "No")
  p2A <- ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = Significance), size = 1.5, alpha = 0.6) +
    scale_color_manual(values = c(Up = "red", Down = "blue",
                                  `Not Significant` = "grey80")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
    labs(title = "A", x = "log2 Fold Change", y = "-log10(FDR)") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 16),
          legend.position = "bottom")
  overlap_data <- volcano_data[volcano_data$FEMT_Overlap == "Yes", ]
  if (nrow(overlap_data) > 0) {
    p2A <- p2A +
      geom_point(data = overlap_data,
                 aes(x = logFC, y = -log10(adj.P.Val)),
                 color = "black", size = 3, shape = 1) +
      geom_text_repel(data = overlap_data, aes(label = Gene), size = 3,
                      max.overlaps = 20)
  }
  venn_data <- list("FEMT Targets" = unique(femt_db$Target),
                    "PASC-CVS DEGs" = deg_results %>%
                      filter(DEG_moderate != "No") %>% pull(Gene))
  p2B <- ggvenn(venn_data, fill_color = c("#E41A1C", "#377EB8"),
                stroke_size = 0.5, set_name_size = 4, text_size = 3) +
    labs(title = "B") +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
  figure2 <- p2A + p2B +
    plot_annotation(title = "Figure 2. FEMT-PASC core gene identification",
                    theme = theme(plot.title = element_text(
                      hjust = 0.5, face = "bold", size = 14)))
  ggsave(file.path(output_dir, "Figure_2.pdf"), figure2, width = 14, height = 7)
  ggsave(file.path(output_dir, "Figure_2.png"), figure2, width = 14,
         height = 7, dpi = 300)
  cat("  - Saved Figure 2\n")
  
  # ---- Figure 3 ----
  cat("\n--- Figure 3 ---\n")
  if (!is.null(pathway_results$KEGG) &&
      nrow(pathway_results$KEGG@result) > 0) {
    top_kegg <- pathway_results$KEGG@result %>%
      arrange(pvalue) %>% head(10) %>%
      mutate(short_desc = str_wrap(Description, width = 40))
    p3A <- ggplot(top_kegg, aes(x = Count,
                                y = reorder(short_desc, -pvalue),
                                size = Count, color = -log10(pvalue))) +
      geom_point(alpha = 0.7) +
      scale_color_gradient(low = "blue", high = "red",
                           name = "-log10(p-value)",
                           guide = guide_colorbar(barwidth = 10,
                                                  barheight = 0.8)) +
      scale_size_continuous(range = c(3, 8), name = "Gene Count") +
      labs(title = "A", x = "Gene Count", y = "") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", size = 16),
            axis.text.y = element_text(size = 8),
            legend.position = "right", legend.box = "vertical")
  } else {
    p3A <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No significant KEGG") +
      theme_void() + labs(title = "A")
  }
  if (!is.null(pathway_results$GO_BP) &&
      nrow(pathway_results$GO_BP@result) > 0) {
    top_go <- pathway_results$GO_BP@result %>%
      arrange(pvalue) %>% head(10) %>%
      mutate(short_desc = str_wrap(Description, width = 40))
    p3B <- ggplot(top_go, aes(x = Count,
                              y = reorder(short_desc, -pvalue),
                              size = Count, color = -log10(pvalue))) +
      geom_point(alpha = 0.7) +
      scale_color_gradient(low = "blue", high = "red",
                           name = "-log10(p-value)",
                           guide = guide_colorbar(barwidth = 10,
                                                  barheight = 0.8)) +
      scale_size_continuous(range = c(3, 8), name = "Gene Count") +
      labs(title = "B", x = "Gene Count", y = "") +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold", size = 16),
            axis.text.y = element_text(size = 8),
            legend.position = "right", legend.box = "vertical")
  } else {
    p3B <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No significant GO") +
      theme_void() + labs(title = "B")
  }
  figure3 <- p3A + p3B +
    plot_annotation(title = "Figure 3. Functional enrichment of core genes",
                    theme = theme(plot.title = element_text(
                      hjust = 0.5, face = "bold", size = 14))) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal")
  ggsave(file.path(output_dir, "Figure_3.pdf"), figure3, width = 18, height = 8)
  ggsave(file.path(output_dir, "Figure_3.png"), figure3, width = 18,
         height = 8, dpi = 300)
  cat("  - Saved Figure 3\n")
}

# ====================================================
# 8. Figure 6: TF and drug-gene networks
# ====================================================
plot_tf_network <- function(overlap_genes) {
  cat("\n--- Figure 6A: TF regulatory network ---\n")
  tf_targets <- list(
    JUN    = c("JUN", "FOS", "EGR1", "EGR2", "EGR3", "IL10", "IL4"),
    REL    = c("IL10", "IL4", "IL2RB", "CD4", "CD8A"),
    STAT3  = c("IL10", "IL4", "IL2RB", "ACE", "ACE2"),
    STAT5  = c("IL2RB", "CD4", "CD8A", "IGF1"),
    NR3C1  = c("FKBP5", "NR3C1"),
    CREB1  = c("CREB1", "JUN", "FOS"),
    MAF    = c("IL10", "IL4")
  )
  edges <- data.frame(); tf_list <- c()
  for (tf in names(tf_targets)) {
    targets <- intersect(tf_targets[[tf]], overlap_genes)
    if (length(targets) > 0) {
      tf_list <- c(tf_list, tf)
      for (target in targets)
        edges <- rbind(edges, data.frame(from = tf, to = target,
                                         stringsAsFactors = FALSE))
    }
  }
  if (nrow(edges) == 0) { cat("  - No TF relationships found\n"); return(NULL) }
  all_nodes <- unique(c(edges$from, edges$to))
  nodes <- data.frame(id = all_nodes, label = all_nodes,
                      node_type = ifelse(all_nodes %in% tf_list, "TF", "Gene"),
                      stringsAsFactors = FALSE)
  g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
  V(g)$color <- ifelse(V(g)$node_type == "TF", "#E41A1C", "#377EB8")
  V(g)$size <- ifelse(V(g)$node_type == "TF", 20, 15)
  set.seed(123); l <- layout_with_fr(g)
  pdf(file.path(output_dir, "Figure_6A.pdf"), width = 12, height = 10)
  plot(g, layout = l, vertex.color = V(g)$color, vertex.size = V(g)$size,
       vertex.label = V(g)$label, vertex.label.cex = 0.8,
       vertex.label.color = "black", edge.color = "grey60", edge.width = 1,
       edge.arrow.size = 0.3,
       main = "Figure 6A. Transcription factor regulatory network")
  legend("topright", legend = c("Transcription Factor", "Target Gene"),
         fill = c("#E41A1C", "#377EB8"), title = "Node Type", cex = 0.8)
  dev.off()
  cat("  - Saved Figure 6A\n")
  return(g)
}

plot_drug_network <- function(overlap_genes) {
  cat("\n--- Figure 6B: Drug-gene interaction network ---\n")
  drug_targets <- list(
    ACE = c("Captopril", "Enalapril", "Lisinopril", "Ramipril"),
    ACE2 = c("DIZE", "Resorcinolnaphthalein"),
    CHRM3 = c("Atropine", "Scopolamine", "Tiotropium", "Ipratropium"),
    SLC12A3 = c("Hydrochlorothiazide", "Chlorthalidone"),
    IL10 = c("IL-10 recombinant protein", "Pegilodecakin"),
    IL4 = c("Dupilumab", "Lebrikizumab", "Tralokinumab"),
    IL2RB = c("Aldesleukin", "Denileukin diftitox"),
    CD4 = c("Ibalizumab", "Zanolimumab"),
    CD8A = c("T-cell therapies"),
    HTR2A = c("Risperidone", "Olanzapine", "Quetiapine", "Clozapine",
              "Nefazodone"),
    MAOA = c("Moclobemide", "Phenelzine", "Tranylcypromine"),
    MAOB = c("Selegiline", "Rasagiline", "Safinamide"),
    FKBP5 = c("SAFit1", "SAFit2", "Tacrolimus"),
    NR3C1 = c("Dexamethasone", "Prednisolone", "Hydrocortisone", "Cortisone"),
    IGF1 = c("Mecasermin", "Mecasermin rinfabate"),
    LEP = c("Metreleptin"),
    ADIPOQ = c("Adiponectin receptor agonists"),
    TMPRSS2 = c("Camostat mesylate", "Nafamostat mesylate"),
    NTRK2 = c("Amitriptyline", "Nortriptyline", "7,8-dihydroxyflavone"),
    CAMK2D = c("KN-93", "KN-62")
  )
  edges <- data.frame(); drug_list <- c()
  for (gene in overlap_genes) {
    if (gene %in% names(drug_targets)) {
      for (drug in drug_targets[[gene]]) {
        drug_name <- paste0("Drug_", gsub(" ", "_", drug))
        drug_list <- c(drug_list, drug_name)
        edges <- rbind(edges, data.frame(from = drug_name, to = gene,
                                         stringsAsFactors = FALSE))
      }
    }
  }
  if (nrow(edges) == 0) { cat("  - No drug interactions found\n"); return(NULL) }
  drug_nodes <- unique(edges$from); gene_nodes <- unique(edges$to)
  all_nodes <- data.frame(
    id = c(drug_nodes, gene_nodes),
    label = c(gsub("Drug_", "", drug_nodes), gene_nodes),
    node_type = c(rep("Drug", length(drug_nodes)),
                  rep("Gene", length(gene_nodes))),
    stringsAsFactors = FALSE)
  g <- graph_from_data_frame(edges, vertices = all_nodes, directed = TRUE)
  V(g)$color <- ifelse(V(g)$node_type == "Drug", "#4DAF4A", "#377EB8")
  V(g)$size <- ifelse(V(g)$node_type == "Drug", 15, 12)
  set.seed(123); l <- layout_with_fr(g)
  pdf(file.path(output_dir, "Figure_6B.pdf"), width = 16, height = 12)
  plot(g, layout = l, vertex.color = V(g)$color, vertex.size = V(g)$size,
       vertex.label = V(g)$label, vertex.label.cex = 0.6,
       vertex.label.color = "black", edge.color = "grey60", edge.width = 0.5,
       edge.arrow.size = 0.2,
       main = "Figure 6B. Drug-gene interaction network")
  legend("topright", legend = c("Drug", "Target Gene"),
         fill = c("#4DAF4A", "#377EB8"), title = "Node Type", cex = 0.8)
  dev.off()
  cat("  - Saved Figure 6B\n")
  return(g)
}

# ====================================================
# 9. Table 1 and Table 2 outputs
# ====================================================
output_table1 <- function(femt_db) {
  cat("\n--- Table 1: Distribution of FEMT targets by tone and organ ---\n")
  table1 <- femt_db %>%
    group_by(Tone, Element, Organ) %>%
    summarise(`Number of Targets` = n(), .groups = "drop") %>%
    arrange(match(Tone, c("Zhi", "Yu", "Jue", "Shang", "Gong")))
  write.csv(table1,
            file.path(output_dir, "10_tables", "Table_1_FEMT_distribution.csv"),
            row.names = FALSE)
  cat("  - Saved Table 1\n")
  print(table1)
  return(table1)
}

output_table2 <- function(femt_db, overlap_genes) {
  cat("\n--- Table 2: Integration of five-element theory with molecular findings ---\n")
  tone_order <- c("Zhi", "Yu", "Jue", "Shang", "Gong")
  element_map <- c(Zhi = "Fire", Yu = "Water", Jue = "Wood",
                   Shang = "Metal", Gong = "Earth")
  organ_map <- c(Zhi = "Heart", Yu = "Kidney", Jue = "Liver",
                 Shang = "Lung", Gong = "Spleen")
  rows <- list()
  for (tone in tone_order) {
    genes_in_tone <- intersect(overlap_genes,
                               femt_db$Target[femt_db$Tone == tone])
    rows[[tone]] <- data.frame(
      `Musical Tone` = tone,
      Element = element_map[[tone]],
      Organ = organ_map[[tone]],
      `Number of Genes` = length(genes_in_tone),
      `Key Genes` = ifelse(length(genes_in_tone) > 0,
                           paste(genes_in_tone, collapse = ", "), "None"),
      check.names = FALSE)
  }
  table2 <- do.call(rbind, rows); rownames(table2) <- NULL
  write.csv(table2,
            file.path(output_dir, "10_tables",
                      "Table_2_Five_element_integration.csv"),
            row.names = FALSE)
  cat("  - Saved Table 2\n")
  print(table2)
  return(table2)
}

# ====================================================
# 10. Main analysis
# ====================================================
main_analysis_enhanced <- function() {
  print_title("NETWORK PHARMACOLOGY ANALYSIS: FEMT FOR PASC-CVS")
  
  cat("STEP 1: Creating FEMT database...\n")
  femt_db <- create_enhanced_FEMT_database()
  
  cat("\nSTEP 2: Main dataset (GSE171110)...\n")
  geo_data_main <- process_GSE171110()
  cat(sprintf("  - Samples: %d\n", ncol(geo_data_main$expr_data)))
  cat(sprintf("  - Groups: %s\n",
              paste(levels(geo_data_main$group), collapse = ", ")))
  
  cat("\nSTEP 3: Validation dataset (GSE157103)...\n")
  geo_data_valid <- process_GSE157103()
  cat(sprintf("  - Samples: %d\n", ncol(geo_data_valid$expr_data)))
  cat(sprintf("  - Groups: %s\n",
              paste(levels(geo_data_valid$group), collapse = ", ")))
  
  cat("\nSTEP 4: Differential expression (main)...\n")
  deg_results_main <- perform_differential_expression(
    geo_data_main$expr_data, geo_data_main$group)
  
  cat("\nSTEP 5: Overlap analysis...\n")
  significant_degs <- deg_results_main %>%
    filter(DEG_moderate != "No") %>% pull(Gene)
  femt_targets <- unique(femt_db$Target)
  overlap_genes <- intersect(femt_targets, significant_degs)
  cat(sprintf("  - FEMT targets: %d\n", length(femt_targets)))
  cat(sprintf("  - Significant DEGs: %d\n", length(significant_degs)))
  cat(sprintf("  - Overlap genes: %d\n", length(overlap_genes)))
  if (length(overlap_genes) > 0)
    cat("  - Overlap genes:", paste(overlap_genes, collapse = ", "), "\n")
  
  cat("\nSTEP 6: Pathway enrichment...\n")
  pathway_results <- if (length(overlap_genes) >= 3)
    perform_pathway_enrichment(overlap_genes,
                               rownames(geo_data_main$expr_data)) else NULL
  
  cat("\nSTEP 7: ROC analysis...\n")
  roc_results_main <- NULL; roc_results_valid <- NULL
  if (length(overlap_genes) >= 2) {
    roc_results_main <- perform_ROC_analysis(
      geo_data_main$expr_data, geo_data_main$group, overlap_genes, "Main")
    valid_genes <- intersect(overlap_genes, rownames(geo_data_valid$expr_data))
    if (length(valid_genes) >= 2) {
      cat(sprintf("  - %d/%d overlap genes in validation dataset\n",
                  length(valid_genes), length(overlap_genes)))
      roc_results_valid <- perform_ROC_analysis(
        geo_data_valid$expr_data, geo_data_valid$group, valid_genes,
        "Validation")
    }
    if (!is.null(roc_results_main)) plot_ROC_curves(roc_results_main, "Main")
    if (!is.null(roc_results_valid)) plot_ROC_curves(roc_results_valid, "Validation")
    if (!is.null(roc_results_main) && !is.null(roc_results_valid))
      generate_auc_table(roc_results_main, roc_results_valid,
                         "Supplementary_Table_2_AUC_comparison.csv")
  }
  
  cat("\nSTEP 8: Network construction...\n")
  network_result <- construct_network(femt_db, overlap_genes,
                                      pathway_results, deg_results_main)
  if (!is.null(network_result)) plot_network(network_result)
  
  cat("\nSTEP 9: Basic visualizations (Figures 1-3)...\n")
  generate_visualizations(femt_db, deg_results_main, overlap_genes,
                          pathway_results, geo_data_main)
  
  cat("\nSTEP 10: Saving results...\n")
  save(femt_db, geo_data_main, geo_data_valid, deg_results_main,
       overlap_genes, pathway_results, roc_results_main, roc_results_valid,
       network_result,
       file = file.path(output_dir, "complete_analysis_results.RData"))
  
  print_line("=", 60)
  cat("ANALYSIS COMPLETED\n")
  print_line("=", 60)
  
  invisible(list(
    femt_database = femt_db, main_data = geo_data_main,
    validation_data = geo_data_valid, deg_results = deg_results_main,
    overlap_genes = overlap_genes, pathway_results = pathway_results,
    roc_main = roc_results_main, roc_valid = roc_results_valid,
    network = network_result))
}

results_enhanced <- main_analysis_enhanced()

# ====================================================
# 11. Post-analysis outputs (tables, Fig 6, S1/S2, S3)
# ====================================================
if (!is.null(results_enhanced)) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("POST-ANALYSIS OUTPUTS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  # Table 1: distribution of FEMT targets by tone and organ
  table1 <- output_table1(results_enhanced$femt_database)
  
  # Table 2: integration of five-element theory with molecular findings
  table2 <- output_table2(results_enhanced$femt_database,
                          results_enhanced$overlap_genes)
  
  # Supplementary Table S3: degree centrality
  if (!is.null(results_enhanced$network)) {
    generate_degree_centrality_table_manual(
      results_enhanced$network, results_enhanced$femt_database,
      "Supplementary_Table_3_Degree_centrality.csv")
  }
  
  # Figure 6A and 6B
  plot_tf_network(results_enhanced$overlap_genes)
  plot_drug_network(results_enhanced$overlap_genes)
  
  # Supplementary Figures S1 and S2: individual ROC curves
  if (!is.null(results_enhanced$roc_main))
    plot_all_roc_curves(results_enhanced$roc_main, "Main", "1")
  if (!is.null(results_enhanced$roc_valid))
    plot_all_roc_curves(results_enhanced$roc_valid, "Validation", "2")
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("ALL OUTPUTS COMPLETED\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
}

cat(rep("=", 70), "\n")
cat("Results directory:", output_dir, "\n")
cat("Key outputs:\n")
cat("  results/Figure_1.pdf                        FEMT target database\n")
cat("  results/Figure_2.pdf                        Volcano + Venn\n")
cat("  results/Figure_3.pdf                        KEGG/GO enrichment\n")
cat("  results/Figure_4A.pdf                       ROC - main dataset\n")
cat("  results/Figure_4B.pdf                       ROC - validation dataset\n")
cat("  results/Figure_5.pdf                        FEMT-PASC-CVS network\n")
cat("  results/Figure_6A.pdf                       TF regulatory network\n")
cat("  results/Figure_6B.pdf                       Drug-gene network\n")
cat("  results/Supplementary_Figure_S1.pdf         Individual ROC - main\n")
cat("  results/Supplementary_Figure_S2.pdf         Individual ROC - validation\n")
cat("  results/10_tables/Table_1_FEMT_distribution.csv\n")
cat("  results/10_tables/Table_2_Five_element_integration.csv\n")
cat("  results/10_tables/FEMT_database_complete.csv\n")
cat("  results/10_tables/Supplementary_Table_2_AUC_comparison.csv\n")
cat("  results/10_tables/Supplementary_Table_3_Degree_centrality.csv\n")