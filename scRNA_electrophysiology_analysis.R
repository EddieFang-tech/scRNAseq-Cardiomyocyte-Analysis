############################################################
# Final R Script for Single-cell RNA-seq Analysis
# Final four-sample dataset
# TN331-11 = hESC-CM
# TN331-12 = hESC-CM
# TN331-13 = Hypoimmunogenic hiPSC-CM
# TN331-14 = Hypoimmunogenic hiPSC-CM

############################################################

# TN331-11 = hESC-CM
# TN331-12 = hESC-CM
# TN331-13 = Hypoimmunogenic hiPSC-CM
# TN331-14 = Hypoimmunogenic hiPSC-CM
############################################################

############################################################
# 0. Load packages
############################################################

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(patchwork)
library(pheatmap)

set.seed(1234)

############################################################
# 1. Sample information
############################################################

sample_names <- c(
  "TN331-11",
  "TN331-12",
  "TN331-13",
  "TN331-14"
)

sample_group <- c(
  "TN331-11" = "hESC-CM",
  "TN331-12" = "hESC-CM",
  "TN331-13" = "Hypo_hiPSC-CM",
  "TN331-14" = "Hypo_hiPSC-CM"
)

cat("\nSample mapping:\n")
print(sample_group)

############################################################
# 2. Select the TOP-LEVEL scRNA folder
############################################################

select_root_dir <- function() {
  
  selected_dir <- NULL
  
  if (
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()
  ) {
    
    selected_dir <- tryCatch(
      rstudioapi::selectDirectory(
        caption = paste(
          "Select the TOP-LEVEL folder containing your scRNA data.",
          "The script will search recursively for TN331-11 to TN331-14."
        )
      ),
      error = function(e) NULL
    )
  }
  
  if (
    is.null(selected_dir) ||
    length(selected_dir) == 0 ||
    is.na(selected_dir) ||
    selected_dir == ""
  ) {
    
    if (.Platform$OS.type == "windows") {
      
      selected_dir <- tryCatch(
        utils::choose.dir(
          caption = "Select the TOP-LEVEL folder containing your scRNA data."
        ),
        error = function(e) NA_character_
      )
    }
  }
  
  if (
    is.null(selected_dir) ||
    length(selected_dir) == 0 ||
    is.na(selected_dir) ||
    selected_dir == ""
  ) {
    stop("No valid folder was selected.")
  }
  
  normalizePath(
    selected_dir,
    winslash = "/",
    mustWork = TRUE
  )
}

root_dir <- select_root_dir()

cat(
  "\n====================================\n",
  "Selected root directory:\n",
  root_dir,
  "\n====================================\n\n",
  sep = ""
)

############################################################
# 3. Output directory
############################################################

output_dir <- file.path(
  root_dir,
  "Fig3_final_4samples"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

cat(
  "Output directory:\n",
  output_dir,
  "\n\n",
  sep = ""
)

############################################################
# 4. Search recursively for 10x matrix belonging to sample
############################################################

find_sample_10x_dir <- function(root_dir, sample_id) {
  
  cat(
    "\n------------------------------------\n",
    "Searching for sample: ",
    sample_id,
    "\n------------------------------------\n",
    sep = ""
  )
  
  all_dirs <- list.dirs(
    root_dir,
    recursive = TRUE,
    full.names = TRUE
  )
  
  matrix_dirs <- all_dirs[
    basename(all_dirs) == "filtered_feature_bc_matrix"
  ]
  
  sample_matrix_dirs <- matrix_dirs[
    grepl(sample_id, matrix_dirs, fixed = TRUE)
  ]
  
  if (length(sample_matrix_dirs) == 1) {
    cat(
      "10x matrix found:\n",
      sample_matrix_dirs,
      "\n",
      sep = ""
    )
    return(sample_matrix_dirs)
  }
  
  if (length(sample_matrix_dirs) > 1) {
    cat(
      "\nMultiple matrix folders found for ",
      sample_id,
      ":\n",
      sep = ""
    )
    print(sample_matrix_dirs)
    cat(
      "\nUsing:\n",
      sample_matrix_dirs[1],
      "\n",
      sep = ""
    )
    return(sample_matrix_dirs[1])
  }
  
  sample_dirs <- all_dirs[
    grepl(sample_id, all_dirs, fixed = TRUE)
  ]
  
  if (length(sample_dirs) > 0) {
    for (sample_dir in sample_dirs) {
      
      candidates <- c(
        file.path(sample_dir, "filtered_feature_bc_matrix"),
        file.path(sample_dir, "outs", "filtered_feature_bc_matrix")
      )
      
      for (candidate in candidates) {
        if (dir.exists(candidate)) {
          cat(
            "\n10x matrix found:\n",
            candidate,
            "\n",
            sep = ""
          )
          return(candidate)
        }
      }
      
      files <- list.files(
        sample_dir,
        full.names = FALSE
      )
      
      has_matrix <- any(
        grepl("^matrix\\.mtx(\\.gz)?$", files)
      )
      
      has_barcodes <- any(
        grepl("^barcodes\\.tsv(\\.gz)?$", files)
      )
      
      has_features <- any(
        grepl("^(features|genes)\\.tsv(\\.gz)?$", files)
      )
      
      if (has_matrix && has_barcodes && has_features) {
        cat(
          "\n10x matrix files found directly in:\n",
          sample_dir,
          "\n",
          sep = ""
        )
        return(sample_dir)
      }
    }
  }
  
  all_files <- list.files(
    root_dir,
    recursive = TRUE,
    full.names = TRUE
  )
  
  matrix_files <- all_files[
    grepl("^matrix\\.mtx(\\.gz)?$", basename(all_files))
  ]
  
  matrix_files <- matrix_files[
    grepl(sample_id, matrix_files, fixed = TRUE)
  ]
  
  if (length(matrix_files) >= 1) {
    matrix_parent <- dirname(matrix_files[1])
    cat(
      "\nMatrix file found for ",
      sample_id,
      ":\n",
      matrix_parent,
      "\n",
      sep = ""
    )
    return(matrix_parent)
  }
  
  stop(
    paste0(
      "\nCould not find 10x data for ",
      sample_id,
      " under:\n",
      root_dir
    )
  )
}

############################################################
# 5. Locate all four datasets BEFORE analysis
############################################################

sample_matrix_paths <- setNames(
  rep(NA_character_, length(sample_names)),
  sample_names
)

for (s in sample_names) {
  sample_matrix_paths[s] <- find_sample_10x_dir(
    root_dir,
    s
  )
}

cat(
  "\n============================================\n",
  "FINAL 10x MATRIX PATHS\n",
  "============================================\n",
  sep = ""
)

print(sample_matrix_paths)

############################################################
# 6. Read individual samples
############################################################

seurat_list <- list()

for (s in sample_names) {
  
  message("\nReading sample: ", s)
  
  matrix_dir <- sample_matrix_paths[s]
  
  counts <- Read10X(
    data.dir = matrix_dir
  )
  
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }
  
  obj <- CreateSeuratObject(
    counts = counts,
    project = s,
    min.cells = 3,
    min.features = 200
  )
  
  obj$sample <- s
  obj$group <- unname(sample_group[s])
  
  obj[["percent.mt"]] <- PercentageFeatureSet(
    obj,
    pattern = "^MT-"
  )
  
  seurat_list[[s]] <- obj
  
  message(
    s,
    ": ",
    ncol(obj),
    " cells before QC"
  )
}

############################################################
# 7. QC summary before filtering
############################################################

qc_before <- bind_rows(
  lapply(
    names(seurat_list),
    function(s) {
      obj <- seurat_list[[s]]
      data.frame(
        sample = s,
        group = obj$group[1],
        cells_before_QC = ncol(obj),
        median_genes = median(obj$nFeature_RNA),
        median_UMI = median(obj$nCount_RNA),
        median_percent_mt = median(obj$percent.mt)
      )
    }
  )
)

write.csv(
  qc_before,
  file.path(output_dir, "QC_before_filtering.csv"),
  row.names = FALSE
)

print(qc_before)

############################################################
# 8. QC filtering
############################################################

# Inclusion criteria:
# nFeature_RNA >= 200
# nFeature_RNA <= 5000
# percent.mt <= 15

for (s in names(seurat_list)) {
  obj <- seurat_list[[s]]
  
  obj <- subset(
    obj,
    subset =
      nFeature_RNA >= 200 &
      nFeature_RNA <= 5000 &
      percent.mt <= 15
  )
  
  seurat_list[[s]] <- obj
}

############################################################
# 9. QC summary after filtering
############################################################

qc_after <- bind_rows(
  lapply(
    names(seurat_list),
    function(s) {
      obj <- seurat_list[[s]]
      data.frame(
        sample = s,
        group = obj$group[1],
        cells_after_QC = ncol(obj),
        median_genes = median(obj$nFeature_RNA),
        median_UMI = median(obj$nCount_RNA),
        median_percent_mt = median(obj$percent.mt)
      )
    }
  )
)

write.csv(
  qc_after,
  file.path(output_dir, "QC_after_filtering.csv"),
  row.names = FALSE
)

cat("\nQC after filtering:\n")
print(qc_after)

############################################################
# 10. Normalize each sample
############################################################

for (s in names(seurat_list)) {
  obj <- seurat_list[[s]]
  
  obj <- NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  
  seurat_list[[s]] <- obj
}

############################################################
# 11. Integration
############################################################

integration_features <- SelectIntegrationFeatures(
  object.list = seurat_list,
  nfeatures = 3000
)

anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  anchor.features = integration_features,
  dims = 1:30
)

obj.integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:30
)

############################################################
# 12. PCA / UMAP / clustering
############################################################

DefaultAssay(obj.integrated) <- "integrated"

obj.integrated <- ScaleData(
  obj.integrated,
  verbose = FALSE
)

obj.integrated <- RunPCA(
  obj.integrated,
  npcs = 30,
  verbose = FALSE
)

obj.integrated <- RunUMAP(
  obj.integrated,
  reduction = "pca",
  dims = 1:30,
  seed.use = 1234
)

obj.integrated <- FindNeighbors(
  obj.integrated,
  dims = 1:30
)

obj.integrated <- FindClusters(
  obj.integrated,
  resolution = 0.4
)

############################################################
# 13. Save integrated object
############################################################

saveRDS(
  obj.integrated,
  file.path(output_dir, "TN331_4samples_integrated.rds")
)

############################################################
# 14. UMAP plots
############################################################

p_sample <- DimPlot(
  obj.integrated,
  reduction = "umap",
  group.by = "sample"
) +
  theme_classic()

ggsave(
  filename = file.path(output_dir, "UMAP_by_sample.png"),
  plot = p_sample,
  width = 7,
  height = 6,
  dpi = 300
)

p_group <- DimPlot(
  obj.integrated,
  reduction = "umap",
  group.by = "group"
) +
  theme_classic()

ggsave(
  filename = file.path(output_dir, "UMAP_by_cell_source.png"),
  plot = p_group,
  width = 7,
  height = 6,
  dpi = 300
)

p_cluster <- DimPlot(
  obj.integrated,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  theme_classic()

ggsave(
  filename = file.path(output_dir, "UMAP_by_Seurat_cluster.png"),
  plot = p_cluster,
  width = 7,
  height = 6,
  dpi = 300
)

############################################################
# 15. Switch to RNA assay
############################################################

DefaultAssay(obj.integrated) <- "RNA"

############################################################
# 16. Identity markers
############################################################

identity_genes <- c(
  "TNNT2",
  "MYH7",
  "PDGFRA",
  "ENG",
  "SOX2"
)

identity_genes <- identity_genes[
  identity_genes %in% rownames(obj.integrated)
]

for (gene in identity_genes) {
  
  p <- FeaturePlot(
    obj.integrated,
    features = gene,
    reduction = "umap",
    order = TRUE
  ) +
    ggtitle(gene) +
    theme_classic()
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0("FeaturePlot_", gene, ".png")
    ),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

############################################################
# 17. HCN expression
############################################################

hcn_genes <- c(
  "HCN1",
  "HCN2",
  "HCN4"
)

hcn_genes <- hcn_genes[
  hcn_genes %in% rownames(obj.integrated)
]

for (gene in hcn_genes) {
  
  p <- FeaturePlot(
    obj.integrated,
    features = gene,
    reduction = "umap",
    order = TRUE
  ) +
    ggtitle(gene) +
    theme_classic()
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0("FeaturePlot_", gene, "_all_samples.png")
    ),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

############################################################
# 18. HCN split by cell source
############################################################

for (gene in hcn_genes) {
  
  p <- FeaturePlot(
    obj.integrated,
    features = gene,
    reduction = "umap",
    split.by = "group",
    order = TRUE
  )
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0("FeaturePlot_", gene, "_by_group.png")
    ),
    plot = p,
    width = 11,
    height = 5,
    dpi = 300
  )
}

############################################################
# 19. HCN-positive percentage
############################################################

calculate_positive_percentage <- function(seurat_obj, gene) {
  
  expression_data <- FetchData(
    seurat_obj,
    vars = gene,
    layer = "data"
  )
  
  df <- data.frame(
    sample = seurat_obj$sample,
    group = seurat_obj$group,
    expression = expression_data[, gene]
  )
  
  result <- df %>%
    group_by(sample, group) %>%
    summarise(
      total_cells = n(),
      positive_cells = sum(expression > 0),
      percent_positive = 100 * positive_cells / total_cells,
      .groups = "drop"
    )
  
  result$gene <- gene
  return(result)
}

HCN_positive_results <- bind_rows(
  lapply(
    hcn_genes,
    function(g) {
      calculate_positive_percentage(
        obj.integrated,
        g
      )
    }
  )
)

write.csv(
  HCN_positive_results,
  file.path(output_dir, "HCN_positive_percentage_by_sample.csv"),
  row.names = FALSE
)

HCN_group_summary <- HCN_positive_results %>%
  group_by(group, gene) %>%
  summarise(
    mean_percent = mean(percent_positive),
    SD = sd(percent_positive),
    n = n(),
    .groups = "drop"
  )

write.csv(
  HCN_group_summary,
  file.path(output_dir, "HCN_positive_percentage_by_group.csv"),
  row.names = FALSE
)

cat("\nHCN-positive results:\n")
print(HCN_positive_results)

cat("\nHCN-positive group summary:\n")
print(HCN_group_summary)

############################################################
# 20. Electrophysiological-state gene modules
############################################################

genes_nodal <- c(
  "HCN4",
  "SHOX2",
  "TBX3",
  "TBX18",
  "ISL1",
  "GATA6"
)

genes_HCN1 <- c(
  "HCN1",
  "PTPN4",
  "LRRC10",
  "TRIP8B"
)

genes_conduction <- c(
  "SCN5A",
  "SCN1B",
  "SCN3B",
  "GJA1"
)

genes_EC <- c(
  "CACNA1C",
  "CACNB2",
  "RYR2",
  "ATP2A2",
  "PLN",
  "CASQ2",
  "CAMK2D",
  "SLC8A1"
)

genes_IK1 <- c(
  "KCNJ2"
)

genes_IKr <- c(
  "KCNH2"
)

genes_IKs <- c(
  "KCNQ1"
)

# Verify before final repository release
genes_atrial <- c(
  "MYL7",
  "NPPA",
  "NR2F2",
  "KCNA5"
)

# Verify before final repository release
genes_ventricular <- c(
  "MYH7",
  "MYL2",
  "IRX4",
  "HAND1"
)

############################################################
# 21. Keep genes actually present
############################################################

keep_present <- function(genes, object) {
  genes[
    genes %in% rownames(object)
  ]
}

genes_nodal <- keep_present(genes_nodal, obj.integrated)
genes_HCN1 <- keep_present(genes_HCN1, obj.integrated)
genes_conduction <- keep_present(genes_conduction, obj.integrated)
genes_EC <- keep_present(genes_EC, obj.integrated)
genes_IK1 <- keep_present(genes_IK1, obj.integrated)
genes_IKr <- keep_present(genes_IKr, obj.integrated)
genes_IKs <- keep_present(genes_IKs, obj.integrated)
genes_atrial <- keep_present(genes_atrial, obj.integrated)
genes_ventricular <- keep_present(genes_ventricular, obj.integrated)

cat("\nNodal genes:\n")
print(genes_nodal)

cat("\nHCN1 genes:\n")
print(genes_HCN1)

cat("\nConduction genes:\n")
print(genes_conduction)

cat("\nE-C coupling genes:\n")
print(genes_EC)

cat("\nIK1 genes:\n")
print(genes_IK1)

cat("\nIKr genes:\n")
print(genes_IKr)

cat("\nIKs genes:\n")
print(genes_IKs)

cat("\nAtrial genes:\n")
print(genes_atrial)

cat("\nVentricular genes:\n")
print(genes_ventricular)

############################################################
# 22. AddModuleScore
############################################################

DefaultAssay(obj.integrated) <- "RNA"

obj.integrated <- AddModuleScore(
  obj.integrated,
  features = list(genes_nodal),
  name = "EP_Nodal_",
  seed = 1234
)

obj.integrated <- AddModuleScore(
  obj.integrated,
  features = list(genes_HCN1),
  name = "EP_HCN1_",
  seed = 1234
)

obj.integrated <- AddModuleScore(
  obj.integrated,
  features = list(genes_conduction),
  name = "EP_Conduction_",
  seed = 1234
)

obj.integrated <- AddModuleScore(
  obj.integrated,
  features = list(genes_EC),
  name = "EP_EC_",
  seed = 1234
)

obj.integrated <- AddModuleScore(
  obj.integrated,
  features = list(genes_IK1),
  name = "EP_IK1_",
  seed = 1234
)

obj.integrated <- AddModuleScore(
  obj.integrated,
  features = list(genes_IKr),
  name = "EP_IKr_",
  seed = 1234
)

obj.integrated <- AddModuleScore(
  obj.integrated,
  features = list(genes_IKs),
  name = "EP_IKs_",
  seed = 1234
)

obj.integrated <- AddModuleScore(
  obj.integrated,
  features = list(genes_atrial),
  name = "EP_Atrial_",
  seed = 1234
)

obj.integrated <- AddModuleScore(
  obj.integrated,
  features = list(genes_ventricular),
  name = "EP_Ventricular_",
  seed = 1234
)

############################################################
# 23. Module score column names
############################################################

score_columns <- c(
  "EP_Nodal_1",
  "EP_HCN1_1",
  "EP_Conduction_1",
  "EP_EC_1",
  "EP_IK1_1",
  "EP_IKr_1",
  "EP_IKs_1",
  "EP_Atrial_1",
  "EP_Ventricular_1"
)

score_labels <- c(
  "Nodal core",
  "HCN1-high",
  "Fast depolarization/conduction",
  "E-C coupling high",
  "IK1 mature (KCNJ2)",
  "IKr (KCNH2)",
  "IKs (KCNQ1)",
  "Atrial identity",
  "Ventricular identity"
)

missing_score_columns <- setdiff(
  score_columns,
  colnames(obj.integrated@meta.data)
)

if (length(missing_score_columns) > 0) {
  stop(
    paste(
      "These module-score columns are missing:",
      paste(missing_score_columns, collapse = ", "),
      "\nPlease rerun AddModuleScore first."
    )
  )
}

############################################################
# 24. Calculate 70th percentile thresholds
############################################################

thresholds <- sapply(
  score_columns,
  function(x) {
    values <- obj.integrated@meta.data[, x]
    quantile(
      values,
      probs = 0.70,
      na.rm = TRUE
    )
  }
)

cat("\n70th percentile thresholds:\n")
print(thresholds)

threshold_table <- data.frame(
  module = score_labels,
  column = score_columns,
  threshold = as.numeric(thresholds)
)

write.csv(
  threshold_table,
  file.path(output_dir, "Electrophysiological_module_thresholds.csv"),
  row.names = FALSE
)

############################################################
# 25. Assign electrophysiological state
############################################################

score_matrix <- as.matrix(
  obj.integrated@meta.data[
    ,
    score_columns,
    drop = FALSE
  ]
)

EP_state <- rep(
  "Unclassified",
  nrow(score_matrix)
)

for (i in seq_len(nrow(score_matrix))) {
  
  scores_i <- score_matrix[
    i,
    ,
    drop = TRUE
  ]
  
  passed <- scores_i > thresholds
  
  if (any(passed, na.rm = TRUE)) {
    
    eligible_scores <- scores_i
    eligible_scores[!passed] <- -Inf
    
    max_index <- which.max(eligible_scores)
    EP_state[i] <- score_labels[max_index]
  }
}

obj.integrated$EP_state <- EP_state

cat("\nElectrophysiological-state counts:\n")
print(table(obj.integrated$EP_state))

############################################################
# 26. Electrophysiological-state UMAP
############################################################

p_EP <- DimPlot(
  obj.integrated,
  reduction = "umap",
  group.by = "EP_state",
  shuffle = TRUE
) +
  theme_classic()

ggsave(
  filename = file.path(output_dir, "UMAP_electrophysiological_states.png"),
  plot = p_EP,
  width = 9,
  height = 7,
  dpi = 300
)

p_EP_group <- DimPlot(
  obj.integrated,
  reduction = "umap",
  group.by = "EP_state",
  split.by = "group",
  shuffle = TRUE
)

ggsave(
  filename = file.path(output_dir, "UMAP_electrophysiological_states_by_group.png"),
  plot = p_EP_group,
  width = 13,
  height = 6,
  dpi = 300
)

############################################################
# 27. EP-state composition
############################################################

EP_composition <- obj.integrated@meta.data %>%
  count(
    sample,
    group,
    EP_state
  ) %>%
  group_by(sample) %>%
  mutate(
    percentage = 100 * n / sum(n)
  ) %>%
  ungroup()

write.csv(
  EP_composition,
  file.path(output_dir, "Electrophysiological_state_composition.csv"),
  row.names = FALSE
)

############################################################
# 28. Export module scores
############################################################

module_score_output <- obj.integrated@meta.data[
  ,
  c(
    "sample",
    "group",
    score_columns,
    "EP_state"
  ),
  drop = FALSE
]

write.csv(
  module_score_output,
  file.path(output_dir, "Electrophysiological_module_scores.csv"),
  row.names = TRUE
)

############################################################
# 29. Heatmap genes
############################################################

heatmap_genes <- c(
  
  # Pacemaker / HCN
  "HCN1",
  "HCN2",
  "HCN4",
  
  # Fast depolarization / conduction
  "SCN5A",
  "SCN1B",
  "SCN3B",
  "GJA1",
  
  # Potassium channels
  "KCNJ2",
  "KCNH2",
  "KCNQ1",
  "KCNE1",
  
  # Calcium handling / E-C coupling
  "CACNA1C",
  "CACNB2",
  "RYR2",
  "ATP2A2",
  "PLN",
  "CASQ2",
  "CAMK2D",
  "SLC8A1",
  
  # Atrial identity
  "MYL7",
  "NPPA",
  "NR2F2",
  "KCNA5",
  
  # Ventricular identity
  "MYH7",
  "MYL2",
  "IRX4",
  "HAND1"
)

heatmap_genes <- heatmap_genes[
  heatmap_genes %in% rownames(obj.integrated)
]

cat("\nGenes included in heatmap:\n")
print(heatmap_genes)

############################################################
# 30. Average expression by sample
############################################################

DefaultAssay(obj.integrated) <- "RNA"

avg_result <- AverageExpression(
  obj.integrated,
  assays = "RNA",
  features = heatmap_genes,
  group.by = "sample",
  layer = "data",
  verbose = FALSE
)

avg_expr <- avg_result$RNA

cat("\nAverage-expression matrix dimensions:\n")
print(dim(avg_expr))

cat("\nSamples present in average-expression matrix:\n")
print(colnames(avg_expr))

############################################################
# 31. Force final sample order
############################################################

desired_sample_order <- c(
  "TN331-11",
  "TN331-12",
  "TN331-13",
  "TN331-14"
)

missing_samples <- setdiff(
  desired_sample_order,
  colnames(avg_expr)
)

if (length(missing_samples) > 0) {
  stop(
    paste(
      "Missing samples from AverageExpression:",
      paste(missing_samples, collapse = ", ")
    )
  )
}

avg_expr <- avg_expr[
  ,
  desired_sample_order,
  drop = FALSE
]

############################################################
# 32. Per-gene Z-score
############################################################

heatmap_z <- t(
  scale(
    t(avg_expr)
  )
)

finite_rows <- apply(
  heatmap_z,
  1,
  function(x) {
    all(is.finite(x))
  }
)

if (any(!finite_rows)) {
  cat(
    "\nGenes removed because of zero variance/non-finite Z-score:\n"
  )
  print(
    rownames(heatmap_z)[!finite_rows]
  )
}

heatmap_z <- heatmap_z[
  finite_rows,
  ,
  drop = FALSE
]

############################################################
# 33. Save heatmap matrices
############################################################

write.csv(
  avg_expr,
  file.path(output_dir, "AverageExpression_by_sample.csv"),
  row.names = TRUE
)

write.csv(
  heatmap_z,
  file.path(output_dir, "AverageExpression_Zscore_by_sample.csv"),
  row.names = TRUE
)

############################################################
# 34. Draw heatmap
############################################################

png(
  filename = file.path(output_dir, "Electrophysiology_heatmap_4samples.png"),
  width = 2200,
  height = 2200,
  res = 300
)

pheatmap(
  heatmap_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  border_color = NA,
  fontsize_row = 9,
  fontsize_col = 10
)

dev.off()

############################################################
# 35. Group-level average expression
############################################################

avg_group_result <- AverageExpression(
  obj.integrated,
  assays = "RNA",
  features = heatmap_genes,
  group.by = "group",
  layer = "data",
  verbose = FALSE
)

avg_expr_group <- avg_group_result$RNA

write.csv(
  avg_expr_group,
  file.path(output_dir, "AverageExpression_by_cell_source.csv"),
  row.names = TRUE
)

############################################################
# 36. Save final object and metadata
############################################################

saveRDS(
  obj.integrated,
  file.path(output_dir, "TN331_final_4samples_with_EP_states.rds")
)

write.csv(
  obj.integrated@meta.data,
  file.path(output_dir, "Final_cell_metadata.csv"),
  row.names = TRUE
)

############################################################
# 37. Record software versions
############################################################

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

############################################################
# 38. Final report
############################################################

cat("\n============================================\n")
cat("ANALYSIS COMPLETED\n")
cat("============================================\n\n")

cat("Samples used:\n")
print(sample_names)

cat("\nSample-to-group mapping:\n")
print(sample_group)

cat("\n10x matrix paths:\n")
print(sample_matrix_paths)

cat("\nNumber of cells per sample after QC:\n")
print(table(obj.integrated$sample))

cat("\nNumber of cells per cell source:\n")
print(table(obj.integrated$group))

cat("\nElectrophysiological states:\n")
print(table(obj.integrated$EP_state))

cat("\nHCN-positive percentage:\n")
print(HCN_positive_results)

cat("\nOutput directory:\n")
print(output_dir)
