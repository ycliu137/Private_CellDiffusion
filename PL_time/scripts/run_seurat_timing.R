#!/usr/bin/env Rscript
#
# Seurat integration with timing measurements
# Records dataset statistics and running time for each step
#

library(Seurat)
library(Matrix)
library(jsonlite)

# Get arguments from Snakemake
input_h5ad <- snakemake@input$h5ad
output_h5ad <- snakemake@output$h5ad
output_timing <- snakemake@output$timing
output_stats <- snakemake@output$stats
params <- snakemake@params

batch_key <- params$batch_key
label_key <- params$label_key
normalized_data <- params$normalized_data

# Initialize timing and stats lists
# Extract dataset name from input path
dataset_name <- basename(dirname(input_h5ad))

timing_list <- list(
    dataset = dataset_name,
    method = "Seurat",
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    steps = list()
)

stats_list <- list(
    dataset = dataset_name,
    method = "Seurat",
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
)

cat("\n=== Seurat Timing Pipeline ===\n")
cat("Input:", input_h5ad, "\n")
cat("Output:", output_h5ad, "\n")

# ===== Load data =====
t0 <- Sys.time()
cat("\n=== Loading data ===\n")

# Read h5ad using Seurat's Read10X_h5 alternative (use rhdf5 for .h5ad)
library(rhdf5)
h5f <- H5Fopen(input_h5ad, flags="H5F_ACC_RDONLY")

# Read expression matrix
X <- h5read(h5f, "X")
if (class(X)[1] == "dgCMatrix" || class(X)[1] == "matrix") {
    expr_matrix <- X
} else {
    expr_matrix <- t(X)
}

# Read observation metadata
obs_names <- h5read(h5f, "obs_names")
var_names <- h5read(h5f, "var_names")
obs_keys <- h5ls(h5f, recursive=FALSE)$name[grep("^obs$", h5ls(h5f, recursive=FALSE)$name)]
obs_data <- list()
for (key in obs_keys) {
    tryCatch({
        obs_data[[key]] <- h5read(h5f, paste0("obs/", key))
    }, error = function(e) {})
}

H5Fclose(h5f)

colnames(expr_matrix) <- obs_names
rownames(expr_matrix) <- var_names

cat("Data shape:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "cells\n")

# Record statistics
stats_list$n_cells <- ncol(expr_matrix)
stats_list$n_genes <- nrow(expr_matrix)

# Get batch info if available
if (batch_key %in% names(obs_data)) {
    batch_col <- obs_data[[batch_key]]
    stats_list$n_batches <- length(unique(batch_col))
    cat("Batches:", stats_list$n_batches, "\n")
} else {
    stats_list$n_batches <- 1
    cat("No batch key found, assuming single batch\n")
}

timing_list$steps$load_data <- as.numeric(difftime(Sys.time(), t0, units="secs"))

# Create Seurat object
t0 <- Sys.time()
seurat_obj <- CreateSeuratObject(counts=expr_matrix, min.cells=3, min.features=200)

if (batch_key %in% names(obs_data)) {
    seurat_obj[[batch_key]] <- obs_data[[batch_key]]
}
if (label_key %in% names(obs_data)) {
    seurat_obj[[label_key]] <- obs_data[[label_key]]
}

timing_list$steps$create_seurat <- as.numeric(difftime(Sys.time(), t0, units="secs"))

# ===== Step 0: Preprocessing =====
cat("\n=== Step 0: Preprocessing ===\n")
t0 <- Sys.time()

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures=params$nfeatures)
seurat_obj <- ScaleData(seurat_obj)

timing_list$steps$preprocessing <- as.numeric(difftime(Sys.time(), t0, units="secs"))
cat("Preprocessing time:", timing_list$steps$preprocessing, "s\n")

# ===== Step 1: PCA =====
cat("\n=== Step 1: PCA ===\n")
t0 <- Sys.time()

seurat_obj <- RunPCA(seurat_obj, npcs=params$dims, verbose=FALSE)

timing_list$steps$pca <- as.numeric(difftime(Sys.time(), t0, units="secs"))
cat("PCA time:", timing_list$steps$pca, "s\n")

# ===== Step 2: Batch Correction (Harmony integration) =====
cat("\n=== Step 2: Harmony Batch Correction ===\n")
t0 <- Sys.time()

library(harmony)
seurat_obj <- RunHarmony(
    seurat_obj,
    group.by.vars=batch_key,
    dims.use=1:params$dims,
    verbose=FALSE
)

timing_list$steps$harmony <- as.numeric(difftime(Sys.time(), t0, units="secs"))
cat("Harmony time:", timing_list$steps$harmony, "s\n")

# ===== Step 3: UMAP =====
cat("\n=== Step 3: UMAP ===\n")
t0 <- Sys.time()

seurat_obj <- RunUMAP(
    seurat_obj,
    reduction="harmony",
    dims=1:params$dims,
    verbose=FALSE
)
seurat_obj[["umap_seurat"]] <- seurat_obj[["umap"]]

timing_list$steps$umap <- as.numeric(difftime(Sys.time(), t0, units="secs"))
cat("UMAP time:", timing_list$steps$umap, "s\n")

# ===== Step 4: Leiden Clustering =====
cat("\n=== Step 4: Leiden Clustering ===\n")
cat("  Resolution:", params$resolution, "\n")
t0 <- Sys.time()

seurat_obj <- FindNeighbors(seurat_obj, dims=1:params$dims, verbose=FALSE)
seurat_obj <- FindClusters(
    seurat_obj,
    resolution=params$resolution,
    algorithm=4,  # Leiden algorithm
    verbose=FALSE
)
seurat_obj$leiden_seurat <- seurat_obj$seurat_clusters

timing_list$steps$leiden <- as.numeric(difftime(Sys.time(), t0, units="secs"))
cat("Leiden time:", timing_list$steps$leiden, "s\n")

# ===== Save results =====
cat("\n=== Saving results ===\n")

dir.create(dirname(output_h5ad), showWarnings=FALSE, recursive=TRUE)

# Save as h5ad using Seurat's SaveH5Seurat
SaveH5Seurat(seurat_obj, filename=sub("\\.h5ad$", "", output_h5ad), overwrite=TRUE)

# Alternatively, save as RDS
rds_path <- sub("\\.h5ad$", ".rds", output_h5ad)
saveRDS(seurat_obj, file=rds_path)
cat("Seurat object saved to:", rds_path, "\n")

# Calculate total time
total_time <- sum(unlist(timing_list$steps))
timing_list$total_time <- total_time
cat("\nTotal Seurat time:", total_time, "s\n")

# Save timing results as JSON
dir.create(dirname(output_timing), showWarnings=FALSE, recursive=TRUE)
write_json(timing_list, output_timing, pretty=TRUE)
cat("Timing saved to:", output_timing, "\n")

# Save statistics as JSON
dir.create(dirname(output_stats), showWarnings=FALSE, recursive=TRUE)
write_json(stats_list, output_stats, pretty=TRUE)
cat("Statistics saved to:", output_stats, "\n")

cat("\n=== Seurat timing pipeline complete! ===\n")
