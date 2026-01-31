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

# List all available datasets/groups
all_keys <- h5ls(h5f)
cat("Available keys in h5ad:\n")
print(all_keys)

# h5ad stores sparse matrices in CSR/CSC format with data, indices, indptr
expr_matrix <- NULL

# Check if X exists and its structure
X_list <- h5ls(h5f)
X_entries <- X_list[X_list$group == "/X", ]

cat("\nFound X, structure:\n")
print(X_entries)

# Check if it's a sparse matrix (has data, indices, indptr)
if ("data" %in% X_entries$name && "indices" %in% X_entries$name && "indptr" %in% X_entries$name) {
    cat("Reading sparse matrix in CSR format\n")
    data <- h5read(h5f, "X/data")
    indices <- h5read(h5f, "X/indices")
    indptr <- h5read(h5f, "X/indptr")
    
    cat("Data length:", length(data), "\n")
    cat("Indices length:", length(indices), "\n")
    cat("Indptr length:", length(indptr), "\n")
    
    # In CSR format:
    # - indptr describes row structure, so n_rows = length(indptr) - 1
    # - indices are column indices, so n_cols = max(indices) + 1
    n_rows <- length(indptr) - 1
    n_cols <- max(indices, na.rm=TRUE) + 1
    
    cat("Inferred sparse matrix shape: ", n_rows, "rows (cells) x", n_cols, "cols (genes)\n")
    
    # Verify dimensions match
    cat("Max index value:", max(indices, na.rm=TRUE), "\n")
    cat("Expected max index:", n_cols - 1, "\n")
    
    # Convert 0-based indices to 1-based for R
    indices <- indices + 1
    
    # Create sparse matrix from CSR format
    library(Matrix)
    # CSR format: indptr describes where each row starts
    # We need to create (row_indices, col_indices, values) for sparse matrix
    row_indices <- rep(1:n_rows, diff(indptr))
    col_indices <- indices
    values <- data
    
    cat("Creating sparse matrix with dimensions:", n_rows, "x", n_cols, "\n")
    cat("Row indices range: ", min(row_indices), "-", max(row_indices), "\n")
    cat("Col indices range: ", min(col_indices), "-", max(col_indices), "\n")
    
    expr_matrix <- sparseMatrix(
        i = row_indices,
        j = col_indices,
        x = values,
        dims = c(n_rows, n_cols),
        giveCsparse = FALSE
    )
} else {
    cat("X is stored as dense matrix or other format\n")
    cat("X entries:", paste(X_entries$name, collapse=", "), "\n")
    if ("data" %in% X_entries$name) {
        cat("Reading from X/data\n")
        expr_matrix <- h5read(h5f, "X/data")
    } else {
        stop("Cannot determine X format")
    }
}

if (is.null(expr_matrix)) {
    stop("Could not read expression matrix from X")
}

cat("Expression matrix class:", class(expr_matrix), "\n")
cat("Expression matrix dimensions:", dim(expr_matrix), "\n")

# Read observation and variable names from h5ad structure
obs_names <- as.character(h5read(h5f, "obs/_index"))
var_names <- as.character(h5read(h5f, "var/_index"))

cat("obs_names length:", length(obs_names), "\n")
cat("var_names length:", length(var_names), "\n")

# Read observation metadata (batch information)
obs_data <- list()
tryCatch({
    obs_contents <- h5ls(h5f, path="/obs")
    obs_keys <- obs_contents$name[obs_contents$name != "_index"]
    cat("Found obs keys:", paste(obs_keys, collapse=", "), "\n")
    
    for (key in obs_keys) {
        tryCatch({
            obs_data[[key]] <- h5read(h5f, paste0("obs/", key))
        }, error = function(e) {
            cat("Warning: Could not read obs/", key, "\n")
        })
    }
}, error = function(e) {
    cat("Warning: Could not read obs metadata\n")
})

H5Fclose(h5f)

# Set matrix names
if (ncol(expr_matrix) == length(obs_names) && nrow(expr_matrix) == length(var_names)) {
    colnames(expr_matrix) <- obs_names
    rownames(expr_matrix) <- var_names
} else if (nrow(expr_matrix) == length(obs_names) && ncol(expr_matrix) == length(var_names)) {
    cat("Transposing matrix to match dimensions\n")
    expr_matrix <- t(expr_matrix)
    colnames(expr_matrix) <- obs_names
    rownames(expr_matrix) <- var_names
} else {
    cat("ERROR: Dimension mismatch!\n")
    cat("Matrix:", nrow(expr_matrix), "x", ncol(expr_matrix), "\n")
    cat("obs_names:", length(obs_names), "var_names:", length(var_names), "\n")
    stop("Cannot match matrix dimensions with index names")
}

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
