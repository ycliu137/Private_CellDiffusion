#!/usr/bin/env Rscript
#
# Seurat integration with timing measurements
# Records dataset statistics and running time for each step
#

library(Seurat)
library(Matrix)
library(jsonlite)
library(rhdf5)

# ========================================
# Helper function: Read h5ad into Seurat
# ========================================
read_h5ad_to_seurat <- function(h5ad_path) {
    # Robustly read h5ad file and convert to Seurat object
    # Handles:
    # - CSC sparse matrix format (AnnData default)
    # - Correct obs/var dimension mapping
    # - Batch and label metadata
    
    h5f <- H5Fopen(h5ad_path, flags="H5F_ACC_RDONLY")
    
    # Cache h5ls() result for performance
    h5_tree <- h5ls(h5f)
    
    # ===== Read X matrix (expression data) =====
    cat("Reading expression matrix X...\n")
    
    X_entries <- h5_tree[h5_tree$group == "/X", ]
    
    # h5ad stores X as CSC (Column Sparse Compressed)
    # NOT CSR (Row Sparse Compressed)
    # Key point: indptr describes COLUMNS, not rows
    
    if ("data" %in% X_entries$name && 
        "indices" %in% X_entries$name && 
        "indptr" %in% X_entries$name) {
        
        cat("  Found CSC sparse matrix\n")
        data <- as.numeric(as.vector(h5read(h5f, "X/data")))
        indices <- as.integer(as.vector(h5read(h5f, "X/indices")))
        indptr <- as.integer(as.vector(h5read(h5f, "X/indptr")))
        
        # CSC format interpretation:
        # indptr: column pointers (length = n_cols + 1)
        # indices: row indices (0-based, length = nnz)
        # data: values (length = nnz)
        
        n_cols <- length(indptr) - 1  # Number of columns (genes in CSC)
        n_rows <- max(indices, na.rm=TRUE) + 1  # Number of rows (cells in CSC)
        
        cat(sprintf("  CSC matrix dims: %d rows (cells) × %d cols (genes)\n", n_rows, n_cols))
        
        # Build sparse matrix from CSC format
        # For sparseMatrix(), we need (i, j) = (row, col) indices
        col_indices <- rep(seq_len(n_cols), diff(indptr))  # Column index for each element
        row_indices <- indices + 1L  # Row index (convert 0-based to 1-based)
        values <- data
        
        expr_matrix <- sparseMatrix(
            i = row_indices,
            j = col_indices,
            x = values,
            dims = c(n_rows, n_cols)
        )
        
    } else {
        stop(sprintf("Cannot find CSC sparse matrix components in %s", h5ad_path))
    }
    
    cat(sprintf("  Expression matrix loaded: %d genes × %d cells\n", 
                nrow(expr_matrix), ncol(expr_matrix)))
    
    # ===== Read obs and var indices =====
    cat("Reading metadata...\n")
    
    obs_names <- as.character(h5read(h5f, "obs/_index"))
    var_names <- as.character(h5read(h5f, "var/_index"))
    
    cat(sprintf("  obs: %d cells\n", length(obs_names)))
    cat(sprintf("  var: %d genes\n", length(var_names)))
    
    # ===== Dimension sanity check =====
    # h5ad convention: X is cells × genes (CSC format)
    # This means nrow(X) = cells, ncol(X) = genes
    
    if (nrow(expr_matrix) != length(obs_names)) {
        stop(sprintf(
            "Dimension mismatch: expr_matrix has %d rows but obs_names has %d\n",
            nrow(expr_matrix), length(obs_names)
        ))
    }
    
    if (ncol(expr_matrix) != length(var_names)) {
        stop(sprintf(
            "Dimension mismatch: expr_matrix has %d cols but var_names has %d\n",
            ncol(expr_matrix), length(var_names)
        ))
    }
    
    # ===== Set dimension names =====
    # Currently: expr_matrix is cells × genes (AnnData convention)
    rownames(expr_matrix) <- obs_names
    colnames(expr_matrix) <- var_names
    
    # Seurat expects: genes × cells (transpose required)
    expr_matrix <- t(expr_matrix)
    
    # ===== Read obs metadata =====
    obs_data <- list()
    tryCatch({
        # Build full paths for checking (h5ls()$name only returns basename)
        all_paths <- paste(h5_tree$group, h5_tree$name, sep="/")
        
        obs_meta_entries <- h5_tree[h5_tree$group == "/obs" & h5_tree$otype == "H5I_GROUP", ]
        
        for (i in seq_len(nrow(obs_meta_entries))) {
            key <- obs_meta_entries$name[i]
            if (key == "_index") next  # Skip index
            
            tryCatch({
                # Try to read codes if it exists (for categorical)
                codes_path <- sprintf("/obs/%s/codes", key)
                cats_path <- sprintf("/obs/%s/categories", key)
                
                if (codes_path %in% all_paths && cats_path %in% all_paths) {
                    # Categorical variable
                    codes <- h5read(h5f, codes_path)
                    categories <- h5read(h5f, cats_path)
                    obs_data[[key]] <- as.factor(categories[codes + 1])
                } else {
                    # Try to read as regular array (e.g., /obs/batch/data)
                    data_path <- sprintf("/obs/%s/data", key)
                    if (data_path %in% all_paths) {
                        obs_data[[key]] <- h5read(h5f, data_path)
                    }
                }
            }, error = function(e) {
                cat(sprintf("  Warning: Could not read obs/%s: %s\n", key, e$message))
            })
        }
    }, error = function(e) {
        cat("  Warning: Could not read obs metadata\n")
    })
    
    H5Fclose(h5f)
    
    # ===== Create Seurat object =====
    # Seurat expects: genes × cells
    seurat_obj <- CreateSeuratObject(
        counts = expr_matrix,
        min.cells = 3,
        min.features = 200
    )
    
    # Add metadata
    for (key in names(obs_data)) {
        if (length(obs_data[[key]]) == ncol(seurat_obj)) {
            seurat_obj[[key]] <- obs_data[[key]]
        }
    }
    
    cat(sprintf("Seurat object created: %d genes × %d cells\n",
                nrow(seurat_obj), ncol(seurat_obj)))
    
    # Return dataset stats along with object
    return(list(
        obj = seurat_obj,
        n_cells = ncol(seurat_obj),
        n_genes = nrow(seurat_obj),
        obs_data = obs_data
    ))
}

# ========================================
# Main script
# ========================================

# Get arguments from Snakemake
input_h5ad <- snakemake@input$h5ad
output_h5ad <- snakemake@output$h5ad  # Will be .h5seurat
output_timing <- snakemake@output$timing
output_stats <- snakemake@output$stats
params <- snakemake@params

batch_key <- params$batch_key
label_key <- params$label_key

# Initialize timing and stats
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

h5ad_data <- read_h5ad_to_seurat(input_h5ad)
seurat_obj <- h5ad_data$obj
obs_data <- h5ad_data$obs_data

stats_list$n_cells <- h5ad_data$n_cells
stats_list$n_genes <- h5ad_data$n_genes

# Determine number of batches
if (batch_key %in% names(obs_data)) {
    batch_col <- obs_data[[batch_key]]
    stats_list$n_batches <- length(unique(batch_col))
    cat("Batches:", stats_list$n_batches, "\n")
} else {
    stats_list$n_batches <- 1
    cat("No batch key found, assuming single batch\n")
}

timing_list$steps$load_data <- as.numeric(difftime(Sys.time(), t0, units="secs"))

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

# ===== Step 2: Batch Correction (Harmony) - GUARDED =====
cat("\n=== Step 2: Batch Correction ===\n")

# Only run Harmony if we have multiple batches
n_batches <- if (batch_key %in% names(obs_data)) {
    length(unique(obs_data[[batch_key]]))
} else {
    1
}

if (n_batches > 1 && batch_key %in% names(obs_data)) {
    t0 <- Sys.time()
    cat("Running Harmony with", n_batches, "batches...\n")
    
    library(harmony)
    seurat_obj <- RunHarmony(
        seurat_obj,
        group.by.vars=batch_key,
        dims.use=1:params$dims,
        verbose=FALSE
    )
    
    timing_list$steps$harmony <- as.numeric(difftime(Sys.time(), t0, units="secs"))
    cat("Harmony time:", timing_list$steps$harmony, "s\n")
    use_reduction <- "harmony"
} else {
    cat("Skipping Harmony: single batch or no batch key\n")
    timing_list$steps$harmony <- 0
    use_reduction <- "pca"
}

# ===== Step 3: UMAP =====
cat("\n=== Step 3: UMAP ===\n")
t0 <- Sys.time()

seurat_obj <- RunUMAP(
    seurat_obj,
    reduction=use_reduction,
    verbose=FALSE
)

timing_list$steps$umap <- as.numeric(difftime(Sys.time(), t0, units="secs"))
cat("UMAP time:", timing_list$steps$umap, "s\n")

# ===== Step 4: Leiden Clustering =====
cat("\n=== Step 4: Leiden Clustering ===\n")
t0 <- Sys.time()

seurat_obj <- FindNeighbors(seurat_obj, reduction=use_reduction, verbose=FALSE)
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

# Save as RDS (more stable than h5seurat)
# RDS is native R format, doesn't require SeuratDisk
cat("Saving Seurat object to:", output_h5ad, "\n")
saveRDS(seurat_obj, file=output_h5ad)

# Calculate total time with na.rm=TRUE to handle skipped steps
timing_list$total_time <- sum(unlist(timing_list$steps), na.rm=TRUE)
cat("\nTotal Seurat time:", timing_list$total_time, "s\n")

# Save timing and stats as JSON
dir.create(dirname(output_timing), showWarnings=FALSE, recursive=TRUE)
write_json(timing_list, output_timing, pretty=TRUE)
cat("Timing saved to:", output_timing, "\n")

dir.create(dirname(output_stats), showWarnings=FALSE, recursive=TRUE)
write_json(stats_list, output_stats, pretty=TRUE)
cat("Statistics saved to:", output_stats, "\n")

cat("\n=== Seurat timing pipeline complete! ===\n")

