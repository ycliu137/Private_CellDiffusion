#!/usr/bin/env Rscript
#
# Seurat integration with timing measurements
# Records dataset statistics and running time for each step
#

library(Seurat)
library(Matrix)
library(jsonlite)

# ========================================
# Helper function: Read 10X data into Seurat
# ========================================
read_10x_to_seurat <- function(data_dir) {
    # Read 10X data directory and associated metadata CSV files
    # Handles:
    # - 10X matrix.mtx, genes.tsv, barcodes.tsv
    # - batch.csv (1 column, no header)
    # - labels.csv (1 column, no header)
    
    cat("Reading 10X data from:", data_dir, "\n")
    
    # ===== Read 10X matrix =====
    seurat_obj <- Read10X(data.dir = data_dir)
    
    # If Read10X returns a list (multiple feature types), take the first
    if (is.list(seurat_obj)) {
        cat("  Multiple feature types detected, using first one\n")
        seurat_obj <- seurat_obj[[1]]
    }
    
    cat(sprintf("  10X matrix loaded: %d genes × %d cells\n", 
                nrow(seurat_obj), ncol(seurat_obj)))
    
    # ===== Create Seurat object =====
    seurat_obj <- CreateSeuratObject(
        counts = seurat_obj,
        min.cells = 3,
        min.features = 200
    )
    
    cat(sprintf("Seurat object created: %d genes × %d cells\n",
                nrow(seurat_obj), ncol(seurat_obj)))
    
    # ===== Read metadata CSVs =====
    obs_data <- list()
    
    # Read batch.csv
    batch_file <- file.path(data_dir, "batch.csv")
    if (file.exists(batch_file)) {
        cat("Reading batch metadata from:", batch_file, "\n")
        batch_data <- read.csv(batch_file, header=FALSE, stringsAsFactors=FALSE)[[1]]
        
        if (length(batch_data) == ncol(seurat_obj)) {
            seurat_obj$batch <- batch_data
            obs_data$batch <- batch_data
            cat("  Batch metadata added:", length(unique(batch_data)), "batches\n")
        } else {
            warning(sprintf("Batch file length (%d) does not match cell count (%d)",
                          length(batch_data), ncol(seurat_obj)))
        }
    } else {
        warning("batch.csv not found in", data_dir)
    }
    
    # Read labels.csv
    labels_file <- file.path(data_dir, "labels.csv")
    if (file.exists(labels_file)) {
        cat("Reading label metadata from:", labels_file, "\n")
        label_data <- read.csv(labels_file, header=FALSE, stringsAsFactors=FALSE)[[1]]
        
        if (length(label_data) == ncol(seurat_obj)) {
            seurat_obj$label <- label_data
            obs_data$label <- label_data
            cat("  Label metadata added:", length(unique(label_data)), "unique labels\n")
        } else {
            warning(sprintf("Labels file length (%d) does not match cell count (%d)",
                          length(label_data), ncol(seurat_obj)))
        }
    } else {
        warning("labels.csv not found in", data_dir)
    }
    
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
data_dir <- snakemake@input$data_dir
output_h5ad <- snakemake@output$h5ad  # Will be .rds
output_timing <- snakemake@output$timing
output_stats <- snakemake@output$stats
params <- snakemake@params

batch_key <- params$batch_key
label_key <- params$label_key

# Initialize timing and stats
dataset_name <- basename(data_dir)

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
cat("Input:", data_dir, "\n")
cat("Output:", output_h5ad, "\n")

# ===== Load data =====
t0 <- Sys.time()
cat("\n=== Loading data ===\n")

data_10x <- read_10x_to_seurat(data_dir)
seurat_obj <- data_10x$obj
obs_data <- data_10x$obs_data

stats_list$n_cells <- data_10x$n_cells
stats_list$n_genes <- data_10x$n_genes

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

