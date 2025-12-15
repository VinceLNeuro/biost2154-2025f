# Load VLMC subset from H5AD and prepare for trajectory analysis
library(tidyverse)
library(Seurat)
library(anndata) #do not use SeuratDisk
library(monocle3)
library(patchwork)

# Path to VLMC subset
h5ad_path <- "/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/data/MTG/SEAAD_MTG_RNAseq_final-nuclei_VLMC_subset.h5ad"
output_dir <- "/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/code/5_trajectoryAnalysis/out_vlmc_trajectory"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading H5AD file using anndata...\n")
adata <- read_h5ad(h5ad_path)

# Extract count matrix
counts <- t(as.matrix(adata$X))
rownames(counts) <- rownames(adata$var)
colnames(counts) <- rownames(adata$obs)

# Create Seurat object
vlmc_seurat <- CreateSeuratObject(
  counts = counts,
  meta.data = as.data.frame(adata$obs),
  project = "VLMC_SEA-AD"
)
vlmc_seurat

# Summary
cat("\nDonor summary:\n")
print(table(vlmc_seurat$`Donor ID`))

cat("\nCognitive Status summary:\n")
print(table(vlmc_seurat$`Cognitive Status`))

# Save R object for downstream analysis
saveRDS(vlmc_seurat, "/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/data/MTG/vlmc_seurat.rds")

