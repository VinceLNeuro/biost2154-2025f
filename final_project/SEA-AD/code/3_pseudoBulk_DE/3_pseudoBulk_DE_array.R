options(echo=TRUE)

# Get array task ID from command line
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])

cat(sprintf("\n========== ARRAY JOB %d ==========\n", task_id))

######## Load libraries ########
library(tidyverse)
library(skimr)
library(BiocParallel)
# Use fewer cores for array jobs
n_cores <- 6
register(MulticoreParam(workers = n_cores))
library(DESeq2)
library(apeglm)

################ Load Clean & Formatted Data ################
data = readRDS("/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/data/MTG/pseudoBulk_20251109/cleaned_formatted_pseudobulk_data.rds")

skim(data$coldata_sub)
# Need to scale numeric variables (due to large mean; indicated by DESeq2)
data$coldata_sub <- data$coldata_sub %>%
  mutate(
    Age_at_Death_scaled = scale(Age_at_Death)[,1],
    PMI_scaled = scale(PMI)[,1],
    RIN_scaled = scale(RIN)[,1],
    Years_of_education_scaled = scale(Years_of_education)[,1]
  )
skim(data$coldata_sub)

# Get unique cell types
cell_types <- unique(data$coldata_sub$Subclass)
cell_types <- sort(cell_types)

cat(sprintf("Total cell types: %d\n", length(cell_types)))
cat(sprintf("Processing cell type %d: %s\n\n", task_id, cell_types[task_id]))

# Select the cell type for this array job
ct <- cell_types[task_id]

######## Analysis Parameters ########
comparison_col = "latenc_st2"
base_formula = "~ latenc_st2 + RIN_scaled + PMI_scaled + Age_at_Death_scaled + Sex + Years_of_education_scaled + APOE_Genotype + Cognitive_Status"
adjusted_formula = "~ latenc_st2 + RIN_scaled + PMI_scaled + Age_at_Death_scaled + Sex + Years_of_education_scaled + APOE_Genotype + Cognitive_Status + Overall_AD_neuropath + Overall_CAA_Score + Highest_Lewy_Body + Total_Microinfarcts + Atherosclerosis + Arteriolosclerosis"
min_counts = 10
fdr_threshold = 0.1
lfc_threshold = 0.25
shrinkage_method = "apeglm"

######## Process Single Cell Type ########
message(paste0("========== Analyzing '", ct, "' =========="))

# Subset metadata for this cell type
meta_subset <- data$coldata_sub %>% as.data.frame() %>%
  filter(Subclass == ct)

# Get corresponding column names in the matrix
sample_cols <- meta_subset$Donor_Subclass

# Extract counts for this cell type
counts <- data$counts_sub[, sample_cols, drop = FALSE]

# Set rownames to match column names
rownames(meta_subset) <- meta_subset$Donor_Subclass
message(sprintf("  N samples: %d", ncol(counts)))

# Check sample distribution
group_counts <- table(meta_subset[[comparison_col]])
message(sprintf("  Sample distribution:"))
for (grp in names(group_counts)) {
  message(sprintf("    %s: %d", grp, group_counts[grp]))
}

#### Filter low-count genes ####
smallestGroupSize <- min(table(meta_subset[[comparison_col]]))
keep <- rowSums(counts >= min_counts) >= smallestGroupSize
counts_filtered <- counts[keep, , drop = FALSE]
message(sprintf("  Retained %d/%d genes after filtering", 
                sum(keep), nrow(counts)))

#### BASE MODEL ####
message("\n=== Running BASE model ===")
res_base_df <- tibble()
dds_base <- NULL

tryCatch({
  dds_base <- DESeqDataSetFromMatrix(
    countData = counts_filtered,
    colData = meta_subset,
    design = as.formula(base_formula)
  )
  
  # Run DESeq2
  dds_base <- DESeq(dds_base, parallel=TRUE)
  
  # Get coefficient name
  coef_name <- resultsNames(dds_base)
  coef_name <- coef_name[grep(comparison_col, coef_name)]
  
  if (length(coef_name) == 0) {
    stop("Could not find coefficient for comparison")
  }
  
  # Extract results
  res_base <- results(dds_base, name = coef_name, 
                      independentFiltering = TRUE,
                      pAdjustMethod = "BH", alpha = fdr_threshold)
  
  # Shrink log2 fold changes
  res_base_shrink <- lfcShrink(
    dds_base,
    coef = coef_name,
    res = res_base,
    type = shrinkage_method,
    quiet = TRUE
  )
  
  summary(res_base_shrink)
  
  res_base_df <- as.data.frame(res_base_shrink) %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    arrange(padj)
  
  message("  Base model: SUCCESS")
  
}, error = function(e) {
  message(sprintf("  Base model: ERROR - %s", e$message))
})


#### ADJUSTED MODEL ####
message("\n=== Running ADJUSTED model ===")
res_adj_df <- tibble()
dds_adj <- NULL

tryCatch({
  dds_adj <- DESeqDataSetFromMatrix(
    countData = counts_filtered,
    colData = meta_subset,
    design = as.formula(adjusted_formula)
  )
  
  # Run DESeq2
  dds_adj <- DESeq(dds_adj, parallel=TRUE)
  
  # Get coefficient name
  coef_name <- resultsNames(dds_adj)
  coef_name <- coef_name[grep(comparison_col, coef_name)]
  
  if (length(coef_name) == 0) {
    stop("Could not find coefficient for comparison")
  }
  
  # Extract results
  res_adj <- results(dds_adj, name = coef_name, 
                     independentFiltering = TRUE,
                     pAdjustMethod = "BH", alpha = fdr_threshold)
  
  # Shrink log2 fold changes
  res_adj_shrink <- lfcShrink(
    dds_adj,
    coef = coef_name,
    res = res_adj,
    type = shrinkage_method,
    quiet = TRUE
  )
  
  summary(res_adj_shrink)
  
  res_adj_df <- as.data.frame(res_adj_shrink) %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    arrange(padj)
  
  message("  Adjusted model: SUCCESS")
  
}, error = function(e) {
  message(sprintf("  Adjusted model: ERROR - %s", e$message))
})

#### Filter for DEGs ####
base_degs <- res_base_df %>%
  filter(!is.na(padj) & 
           padj < fdr_threshold & 
           abs(log2FoldChange) > lfc_threshold)

adjusted_degs <- res_adj_df %>%
  filter(!is.na(padj) &
           padj < fdr_threshold &
           abs(log2FoldChange) > lfc_threshold)

message(sprintf("\n  Base model DEGs: %d", nrow(base_degs)))
message(sprintf("  Adjusted model DEGs: %d", nrow(adjusted_degs)))

#### Save Results ####
output_dir <- "/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/code/3_pseudoBulk_DE/results_by_celltype"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Clean cell type name for filename
ct_clean <- gsub("/", "_", ct)
ct_clean <- gsub(" ", "_", ct_clean)

result_obj <- list(
  cell_type = ct,
  dds_base = dds_base,
  dds_adjusted = dds_adj,
  base_model = res_base_df,
  adjusted_model = res_adj_df,
  base_degs = base_degs,
  adjusted_degs = adjusted_degs,
  metadata = meta_subset,
  n_samples = ncol(counts_filtered),
  parameters = list(
    comparison_col = comparison_col,
    base_formula = base_formula,
    adjusted_formula = adjusted_formula,
    min_counts = min_counts,
    fdr_threshold = fdr_threshold,
    lfc_threshold = lfc_threshold,
    shrinkage_method = shrinkage_method
  )
)

output_file <- file.path(output_dir, sprintf("%s_results.rds", ct_clean))
saveRDS(result_obj, output_file)
message(sprintf("\n=== Results saved to: %s ===\n", output_file))

cat("\n========== ARRAY JOB COMPLETED ==========\n")
