######## Load libraries ########
rm(list=ls())
library(tidyverse)

######## Combine Array Job Results ########
results_dir <- "/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/code/3_pseudoBulk_DE/results_by_celltype"

# List all result files
result_files <- list.files(results_dir, pattern = "_results.rds$", full.names = TRUE)
cat(sprintf("Found %d result files\n", length(result_files)))

# Load all results
all_results <- list()
for (f in result_files) {
  ct_name <- basename(f) %>% str_remove("_results.rds")
  cat(sprintf("Loading: %s\n", ct_name))
  all_results[[ct_name]] <- readRDS(f)
}

######## Create Summary Tables ########

# Summary of DEG counts
deg_summary <- map_df(names(all_results), function(ct) {
  res <- all_results[[ct]]
  tibble(
    cell_type = res$cell_type,
    n_samples = res$n_samples,
    base_degs = nrow(res$base_degs),
    adjusted_degs = nrow(res$adjusted_degs)
  )
}) %>%
  arrange(desc(adjusted_degs))

print(deg_summary)

# # Save combined results
# output_file <- "/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/code/3_pseudoBulk_DE/de_results_lateNC_combined.rds"
# saveRDS(all_results, output_file)
# cat(sprintf("\nCombined results saved to: %s\n", output_file))
# 
# # Save summary table
# summary_file <- "/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/code/3_pseudoBulk_DE/deg_summary.csv"
# write_csv(deg_summary, summary_file)
# cat(sprintf("Summary table saved to: %s\n", summary_file))

######## Create Summary Barplot ########
# Filter out cell types with 0 DEGs in both models
deg_summary_filtered <- deg_summary %>%
  filter(base_degs > 0 | adjusted_degs > 0)

deg_summary_long <- deg_summary_filtered %>%
  pivot_longer(cols = c(base_degs, adjusted_degs),
               names_to = "model",
               values_to = "n_degs") %>%
  mutate(model = factor(model, 
                       levels = c("base_degs", "adjusted_degs"),
                       labels = c("Base Model", "Adjusted Model")))

p <- ggplot(deg_summary_long, aes(x = reorder(cell_type, n_degs), y = n_degs, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n_degs), 
            position = position_dodge(width = 0.9), hjust = -0.2, size = 2) +
  coord_flip() +
  labs(title = "Number of DEGs by Cell Type (LATE-NC)",
       subtitle = paste0("FDR < 0.1, |log2FC| > 0.25"),
       x = "Cell Type",
       y = "Number of DEGs",
       fill = "Model") +
  scale_fill_manual(values = c("Base Model" = "steelblue", 
                               "Adjusted Model" = "darkblue")) +
  scale_y_sqrt() +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p)

# ggsave("final_project/SEA-AD/code/3_pseudoBulk_DE/deg_counts_barplot.png",
#        plot = p, width = 5, height = 6, dpi = 500)


######## Combine all DEG results from both base and adjusted models ########
# Base model results
all_degs_base <- map_df(names(all_results), function(ct) {
  res_df <- all_results[[ct]]$base_degs
  res_df %>%
    mutate(cell_type = ct, model = "base", .before = 1)
}) %>%
  arrange(cell_type, padj)

# Adjusted model results
all_degs_adjusted <- map_df(names(all_results), function(ct) {
  res_df <- all_results[[ct]]$adjusted_degs
  res_df %>%
    mutate(cell_type = ct, model = "adjusted", .before = 1)
}) %>%
  arrange(cell_type, padj)

# Combine both
all_degs <- bind_rows(all_degs_base, all_degs_adjusted)

# # Save complete DEG tables
# degs_file <- "/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/code/3_pseudoBulk_DE/all_degs_base_and_adjusted.csv"
# write_csv(all_degs, degs_file)
# cat(sprintf("All DEGs (base + adjusted) saved to: %s\n", degs_file))

all_degs %>% filter(padj < 0.05) %>% arrange(padj)


######## FDR < 0.05 barplot ########
deg_summary_005 <- map_df(names(all_results), function(ct) {
  res <- all_results[[ct]]
  tibble( #regain the summary information for alpha = 0.05
    cell_type = res$cell_type,
    base_degs = res$base_model  %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.25) %>% nrow(),
    adjusted_degs = res$adjusted_model %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.25) %>% nrow()
  )
}) %>%
  filter(base_degs > 0 | adjusted_degs > 0) %>%
  pivot_longer(cols = c(base_degs, adjusted_degs),
               names_to = "model", values_to = "n_degs") %>%
  mutate(model = factor(model,
                        levels = c("base_degs", "adjusted_degs"),
                        labels = c("Base Model", "Adjusted Model")))

p_005 <- ggplot(deg_summary_005,
                aes(x = reorder(cell_type, n_degs), y = n_degs, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n_degs),
            position = position_dodge(width = 0.9), hjust = -0.2, size = 3) +
  coord_flip() +
  labs(title = "Number of DEGs by Cell Type (LATE-NC)",
       subtitle = "FDR < 0.05, |log2FC| > 0.25",
       x = "Cell Type", y = "Number of DEGs", fill = "Model") +
  scale_fill_manual(values = c("Base Model" = "steelblue", "Adjusted Model" = "darkblue")) +
  # scale_y_sqrt() +
  theme_minimal() +
  theme(legend.position = "bottom")
p_005

# ggsave("/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/code/3_pseudoBulk_DE/deg_counts_barplot_fdr0.05.png",
#        plot = p_005, width = 5, height = 6, dpi = 500)


######## Check L2/3 IT specifically (non-convergence) ########
library(DESeq2)
dds_l23 <- all_results[["L2_3_IT"]]$dds_adjusted
# png(filename = "final_project/SEA-AD/code/3_pseudoBulk_DE/weird_res_L2_3_IT/plot_dispersion.png",
#     width=6,height=6,units = 'in',
#     res=500)
plotDispEsts(dds_l23)
# dev.off()

sum(!mcols(dds_l23)$betaConv) #100 non-convergent genes

# png(filename = "final_project/SEA-AD/code/3_pseudoBulk_DE/weird_res_L2_3_IT/plot_pval_distribution.png",
#     width=6,height=6,units = 'in',
#     res=500)
hist(all_results[["L2_3_IT"]]$adjusted_model$pvalue, breaks = 50) #p-value distribution looks fine
# dev.off()

# Get non-convergence genes and compare with sig.genes
non_conv_genes <- rownames(dds_l23)[!mcols(dds_l23)$betaConv]
# Check overlap with DEGs
deg_genes <- all_results[["L2_3_IT"]]$adjusted_degs$gene
sum(deg_genes %in% non_conv_genes)
deg_genes[deg_genes %in% non_conv_genes]


#### Double check the counts between LATE #### 
# Function to plot normalized counts
gg_plotCounts = function(dds, gene, intgroup){
  library(ggplot2)
  d = plotCounts(dds = dds, gene = gene, intgroup = intgroup, 
                 normalized = TRUE, transform = TRUE, #default
                 returnData = TRUE)
  
  p = ggplot(d, aes(x=.data[[intgroup]], y=count, fill=.data[[intgroup]])) + 
    geom_violin(trim=FALSE, alpha=0.5, color=NA) +   # violin plot
    geom_boxplot(width=0.1, outlier.shape=NA) +  # boxplot inside violin
    geom_jitter(width=0.1, height=0, size=1.8, color="black", alpha = 0.5) +  # overlay points
    scale_y_log10() +
    scale_fill_manual(values = c("0" = "steelblue", "1" = "coral"),
                      labels = c("0" = "LATE-NC Negative", "1" = "LATE-NC Positive")) +
    theme_minimal(base_size = 13) +
    labs(y="Log normalized counts", title = gene) +
    theme(legend.position="none")
  
  return(p)
}
sig_gene = all_results[["L2_3_IT"]]$adjusted_model %>% 
  filter(padj < 0.05) %>% pull(gene)

for (gene in sig_gene){
  p = gg_plotCounts(dds_l23, gene=gene, intgroup="latenc_st2")
  print(p)
  ggsave(paste0("final_project/SEA-AD/code/3_pseudoBulk_DE/weird_res_L2_3_IT/",gene,".png"),
         plot = p, width = 5, height = 6, dpi = 500)
}


