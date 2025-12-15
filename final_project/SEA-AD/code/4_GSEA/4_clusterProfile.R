#### GSEA on adjusted-model results (one run per cell type) ####
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Paths
base_dir <- "/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD"
res_dir  <- file.path(base_dir, "code/3_pseudoBulk_DE")
out_dir  <- file.path(base_dir, "code/4_GSEA/out_clusterProfile")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load combined DESeq2 results
all_results <- readRDS("/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/code/3_pseudoBulk_DE/de_results_lateNC_combined.rds")

#### Run GOBP enrichment for each cell type (ranked by log2FC) ####
run_gsea_by_celltype <- function(de_results, organism = "org.Hs.eg.db") {
  
  gsea_results <- map(names(de_results), function(ct) { #map to all celltypes
    
    # Prepare ranked (by log2FoldChange) gene list from adjusted model
    res <- de_results[[ct]]$adjusted_model %>%
      filter(!is.na(log2FoldChange)) %>%
      arrange(desc(log2FoldChange))
    
    gene_list <- res$log2FoldChange
    names(gene_list) <- res$gene
    
    # GSEA
    gsea_go <- gseGO(geneList = gene_list,
                     OrgDb = organism,
                     ont = "BP",
                     keyType = "SYMBOL",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH")
    
    return(gsea_go)
  })
  
  names(gsea_results) <- names(de_results)
  return(gsea_results)
}

gsea_results = run_gsea_by_celltype(all_results)
names(gsea_results)

#### Enrichment map showing pathway inter-relationships ####
plot_gsea_emaplot <- function(gsea_result,
                              cell_type = NULL,
                              top_n = 30,
                              min_edge = 0.2) {
  
  if (is.null(gsea_result) || nrow(gsea_result@result) == 0) {
    message("No significant pathways to plot")
    return(NULL)
  }
  
  # Pairwise term similarity
  gsea_result <- pairwise_termsim(gsea_result)
  
  p <- emapplot(gsea_result,
                showCategory = top_n,
                layout = "kk")
  
  if (!is.null(cell_type)) {
    p <- p + 
      ggtitle(paste("Pathway Enrichment Map -", cell_type)) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
  }
  
  print(p)
}

#### Cnetplot for a single cell type ####
plot_gsea_cnetplot <- function(gsea_result,
                               cell_type = NULL,
                               top_n = 5,
                               node_label_size = 3,
                               colorEdge = TRUE,
                               foldChange = NULL) {
  
  if (is.null(gsea_result) || nrow(gsea_result@result) == 0) {
    message("No significant pathways to plot")
    return(NULL)
  }
  
  # If foldChange vector is provided, use it for node coloring
  p <- cnetplot(gsea_result,
                showCategory = top_n,
                colorEdge = colorEdge,
                foldChange = foldChange)
  
  if (!is.null(cell_type)) {
    p <- p + 
      ggtitle(paste("Gene-Concept Network -", cell_type)) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
  }
  
  print(p)
}


ls_ct = c("Pvalb", "L4_IT", "L2_3_IT", "VLMC", "L6_CT", "Sst_Chodl", "Sst")

# pdf("final_project/SEA-AD/code/4_GSEA/out_clusterProfile/gsea_emaps_cnetplots_all.pdf", width = 12, height = 10)
for (ct in ls_ct){
  # print(ct)
  gsea_result = gsea_results[[ct]]
  plot_gsea_emaplot(gsea_result, ct)
  plot_gsea_cnetplot(gsea_result, ct, foldChange = gsea_result@geneList)
}
# dev.off()



# Compare enriched pathways across cell types
compare_pathways <- function(gsea_results, top_n = 20) {
  # Extract top pathways per cell type
  top_pathways <- map_dfr(names(gsea_results), function(ct) { #bind rows
    if (!is.null(gsea_results[[ct]]) && nrow(gsea_results[[ct]]@result) > 0) {
      gsea_results[[ct]]@result %>%
        as_tibble() %>%
        arrange(pvalue) %>%
        head(top_n) %>%
        mutate(cell_type = ct)
    }
  })
  
  # Heatmap of -log10(pvalue) for shared pathways
  pathway_matrix <- top_pathways %>%
    dplyr::select(cell_type, Description, pvalue) %>%
    mutate(neg_log_p = -log10(pvalue)) %>%
    pivot_wider(id_cols = Description, 
                names_from = cell_type, 
                values_from = neg_log_p,
                values_fill = 0)
  
  return(pathway_matrix)
}

compare_pathways(gsea_results)
