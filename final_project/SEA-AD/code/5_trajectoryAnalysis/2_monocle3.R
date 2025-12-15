rm(list=ls())
library(monocle3)
library(Seurat)
library(harmony)
library(clustree)
library(tidyverse)
library(patchwork)

vlmc_seurat = readRDS("/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/data/MTG/vlmc_seurat.rds")
output_dir <- "/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/code/5_trajectoryAnalysis/out_vlmc_trajectory"

#### Step 1: Prepare metadata ####
colnames(vlmc_seurat@meta.data)
unique(vlmc_seurat@meta.data$`Donor ID`)

vlmc_subset = subset(
  vlmc_seurat, subset = 
         `Race (choice=White)` == "Checked" & 
         LATE %in% c("Not Identified", "LATE Stage 2", "LATE Stage 1", "LATE Stage 3")
)

vlmc_subset@meta.data = vlmc_subset@meta.data %>% mutate(
  latenc_st2 = if_else(LATE %in% c("Not Identified", "LATE Stage 1"), true = 0, false = 1)
)

colnames(vlmc_subset@meta.data)
unique(vlmc_subset@meta.data$`Donor ID`)

#### Step 2: Standard Seurat processing ####
vlmc_subset <- NormalizeData(vlmc_subset, 
                            normalization.method = "LogNormalize",
                            scale.factor = 10000,) #UP10K
vlmc_subset <- FindVariableFeatures(vlmc_subset, 
                                    selection.method = "vst",
                                    nfeatures = 3000) 
vlmc_subset <- ScaleData(vlmc_subset, features = VariableFeatures(vlmc_subset)) #default
vlmc_subset <- RunPCA(vlmc_subset, features = VariableFeatures(vlmc_subset), 
                       npcs = 50,
                       ndims.print = 1:5, nfeatures.print = 5)
# ElbowPlot(vlmc_subset, ndims = 50)
vlmc_subset <- RunHarmony(vlmc_subset, group.by.vars = c("Donor ID"))
vlmc_subset <- FindNeighbors(vlmc_subset, reduction = "harmony", dims = 1:30,
                             graph.name = c("harmony_nn", "harmony_snn")) #no more reduction in variance; ensure name change
vlmc_subset <- RunUMAP(vlmc_subset, reduction = "harmony", dims = 1:30,
                       reduction.name = "harmony_umap")
# Clustering
for (res in seq(0.1, 2, 0.1)){
  vlmc_subset <- FindClusters(vlmc_subset,
                              graph.name = "harmony_snn",
                              resolution = res, 
                              algorithm = 4, random.seed = 42) #leiden
}
colnames(vlmc_subset@meta.data)


#### Step 3: Visualization ####
#### clustree ####
# png(filename = file.path(output_dir,"clustree.png"),
#     width = 8, height = 12, units = "in", res = 500)
clustree(vlmc_subset, prefix = "harmony_snn_res.", 
         show_axis = TRUE, #show resolution
         # node_label = ,
         # use_core_edges = FALSE
         # node_colour = "sc3_stability", # stability of clusters across resolutions
)
# dev.off()

#### sil_scores ####
sil_scores <- sapply(seq(0.1, 2, 0.1), function(res) {
  res_col <- paste0("harmony_snn_res.",res)
  clusters <- vlmc_subset@meta.data[[res_col]]
  
  emb <- Embeddings(vlmc_subset[["harmony"]])[, 1:30] #clustering input
  # # exact
  # dist.matrix <- dist(x = emb)
  # sil <- cluster::silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
  
  bluster::approxSilhouette(emb, clusters) %>% 
    as.data.frame() %>% pull(width) %>% mean()
})

sil_scores

# png(filename = file.path(output_dir,"Silhouette.png"),
#     width = 6, height = 6, units = "in", res = 500)
plot(seq(0.1, 2, 0.1), sil_scores, type = "b",
     xlab = "Resolution", ylab = "Avg Silhouette")
# dev.off()

#### all resolution clustering patterns ####
for (res in seq(0.1, 2, 0.1)){
  res_col = paste0("harmony_snn_res.",res)
  message(res_col)
  # check number
  print( levels(vlmc_subset[[res_col]][,1]) %>% length() )
  # visualize
  p_res = DimPlot(vlmc_subset, reduction = "harmony_umap", label = TRUE,
                  group.by = res_col) + labs(title = res_col)
  print(p_res)
}


#### "harmony_snn_res.0.1" ####
Idents(vlmc_subset) <- "harmony_snn_res.0.1"
p1 <- DimPlot(vlmc_subset, reduction = "harmony_umap", 
              group.by = "harmony_snn_res.0.1",
              label = TRUE) +
  ggtitle("VLMC Clusters") +
  theme_classic()
p1

vlmc_subset@meta.data$latenc_st2 = factor(vlmc_subset@meta.data$latenc_st2)
p2 <- DimPlot(vlmc_subset, reduction = "harmony_umap", group.by = "latenc_st2") +
    scale_color_manual(values = c("0" = "blue", "1" = "red")) +
    ggtitle("VLMC by LATE-NC +/-") +
    theme_classic()
# DimPlot(vlmc_subset, reduction = "harmony_umap", split.by = "latenc_st2",
#         group.by = "latenc_st2")

p3 <- DimPlot(vlmc_subset, reduction = "harmony_umap", group.by = "Cognitive Status") +
    ggtitle("VLMC by Cognitive Status") +
    scale_color_manual(values = c("No dementia" = "skyblue", "Dementia" = "gold")) + 
    theme_classic()

p4 <- FeaturePlot(vlmc_subset, reduction = "harmony_umap", features = "VEGFA") +
  ggtitle("VEGFA Expression") +
  theme_classic()

combined_seurat <- (p1 | p4) / (p2 | p3)
# ggsave(file.path(output_dir, "umaps_metaInfo.pdf"),
#        combined_seurat, width = 14, height = 10)


# Create 4 groups to understand the pattern
analysis_df <- vlmc_subset@meta.data %>%
  mutate(
    # Group by LATE-NC and Dementia status
    disease_group = case_when(
      latenc_st2 == 0 & `Cognitive Status` == "No dementia" ~ "1. LATE- No Dementia",
      latenc_st2 == 1 & `Cognitive Status` == "No dementia" ~ "2. LATE+ No Dementia",
      latenc_st2 == 0 & `Cognitive Status` != "No dementia" ~ "3. LATE- Dementia",
      latenc_st2 == 1 & `Cognitive Status` != "No dementia" ~ "4. LATE+ Dementia",
      TRUE ~ "Other"
    )
  )
# Add VEGFA expression (logNorm)
vegfa_expr <- GetAssayData(vlmc_subset, slot = "data")["VEGFA", ]
analysis_df$VEGFA <- vegfa_expr[rownames(analysis_df)]
# Check distribution
print(table(analysis_df$disease_group))

# Plot VEGFA by the 4 groups
p_vegfa_groups <- ggplot(analysis_df %>% filter(disease_group != "Other"),
                         aes(x = disease_group, y = VEGFA, fill = disease_group)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggpubr::stat_compare_means(comparisons = list(
    c("2. LATE+ No Dementia", "3. LATE- Dementia")  # KEY comparison
  ), method = "wilcox.test") +
  scale_fill_manual(values = c(
    "1. LATE- No Dementia" = "#2E86AB",    # Blue - healthy
    "2. LATE+ No Dementia" = "#A23B72",    # Purple - LATE but resilient
    "3. LATE- Dementia" = "#F18F01",       # Orange - YOUR OBSERVATION
    "4. LATE+ Dementia" = "#C73E1D"        # Red - both pathologies
  )) +
  labs(
    title = "VEGFA Expression by LATE-NC and Dementia Status",
    subtitle = "Note: 'LATE- Dementia' (orange) has HIGH VEGFA",
    x = "Group",
    y = "VEGFA Expression"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
print(p_vegfa_groups)
# ggsave(file.path(output_dir,"VEGFA_four_groups_explanation.pdf"), 
#        p_vegfa_groups, width = 10, height = 7)

# VEGFA Mean by Group
vegfa_summary <- analysis_df %>%
  filter(disease_group != "Other") %>%
  group_by(disease_group) %>%
  summarise(
    n_cells = n(),
    mean_VEGFA = mean(VEGFA, na.rm = TRUE),
    sd_VEGFA = sd(VEGFA, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_VEGFA))

print(vegfa_summary)



#### Step 4: Convert to Monocle3 ###
# Extract data
expression_matrix <- GetAssayData(vlmc_subset, slot = "data")
cell_metadata <- vlmc_subset@meta.data %>% select('Donor ID', 'harmony_snn_res.0.1', 'Cognitive Status', 'latenc_st2',
                                                  'Age at Death','Sex','Years of education', #demographics
                                                  'APOE Genotype', 'Overall AD neuropathological Change', 'Overall CAA Score', 'Highest Lewy Body Disease', 'Total Microinfarcts (not observed grossly)', 'Atherosclerosis', 'Arteriolosclerosis', 'LATE')  
gene_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix), #required by monocle3
  row.names = rownames(expression_matrix)
)

# Create CDS
cds <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

#### Monocle3 preprocessing ####
# Transfer corrected embeddings + UMAP
reducedDims(cds)$PCA <- Embeddings(vlmc_subset, "harmony")
reducedDims(cds)$UMAP <- Embeddings(vlmc_subset, "harmony_umap")

#### cluster cells ####
cds <- cluster_cells(cds, reduction_method = "UMAP")
cat(sprintf("Found %d clusters\n", 
            length(unique(clusters(cds)))))
# Visualize
p_cluster <- plot_cells(cds, show_trajectory_graph = FALSE,
                        color_cells_by = "cluster",
                        label_cell_groups = FALSE,
                        cell_size = 0.5) +
  ggtitle("Monocle3 Clusters")

plot_cells(cds, show_trajectory_graph = FALSE,
           color_cells_by = "latenc_st2",
           label_cell_groups = FALSE,
           cell_size = 0.5)

#### Learn trajectory ####
cds <- learn_graph(cds)

#   The numbered points are the branching points and endpoints of this graph
#   (not meaningful yet)

p_traj <- plot_cells(cds, show_trajectory_graph = TRUE,
                     color_cells_by = "cluster",
                     label_cell_groups = FALSE,
                     label_leaves = FALSE,
                     label_branch_points = FALSE,
                     # graph_label_size = 4,
                     cell_size = 0.5) +
  ggtitle("Monocle3 Clusters with Trajectory Graph")
p_traj

# ggsave(file.path(output_dir,"monocle3_clusters_trajectoryGraph.pdf"),
#        p_traj, width = 8, height = 7)


#### order cells ####
# Choose root cells (from healthy/LATENC- donors)
root_cells <- colnames(cds)[colData(cds)$latenc_st2 == "0"]
cat(sprintf("Using %d LATENC- cells as root\n", length(root_cells)))

cds <- order_cells(cds, root_cells = root_cells)

# Plot pseudotime
p_pseudo1 <- plot_cells(cds, 
                        show_trajectory_graph = TRUE,
                        trajectory_graph_color = "red",
                        color_cells_by = "pseudotime",
                        label_cell_groups = FALSE,
                        label_leaves = FALSE,
                        label_branch_points = FALSE,
                        label_roots = FALSE,
                        cell_size = 0.5) +
  # scale_color_viridis_c() +
  ggtitle("Pseudotime")


p_pseudo2 <- plot_cells(cds,
                        show_trajectory_graph = TRUE,
                        trajectory_graph_color = "red",
                        color_cells_by = "latenc_st2",
                        label_cell_groups = FALSE,
                        label_leaves = FALSE,
                        label_branch_points = FALSE,
                        label_roots = FALSE,
                        cell_size = 0.5) +
  ggtitle("Latent Status")

p_pseudo_combined <- p_pseudo1 | p_pseudo2
# ggsave(file.path(output_dir,"pseudotime.pdf"),
#        p_pseudo_combined, width = 14, height = 6)


## Find trajectory genes (Finding genes changing along trajectory)
pr_graph_test_res <- graph_test(cds,
                                neighbor_graph = "principal_graph",
                                cores = 12)

pr_deg_ids <- subset(pr_graph_test_res, q_value < 0.05/nrow(pr_graph_test_res) & 
                                        status =='OK' &
                                        morans_I > 0.05) %>%
  arrange(desc(morans_I))

"VEGFA" %in% rownames(pr_deg_ids)

# # Save results
# # write_csv(pr_graph_test_res, 
# #           file.path(output_dir, "trajectory_genes_all.csv"))
# 
# write_csv(pr_deg_ids,
#           file.path(output_dir, "trajectory_genes_significant.csv"))

#### Examine VEGFA specifically ####
vegfa_result <- pr_graph_test_res %>%
  filter(gene_short_name == "VEGFA")


# Plot VEGFA
p_vegfa <- plot_genes_in_pseudotime(
  cds["VEGFA",],
  color_cells_by = "latenc_st2",
  min_expr = 0
) +
  ggtitle("VEGFA Expression Along Pseudotime") +
  theme_classic()

# ggsave(file.path(output_dir, "06_VEGFA_pseudotime.pdf"),
#          p_vegfa, width = 4, height = 4)
  
# Calculate correlation
vegfa_expr <- exprs(cds)["VEGFA",]
pseudo_vals <- pseudotime(cds)

complete_cases <- complete.cases(vegfa_expr, pseudo_vals)
cor_test <- cor.test(vegfa_expr[complete_cases], 
                     pseudo_vals[complete_cases],
                     method = "spearman")

cat(sprintf("  VEGFA-pseudotime correlation: ρ = %.4f, p = %.2e\n",
            cor_test$estimate, cor_test$p.value))


#### Other top Traj genes
top_genes <- pr_deg_ids %>%
  head(12) %>%
  pull(gene_short_name)

# Make sure VEGFA is included
if (!"VEGFA" %in% top_genes && "VEGFA" %in% rownames(cds)) {
  top_genes <- c("VEGFA", top_genes[1:11])
}

p_top_genes <- plot_genes_in_pseudotime(
  cds[top_genes,],
  color_cells_by = "latenc_st2",
  # min_expr = 0.5,
  ncol = 4
) + ggtitle("Other Top Trajectory Genes Expression Along Pseudotime") + 
  theme_classic()

# ggsave(file.path(output_dir, "07_top_trajectory_genes.pdf"),
#        p_top_genes, width = 12, height = 10)


# #### Export pseudotime values ####
# pseudotime_df <- tibble(
#   cell_id = colnames(cds),
#   pseudotime = pseudotime(cds),
#   cluster = clusters(cds)
# )
# # Add metadata
# pseudotime_df$latent_status <- colData(cds)$latenc_st2
# pseudotime_df$`Donor ID` <- colData(cds)$`Donor ID`
# 
# write_csv(pseudotime_df, file.path(output_dir, "cell_pseudotime_values.csv"))
# 
# # Save CDS object
# saveRDS(cds, file.path(output_dir, "monocle3_cds.rds"))
# saveRDS(vlmc_subset, file.path(output_dir, "vlmc_seurat.rds"))
