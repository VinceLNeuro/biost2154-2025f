# Workflow - SEA-AD LATE-NC Analysis

snRNA dataset: `SEA-AD/data/MTG/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad`
- 20 cores for backend mode
- 23 cores for full load mode (DE) (11 mins to load)

<br>

__Table of Contents__
- [Workflow - SEA-AD LATE-NC Analysis](#workflow---sea-ad-late-nc-analysis)
    - [0.1. EDA and Data Conversion](#01-eda-and-data-conversion)
    - [0.2. More Exploration for the data and Clustering QC](#02-more-exploration-for-the-data-and-clustering-qc)
        - [Seurat](#seurat)
        - [Scanpy](#scanpy)
        - [ SEA-AD data introduction](#-sea-ad-data-introduction)
    - [1. Batch QC](#1-batch-qc)
        - [Why without integration, UMAP (without cluster annotation) will look problematic for different batches?](#why-without-integration-umap-without-cluster-annotation-will-look-problematic-for-different-batches)
    - [2. Clustering QC](#2-clustering-qc)
        - [Manual Annotation Approach (Nature Mathys 2024 \> Methods):](#manual-annotation-approach-nature-mathys-2024--methods)
- [Marker genes list database (cross-atlas consensus)](#marker-genes-list-database-cross-atlas-consensus)
    - [Reference-based mapping](#reference-based-mapping)
    - [Continue 2. Clustering QC](#continue-2-clustering-qc)
    - [3. Subsetting Data and Output for PseudoBulk](#3-subsetting-data-and-output-for-pseudobulk)
        - [3.1. Prefiltering \& Aggregation](#31-prefiltering--aggregation)
            - [Filter (*Jiebiao\_07* pseudobulk aggregation)](#filter-jiebiao_07-pseudobulk-aggregation)
            - [Aggregation uses `sc.get.aggregate`](#aggregation-uses-scgetaggregate)
    - [4. PseudoBulk in DESeq2](#4-pseudobulk-in-deseq2)
    - [5. Continue output nebula matrices in `/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/code/2_subset4pseudoBulk/2_subset4pseudoBulk.ipynb`](#5-continue-output-nebula-matrices-in-ihomehparktil177githubbiost2154-2025ffinal_projectsea-adcode2_subset4pseudobulk2_subset4pseudobulkipynb)

<br>

## 0.1. EDA and Data Conversion

- Tried `zellkonverter::readH5AD` in `code/0_eda/0_eda_scRNA.R`: failed; error within wrapper function

- Tried `schard` in `code/0_eda/0_eda_scRNA.R` for direct h5ad-data-read into seurat object. However, <mark> error message: </mark>
    - ```
        > data = schard::h5ad2seurat(h5ad_file)
        
        data = schard::h5ad2seurat(h5ad_file) Error in h5ad2Matrix(filename, expname, use_spam = use_spam) : The object you are trying to load is too large for Seurat and R in general: it has more than (2^31 -1) non-zero values in expression matrix. Consider setting use_spam=TRUE or use python. For more information please check: 1. https://github.com/cellgeni/schard/issues/1 2. https://github.com/chanzuckerberg/cellxgene-census/issues/1095 In addition: Warning message: In h5ad2data.frame(filename, "obs") : unexpected data.frame format, some columns can be missed for this kind of data
      ```
    - `sceasy` usage: https://github.com/cellgeni/sceasy?tab=readme-ov-file
    - `schard` usage: https://github.com/cellgeni/schard, https://github.com/cellgeni/schard/issues/1

<br>

__Files/Dirs in the `code/0_eda/`__

- `env_zellkonverter/`: `zellkonverter` conda environment specs
- `log/`              : All errors in the code ran above
- `0_eda_scRNA.R`     : _(archived)_ for data conversion --> will be refined as `0_h5ad2seuratRDS.R` in the future

- **`run.sh`          : As a template for resource allocation, and submit _Rscript_ in slurm job
- **`0_eda.ipynb`     : Python exploration for snRNAseq data <mark> before clustering QC </mark>

<br>

__Conclusion__

1. The large size of the sparse matrix limits the pre-DE analysis in `R` —> do these in `scanpy`

2. For statistical DE analysis, <mark> split h5ad by cell type (-> can convert to seuratObj -> `nebula` cell-level analysis) & by peudobulk (-> `nebula`, `DESeq2`, `monocle3`, etc) </mark>

<br>


## 0.2. More Exploration for the data and Clustering QC

__Follow [clustering QC from Harvard](https://github.com/hbctraining/Intro-to-scRNAseq/blob/master/lessons/08_SC_clustering_quality_control.md)__

- Go through the QC -> Normalization (regress out unwanted var) -> (integration) -> Dimension Reduction PCA + Clustring pipeline -> find marker (checking) 

- Compare statistical aspects with seurat pipeline

<br>

### Seurat

1. QC _[Jiebiao_05]_
    1. detected genes (high/low)
    2. library size
    3. mito fraction
2. Normalization _[Jiebiao_05]_
    1. __[linear_reg]__ LogNormalize (scale.factor = 10,000, margin = 1) -> FindVariableFeatures -> [(Cell-Cycle Scoring (diff?))](https://satijalab.org/seurat/articles/cell_cycle_vignette.html) -> [ScaleData](https://satijalab.org/seurat/reference/scaledata)
    2. __or [NB]__ [SCTransform](https://satijalab.org/seurat/reference/sctransform) --> pearson residual for PCA
3. PCA _[Jiebiao_05]_ + clustering _[Jiebiao_06]_
    1. RunPCA, (IntegrateLayers), FindNeighbors (GKNN), FindClusters (Leiden clustering), RunUMAP
        - PrepSCTFindMarkers [https://satijalab.org/seurat/articles/integration_introduction#perform-analysis-without-integration]

    <br>

    <mark>Notes for dimension reduction & clustering</mark>

    - Step 1: Build the <span style="color:red">kNN graph</span> in <span style="color:red">lower-dimensional __(PC1:n)__ space</span>, connecting each cell to its <span style="color:red">__k nearest neighbors__</span> (e.g., n_neighbors=15)

    - Step 2: <span style="color:red">Leiden</span> is where clustering really happens (operate on neighbor graph), optimizing (moving the nodes) the KNN graph based on <span style="color:red">__resolution__</span>, to find the best clustering/community result <span style="color:red">in PC space</span>. 

        - The resolution parameter controls the granularity:
            - Low resolution → fewer, larger clusters
            - High resolution → more, smaller clusters

    - Step 3: UMAP for visualization

        - We reduce dimension further from <span style="color:red">the PC space (need the "highD graph" based on n_neighbors)</span> --> <span style="color:red">embed into 2D space ("lowD graph")</span> by UMAP, retaining the similar strucutre using topological representation.
        
            - <mark>Thus, UMAP is only affected by nPC & n_neighbour & distance function</mark>
        
        - Clusters from Leiden can then be overlaid on the UMAP for visual inspection, but UMAP itself does not define clusters.

    - Step 4: Annotation

        - Clusters are then annotated using:
            - Marker genes
            - Module scores
            - Other metadata (cell type, sample, etc.)


### Scanpy

0. Read count matrix into `anndata.AnnData` object -> direct print for summary info
    - `AnnData object with n_obs × n_vars` (cell by gene matrix)

        Notes: 
        ```python
        #check number of donors = 89
        len(adata.obs['Donor ID'].unique())

        #similar to `merge` in R 
        ad.concat(adatas, label="sample")
        ```

    *** Benchmarking: __~8 mins__ to load SEAAD snRNAseq data (`adata = sc.read_h5ad("SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad", backed='r')`)


1. QC/filter (per cell) - SAME as R
    1. detected genes (high/low)
    2. library size (total counts per cell)
    3. mito fraction
    4. __*doublet detection__
    
2. Normalization and Transformation
    1. linear as default
    2. [Pearson residual](https://scanpy.readthedocs.io/en/stable/tutorials/experimental/pearson_residuals.html#)

3. PCA > KNN graph > calculate UMAP > plot UMAP > _(check if need to integration, if so switch to `scVI`; skipped norm+scale+pca (use raw counts) and redo sc.pp.neighbors and so on)_ > Leiden clustering (multiple resolution) > replot UMAP
    - [`scvi-tools`](https://scvi-tools.org) for batch integration.


<br>

### <mark> SEA-AD data introduction
- __1.3M nuclei, 36k genes__ (n_obs × n_vars = 1378211 × 36601) from __89 donors__

    - Raw UMI counts stored in `adata.layers['UMIs']`

    - Noncoding RNAs
        - hsa-mir-1253, hsa-mir-423
        - Number of lnc_MIRHG: 64 
        - Number of long intergenic non-coding: 2086
        - Also, contains AS lncRNA

    - QC metric (inferred from `obs`): 
        - keep `pct_counts_mt < 5`
        - min(n_gene) > 1000
        - scrublet doublet_score < 0.7 (Median = 0.0408)
        - Others not sure...

        |        | n_genes_by_counts | total_counts | pct_counts_mt 
        | :----  | :---------------- | :----------- | :------------
        | median | 5428.0            | 18642.0      |	0.138789
        | min    | 1001.0            | 1221.0       |	0.000000
        | max    | 15607.0           | 511016.0     |	4.999259


    - Integration: `scVI`
        - batch_key: `uns['batch_condition']` -> 'Specimen ID' (likely 'Donor ID')
        
        - confirmed by `adata.uns["neighbors"]["params"]`

<br>
<hr>

## 1. Batch QC

```
Main **technical variables** (batch effects) <- from colnames of `important_cols`

['Donor ID', 'method', 'Doublet score', 'Number of UMIs', 'Genes detected', 'Fraction mitochondrial UMIs']
```

### Why without integration, UMAP (without cluster annotation) will look problematic for different batches? 

If batch effect is not removed (without integration), PC space will capture the batch effects as top/main variance (confounded!) 

--> KNN graph built on the PCs will relfect batch similarity 

--> UMAP is built on gKNN, which will also be problematic!!! (not yet even clustering)

- Visualization Example: https://satijalab.org/seurat/articles/integration_introduction#perform-integration

<br>

Here, we will check if the batch has been well removed by `scVI`. <span style="color:lightgreen">Whether the same UMAP structure appears across donors/methods, or unique community pattern?</span>

```python
for cat in ['Donor ID','method']:
    sc.pl.umap(adata, color = cat) #altogether plot -> look for sharing pattern
    current_data = adata.obs[cat]
    for split in adata.obs[cat].unique().tolist():
        sc.pl.umap(adata[current_data == split], color = "Subclass", #use Subclass for easier visualization to ensure: every community are there in each sample
                   title=split)

for num in ['Doublet score','Number of UMIs','Genes detected','Fraction mitochondrial UMIs']:
    sc.pl.umap(adata, color = num, color_map = 'magma') #look for warmer color (higher value)
```

- No batch effects across donors
- No batch effects across (sequencing) methods
- No problematic QC metrics

<br>


## 2. Clustering QC

Starting with `code/1_clusteringQC`: Further clustering QC > Marker Identification & Manual Annotation


### Manual Annotation Approach (Nature Mathys 2024 > [Methods](https://www.nature.com/articles/s41586-024-07606-7#Sec10)):

>Third, we determined cell type marker genes based on data [published by the Allen Brain Institute](https://www.nature.com/articles/s41586-019-1506-7#Sec2) using the FindAllMarkers function from Seurat (Wilcoxon rank-sum test with Bonferroni correction for multiple testing; Padj < 0.05) and __computed module scores for each cell type marker gene set__ across all neuronal cells analysed in this study using the AddModuleScore function of Seurat.

- Check Mathys et al. (2024) > SupplFigs > __Fig.S1__ (known cell type marker genes in each major cell type) + __Tab.S2__ 

- Check [Allen Brain Institute paper (Nature 2019)](https://www.nature.com/articles/s41586-019-1506-7#Sec2) > Methods > "Combinatorial cell-type marker genes" > S2

<br>
<hr>

# Marker genes list database (cross-atlas consensus)

- [Azimuth Human - Motor Cortex](https://azimuth.hubmapconsortium.org/references/#Human%20-%20Motor%20Cortex) > Annotation Details > subclass

- [CellMarker](http://117.50.127.228/CellMarker/CellMarkerSearch.jsp?index_species=Human&index_tissue=Brain&index_cellname=Astrocyte&index_key=1#cloud_word) > Cell Tools > Cell annotation (-> Tissue = Cortex/Neocortex)

- Above, Allen brain Institute MTG combinatorial marker > S2
- Above, Mathys et al. (2024) > SupplFigs > __Fig.S1__ (known cell type marker genes in each major cell type) + __Tab.S2__***

## Reference-based mapping

- [BrainMap-MapMyCells](https://knowledge.brain-map.org/mapmycells/process?refTaxonomyId=10x_whole_human_brain)

- `Azimuth`

- `CellTypist` (Immune cells & Embryonic cells)
    ```python
    import scanpy as sc
    import celltypist
    model = celltypist.models.download_model('Immune_All_Humans')
    predictions = celltypist.annotate(adata, model=model)
    ```

<hr>

## Continue 2. Clustering QC

- Focused on `Subclass` cell types
- Applied: 
    1. `Mathys_mix_marker` (single gene)
        
        - Dot plot & UMAP well matched
           
            - __Well-clustered based on marker genes__: Ext and Inh (too detailed in SEA-AD), Ast, Oli, OPC, [Mic/Immune](https://www.nature.com/articles/s41467-021-25595-3/figures/1), End, VLMC (Vascular Leptomeningeal Cells) (=Per+SMC+Fib, makes sense)
            - T cell marker is distirbuted all over neuronal types
            - Csf-realted cell marker is NOT performing well. 
    
    2. [subset 1%] `module_Azimuth_motor`, `module_CellMarker_brain`
        
        - Rationale: Since a large focus on the downstream analysis will focus on __neuron and glial cells (also Endo__ given strong relevance and great annotation evidence above), we can use module score to further validate the clustering and annotation.
        
        - Conclusion: 
            - Azimuth MS for ----> Ext & Inh
            - Both MS good for --> Ast, Oli, OPC, Micro, Endo

<br>

## 3. Subsetting Data and Output for PseudoBulk

Starting with `code/2_subset4pseudoBulk`

### 3.1. Prefiltering & Aggregation

#### Filter (*Jiebiao_07* pseudobulk aggregation)
    
* __Raw UMIs__

* __[Filtered]__ ≥ 20 cells per sample/cell type
    * see `SEA-AD/data/MTG/pseudoBulk_20251109/valid_groups__per_DonorSubclass.csv` for later filtering at PB level
    
    * ```
      >pb_counts
      AnnData object with n_obs × n_vars = 2121 × 36601
        obs: 'Donor ID', 'Subclass'
        var: 'gene_ids'
        layers: 'sum'
      ```

    * n_obs = 2121 -> 1933 if use the filter (89*24 = 2136, but seems that some donors lack certain cell subtypes)

* __[All-Passed]__ ≥ 3 biological replicates per condition

<br>

#### Aggregation uses `sc.get.aggregate`

* Outputs:
    * `SEA-AD/data/MTG/pseudoBulk_20251109/valid_groups__per_DonorSubclass.csv`
    * `SEA-AD/data/MTG/pseudoBulk_20251109/pseudobulk_counts.csv`
    * `SEA-AD/data/MTG/pseudoBulk_20251109/pseudobulk_metadata.csv`

<br>
<hr>


## 4. PseudoBulk in DESeq2

Starting with `code/3_pseudoBulk_DE`

* *Jiebiao_07* practical considerations for filtering during pseudobulk DE
* 


## 5. Continue output nebula matrices in `/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/code/2_subset4pseudoBulk/2_subset4pseudoBulk.ipynb`