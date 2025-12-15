"""
Subset VLMC cells from SEA-AD MTG dataset and save as H5AD for downstream trajectory analysis in R.
"""

from __future__ import annotations #default now for name.error issue
import os
import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
SEED = 1234
random.seed(SEED)
np.random.seed(SEED)

import scanpy as sc
import anndata as ad
from collections import Counter #count table

sc.set_figure_params(dpi=300, color_map="viridis_r", facecolor="white", )
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()

# Path to original h5ad file
h5ad_path = "/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/data/MTG/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"

# Load the full dataset
print("Loading full MTG dataset...")
adata = sc.read_h5ad(h5ad_path, backed='r')

print(f"Total cells: {adata.shape[0]}")
print(f"Total genes: {adata.shape[1]}")
print(f"\nCell types available:")
print(adata.obs['Subclass'].value_counts())

#### Subset VLMC cells ####
vlmc_cells = adata[adata.obs['Subclass'] == 'VLMC'].to_memory()

print(f"\nVLMC cells found: {vlmc_cells.shape[0]}")
print(f"Genes: {vlmc_cells.shape[1]}")

# Check donor distribution
print("\nDonor distribution in VLMC subset:")
print(vlmc_cells.obs['Donor ID'].value_counts().head(10))

# Check cognitive status if available
if 'Cognitive Status' in vlmc_cells.obs.columns:
    print("\nCognitive Status distribution:")
    print(vlmc_cells.obs['Cognitive Status'].value_counts())


#### Save VLMC subset as H5AD ####
# Sanitize obs column names (HDF5 disallows forward slashes in keys)
vlmc_cells.obs.columns = [c.replace('/', '_') for c in vlmc_cells.obs.columns]

output_path = "/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/data/MTG/SEAAD_MTG_RNAseq_final-nuclei_VLMC_subset.h5ad"
vlmc_cells.write_h5ad(output_path)
print(f"\nVLMC subset saved to: {output_path}")

# Also save as CSV for quick inspection
metadata_output = "/ihome/hpark/til177/GitHub/biost2154-2025f/final_project/SEA-AD/data/MTG/VLMC_metadata.csv"
vlmc_cells.obs.to_csv(metadata_output)
print(f"VLMC metadata saved to: {metadata_output}")
