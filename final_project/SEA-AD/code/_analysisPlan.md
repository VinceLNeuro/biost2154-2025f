# Analysis Plan

## Background

The middle temporal gyrus (MTG) is part of the middle temporal cortex, a brain region affected at later stages of stage 3 LATE-NC (Limbic-predominant age-related TDP-43 encephalopathy – neuropathologic change). MTG tissue is therefore of interest for studying single-cell transcriptomic differences between LATE-NC positive (LATE-NC+) and negative (LATE-NC−) patients. Although LATE-NC is one of multiple pathological factors contributing to dementia, [??]_advanced LATE-NC stages correlate strongly with dementia risk, independent of amyloid or tau pathology_, highlighting its significance in neurodegeneration. However, the molecular and cellular transcriptional landscape of LATE-NC remains poorly characterized, particularly its independent effects. Identifying LATE-NC–associated transcriptional changes at the single-cell level, independent of other pathologies (and cognitive status), may reveal cell-type-specific risk genes and pathways contributing to disease progression, as well as potential targets (e.g. cell types or hub of molecular networks) for intervention, potentially slowing cognitive decline or protecting vulnerable neuronal populations.

## Key analysis

1. Cell-type-specific transcriptional changes
    - Are certain MTG cell types (e.g., excitatory neurons, inhibitory neurons, astrocytes, microglia, oligodendrocytes) differentially affected in __LATE-NC+ versus LATE-NC-__?
    - Approach: Control for cognitive status and co-pathologies, then perform differential expression (DE) analysis at the cell-type level.
    - Goal: Identify transcriptional programs specifically associated with LATE-NC, independent of other neuropathologies.

2. Glial activation trajectories (pseudotime / RNA velocity)
    - Question: Do glial cells in MTG show a _progression toward activated states_ in LATE-NC?
    - Approach: Subset microglia or astrocytes --> Infer pseudotime trajectories (e.g., diffusion pseudotime in Scanpy); Optionally overlay RNA velocity to infer directionality of transcriptional changes

3. Cell-cell communication changes
    - Question: Does LATE-NC alter interactions between neurons and glia in MTG?
    - Approach: Use tools such as _CellChat_, NicheNet, or Squidpy to compare ligand-receptor interactions between __LATE-NC+ and LATE-NC−__ cells.
    - Hypothesis: Increased glial signaling toward neurons (stress/inflammation)


4. miRNA-mediated regulation in aging populations with LATE-NC

   - Identified brain miRNAs DE in aging population with LATENC, independent of other neuropathologies, in previous studies. Now, we try to identify DE genes across cell-types. See how many significant _downregulated_ transcripts can be overlapped with found DE miRNAs targes within each cell-types (--> those transcripts are more likely the targets of miRNAs).

    - see their envolved pathways

    - Goal: validate prior results, and further provide mechanistic explanation of the post-transcriptional regulation that may contribute to LATE-NC.


5. Maybe can look into dementia + LATENC vs normal + LATENC (adjusting for other pathologies) --> which upregulated genes have protective effects on cognition. 

```python
# Subset cells from LATENC+ individuals:
adata_latenc = adata_mem[adata_mem.obs['LATE'] == 'positive', :].copy()

# Add covariates for modeling:
# - Dementia status (Cognitive Status) → outcome
# - Other pathologies → covariates

# Differential expression (DE) analysis:
# - Pseudo-bulk approach (recommended for multi-donor single-cell data):
# - Aggregate counts per donor per cell type → reduces donor-level confounding.
# - **Compare dementia vs normal within LATENC+.
# - Model: Expression ~ Cognitive_Status + Other_Pathologies + Age + Sex + PMI

# Do cell-type specific analysis
```

<br>

## Additional notes

dataset: `SEA-AD/data/MTG/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad`
