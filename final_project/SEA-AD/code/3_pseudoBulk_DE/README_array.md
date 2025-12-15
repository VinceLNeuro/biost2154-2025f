## After all jobs complete, combine results:
```bash
Rscript combine_array_results.R
```

This creates:
- `de_results_lateNC_combined.rds` - All results in one object
- `deg_summary.csv` - Summary table
- `deg_counts_barplot.png` - Visualization

## Output Files:

- `array_logs/` - Individual log files for each cell type
- `results_by_celltype/` - Individual RDS files per cell type
- `de_results_lateNC_combined.rds` - Final combined results

## Troubleshooting:

If a specific cell type fails:
```bash
# Rerun just that task (e.g., task 5)
sbatch --array=5 run_array.sh
```
