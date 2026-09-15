# analise-recalque
Análise de deslocamentos e probabilidade de falha para recalques

## MATLAB subset simulation diagnostics

This repository now includes:

- `A9_subset_simulation_pf.m` — subset simulation driver that preserves the usual `Pf`, `beta` and `CoV` outputs and by default saves per-level diagnostic data to `C:\Kazuo-Script\out_incremental\A9_subset_simulation_levels.mat` (the base folder can be overridden with `opts.workDir` or `opts.outDir`, diagnostic persistence can be disabled with `opts.saveDiagnostics = false`, result persistence can be disabled with `opts.saveResult = false`, and the built-in conditional sampler assumes a standard-normal base space)
- `A9_plot_subset_simulation_diagnostics.m` — publication-style plotting utility for SS level clouds and the estimated failure boundary, with configurable background color, font sizes, marker sizes and line widths
