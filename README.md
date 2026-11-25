# SynechocystisSingleCell

Code and data for the single-cell metabolic flux analysis of *Synechocystis* sp. PCC 6803. Single-cell ¹³C and ¹⁵N incorporation (stable isotope labelling measured by NanoSIMS) is converted into per-cell uptake and growth rates. These rates then constrain genome-scale (iJN678) and spatial exchange metabolic models that are sampled to infer each cell's flux distribution.

The analysis is organised as numbered Jupyter notebooks (main text) and `S*` notebooks (supplementary information). They are meant to be run from the repository root, since the paths inside them are relative to it.

## Notebooks

| File | Description |
|---|---|
| `1_peak_movement_analysis_uptake_rates.ipynb` | Loads the NanoSIMS isotope data (blank, 4 h, 12 h, 24 h), splits the cells into low/high-labelling peaks (Gaussian mixture), and computes single-cell ¹³C and ¹⁵N incorporation rates by peak-relative projection. Writes the `peak_movement_rates_4to12h_*` files and the cluster ratio table. |
| `2_growth_rates_maxent.ipynb` | Uses the experimental growth data to estimate the average and maximum growth rate, corrects for the truncated distribution, applies maximum entropy to estimate single-cell growth rates, and turns replacement rates into C and N uptake rates (`results_std_corrected_12h_PEAKS.dat`, `max_growths_std_12h_PEAKS.dat`). |
| `3_metabolic_network_modeling.ipynb` | Builds the C/N uptake feasibility region of the iJN678 model (`feasibility_std*.dat`) and places the measured single cells inside it (2D/3D polytope plots). |
| `4_sampling_metabolism.ipynb` | Sets up the model (NH₄ diffusion, standard medium) and samples the genome-scale flux space of every cell with its measured uptake rates (1000 samples per cell, `std_1000samples_cell*.dat`). Also computes the biomass C:N ratio. |
| `5_genome_scale_results_no_exchange.ipynb` | Collects the per-cell samples into averages/modes/stds and analyses the genome-scale results without metabolite exchange (flux distributions, secretion patterns, C fluxes, UMAP, cluster statistics). |
| `6_acid_ex_sampling.ipynb` | Processes the output of the spatial exchange model (`exchange_model/sample_ox.dat`) into per-cell averages and errors (`results_std_corrected_12_toy_model_*`) and analyses exchange of organic acids/amino acids between cells. |
| `7_genome_scale_results_exchange.ipynb` | Genome-scale results when metabolite exchange is allowed: metabolic maps, organic carbon vs. growth, and polytopes with and without exchange. |
| `8_scatter_plots_fig4.ipynb` | Scatter plots and flux maps for Figure 4 (organic C exchange, fluxes with and without exchange, representative cells). |
| `9_cluster_comparison.ipynb` | Compares the two metabolic subpopulations of cells on their sampled flux profiles: per-reaction mean ± std by cluster, ranking of the reactions and subsystems that differ most between clusters. |
| `S1_peak_movement_analysis_uptake_rates_rep1_SI.ipynb` | Same analysis as notebook 1 applied to replicate 1 (Supplementary Information). |
| `S2_toy_model_correlations.ipynb` | Correlations between the toy-model fluxes (C, O₂, photons, acetate, growth, N), derived carbon and energy balance metrics (ATP supply, net fixed carbon, yields), compared between clusters. |
| `S3_interaction_network.ipynb` | Builds and plots the spatial interaction network between cells (who exchanges what with whom) from the cell positions and the exchange-model results. |

## Folders

### `exchange_model/`
Spatial metabolic exchange model in which cells can share metabolites with their neighbours.

| File | Description |
|---|---|
| `create_matrix.ipynb` | Builds the constraint matrix and flux bounds of the multi-cell exchange model from the cell positions and the single-cell uptake rates. |
| `autotroph_matrix.dat` | Constraint (stoichiometric) matrix of the exchange model, produced by `create_matrix.ipynb`. |
| `bounds.dat` | Flux bounds for the exchange model. |
| `cells_pos_3D.dat`, `cells_pos_3D_2.dat` | 3D positions of the cells (spatial configurations). |
| `sampleHRnew_spatial.cpp` | Hit-and-run sampler of the exchange model polytope. Adapted from the software of De Martino & Tkačik (*Statistical mechanics of the E. coli metabolic network in stationary growth*), distributed under GPL v3. |
| `rowechelon.cpp` | Row-echelon reduction routine used by the sampler. |
| `sample_ox.dat` | Output of the sampler. |
| `README.txt` | How to run the sampling. |

To run the sampling:

```bash
# 1) run create_matrix.ipynb
g++ -O9 sampleHRnew_spatial.cpp
./a.out > sample_ox.dat
```

### `experiments/`
Experimental data.

- `PCC6803_CN_Content_Lucia.xlsx`: C and N content measurements of *Synechocystis* PCC 6803.
- `new_data/`
  - `April2026.xlsx`, `April2026_H2.xlsx`, `April2026_ion_counts.xlsx`: raw NanoSIMS measurements (isotope ratios, H and ion counts).
  - `isotope_data.dat/.ods`, `isotope4h.dat`, `isotope12h.dat`: single-cell isotope labelling data, all timepoints and per timepoint.
  - `growth_data.dat/.ods`: OD730 and cell counts over time, used for the growth rate.
  - `OD_high_init_dens.dat`, `OD_low_init_dens.dat`: growth curves at high and low initial density.
  - `cluster_ratios_table.tex`: LaTeX table of the cluster ratios produced by notebook 1.
  - `replicates/`: per-replicate NanoSIMS data for each timepoint (`4h`, `12h`, `24h`, `rep1`/`rep2`, as `.ods` and filtered `.csv`), blanks (`blank.ods`, `blank_H.ods`) and H data (`12h_H.ods`).

### `models/`
Genome-scale metabolic models of *Synechocystis* PCC 6803.

- `iJN678.json`: original iJN678 model.
- `iJN678_std.json`: iJN678 adapted to the standard experimental medium.
- `iJN678_std_aminoacids.json`: standard model with amino acid exchange reactions.
- `iJN678_std_acid_exchange.json`: standard model with organic acid exchange reactions.
- `syn2.xml`: SBML model of *Synechocystis*.

### `results_new/`
Results produced by the notebooks.

- `peak_movement_rates_4to12h_{C,N}_{NEW,PERTURBED}.dat`: single-cell C and N incorporation rates (notebook 1).
- `clusters_4to12h_N_NEW.dat`: cluster (low/high labelling peak) assigned to each cell (notebook 1), used by most of the other notebooks to split the population into two subpopulations.
- `results_std_corrected_12h_PEAKS.dat`, `results_corrected_12h_PEAKS_PERTURBED.dat`: per-cell rates (`r_C`, `r_N`), uptake rates (`u_C`, `u_N`) and growth rate (notebook 2).
- `max_growths_std_*.dat`: maximum growth rates per cell.
- `feasibility_std.dat`, `feasibility_std_NEW.dat`: feasibility region in uptake space (notebook 3).
- `results_std_corrected_12_toy_model_{avg,std}_PEAKS[_ox|_4ox].dat`: per-cell average and std of the exchange-model fluxes (C, O₂, photons, acetate, growth...) for the different configurations (notebook 6).
- `genome_scale_samples_no_exchange_PEAKS/`: average, mode and std of the genome-scale flux samples of each cell without exchange (`std_sampling_exchange_{averages,modes,stds}.dat`).
- `genome_scale_samples_PEAKS_RuBisCO_ox/`: the same with exchange and RuBisCO oxygenation, plus the matching exchange-model results.
- `metabolic_map.csv`, `template_map.csv`, `3values.csv`: layout and values used to draw the metabolic maps.
- `map_7.png`, `metabolic_network.graphml_7.png`: metabolic map figures.

**Note:** the raw per-cell flux samples (`genome_scale_samples_no_exchange_PEAKS/std_1000samples_cell*.dat`, 1000 samples per cell, about 5 GB in total) are not included because of GitHub's size limits. They are produced by `4_sampling_metabolism.ipynb` and are available on request. Their summaries (averages, modes, stds) are included.

## Requirements

Python 3 with `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`, `scikit-learn`, `statsmodels`, `cobra`, `networkx`, `plotly`, `umap-learn`, `kneed` and `netgraph`, plus a C++ compiler (`g++`) for the exchange-model sampler.
