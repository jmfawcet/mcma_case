# Bayesian meta-analysis of disorder prevalence with misclassification correction

This repository contains the simulation code, outputs, case-study analysis and figures associated with:

**Fawcett, Jonathan M., Whitridge, Jedidiah, Bartoš, František, & Fawcett, Emily J.**  
_Bayesian meta-analysis of disorder prevalence with misclassification correction_  
(Manuscript submitted for publication; a preprint link will be added when available)

The methods are implemented in the companion R package [**mcma**](https://github.com/jmfawcet/mcma).

---

## Study Overview

Meta-analyses of disorder prevalence routinely pool studies that used diagnostic interviews with 
studies that used screening questionnaires, yet rarely account for the imperfect sensitivity (*Se*) 
and specificity (*Sp*) of the screening instruments. Because false positives dominate when a disorder 
is rare, this practice inflates pooled prevalence. The paper develops a Bayesian hierarchical framework 
that embeds the Rogan–Gladen misclassification correction within a multilevel prevalence model, using 
informative priors on *Se* and *Sp* (or, when gold-standard studies are available, estimating them from 
the data) to propagate diagnostic uncertainty into the pooled estimate. The framework is evaluated 
in eight simulation studies and applied to a case study of depression prevalence among public safety 
personnel.

---

## Key Findings

- Naive meta-analysis produced severe upward bias in prevalence, with near-zero coverage, 
across the screening-only simulations
- The recommended corrected model yielded near-unbiased prevalence with appropriate coverage, even 
when a single field-level *Se*/*Sp* prior was used for all screening instruments
- Systematic overestimation of *Sp* was the primary vulnerability; including gold-standard 
studies partially mitigated it
- When gold-standard studies were available, screening *Sp* was recoverable from the data 
under a transportability assumption (i.e., studies using gold standard measures and screening 
measures share a common aggregate prevalence within a typical study)
- In a moderator analysis, the corrected model achieved greater statistical power than either a 
naive model or a model restricted to the available gold standard studies
- Re-analysing depression prevalence among public safety personnel reduced the naive estimate from 
12.2% to 2.7%, converging with the gold-standard-only estimate

---

## Repository Structure

```
mcma_case/
├── scripts/
│   ├── 0. Setup/               Model definitions, prior builders, data-generating
│   │                           processes, runners and diagnostics (0.0–0.5)
│   ├── 1. Simulations/         1.0 run the simulations; 1.1 summarise saved
│   │                           outputs; 1.2 figures and tables
│   ├── 2. Case Study/          2.0 case study fitted with the mcma package;
│   │                           2.1 concealment surface (run after 2.0)
│   └── 3. Motivating Heatmap/  Figure 1 (bias under misclassification)
├── stan/                       Stan programs for the supplementary M7 model
├── data/                       Downloaded case-study data (see Data below)
├── simulation_data/            ZIP archives of simulation outputs and generated
│                               datasets, plus summary spreadsheets
├── figures/                    Rendered figures for the paper and supplement
└── models/                     Not included: fitted brms objects are cached
                                here when the case-study script is run
```

---

## Reproducing the Analyses

Open `mcma_case.Rproj`. Extract the simulation archives as shown below before running 
the scripts in numbered order.

**Requirements.** R (the paper used R 4.5.0) with `brms`, `cmdstanr` and a working CmdStan 
installation, `posterior`, `rstan`, `psych`, the `tidyverse` packages, `patchwork`, `gghalves`, 
`scales`, `MASS`, `openxlsx`, and the companion package `mcma`:

```r
# install.packages("remotes")
remotes::install_github("jmfawcet/mcma")
```

**Simulations.** Saved per-replicate model summaries and generated datasets are
provided as ZIP archives under `simulation_data/`. Extract them once from the
repository root:

```r
archives <- list.files("simulation_data", pattern = "\\.zip$", full.names = TRUE)

for (archive in archives) unzip(archive, exdir = "simulation_data")
```

The tables and figures can then be regenerated without refitting:
`1.1 Process Simulation Results.R` rebuilds the summary spreadsheet, and
`1.2 Simulation Prevalence and Sp Figures.R` rebuilds the simulation figures,
including the supplementary tau figures, and prints the table summaries.
With convergence exclusions enabled, a failed model excludes its replicate
from all model summaries within that simulation.

**Submission figures.** Each figure script also writes a journal format copy of
every figure to `figures/submission/`

Refitting from scratch is done through `1.0 Run Simulations.R`; the generated
datasets are stored alongside the outputs (`sim*_dat.rds`). A full refit takes
days to weeks.

The code base was written by hand in stages and later refactored, with help from
Claude and ChatGPT, to simplify it. This was refactored again to create the mcma
package. The refactored code produces output identical to the original, except that 
the kernel-weighted grid averaging (KWGA) estimator was corrected afterwards (it now 
scores the unclamped Rogan-Gladen correction and resamples draws jointly with the 
*Se*/*Sp* grid); the KWGA-based Simulation 8 results and the case-study KWGA estimate 
were recomputed with the corrected estimator. All code was also reviewed and
re-tested by hand.

**Case study.** 

`2.0 Fit Case Study.R` fits the naive, gold-only, corrected and joint 
models and the *Se*/*Sp* sensitivity grid using `mcma`, caching each fitted object under 
`models/` (not included because of size). The fits take several hours on a desktop machine.

`2.1 Concealment Surface.R` then refits M8 over a grid of assumed concealment values (*c*, *γ*)
and draws the corrected prevalence surface (Figure S7).

---

## Data

The case study re-analyses the study-level depression prevalence data compiled by White et al. 
(2025) and shared on the Open Science Framework (https://osf.io/9y2xj/). The case-study script
downloads the data as `data/dep_psp.csv` when absent; internet access is required for the first download.

> White, N., Wagner, S. L., Matthews, L. R., Randall, C., Regehr, C., White, M., Alden, L. E., Buys, N., 
Carey, M. G., Corneil, W., Fyfe, T., Fraess-Phillips, A., & Krutop, E. (2025). Methodological correlates 
of variability in depressive and anxiety symptoms in public safety personnel: A systematic review and 
metaregression. *Traumatology, 31*(2), 346–361. https://doi.org/10.1037/trm0000538

---

## Citation

Fawcett, J. M., Whitridge, J., Bartoš, F., & Fawcett, E. J. 
*Bayesian meta-analysis of disorder prevalence with misclassification correction.* 
Manuscript submitted for publication.

---

## License

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

The code in this repository is released under the [MIT License](LICENSE). The data downloaded
from OSF remain subject to the terms of their original authors.
