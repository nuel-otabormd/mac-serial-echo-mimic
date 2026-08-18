# Onset, progression and consequences of mitral annular calcification on serial clinical echocardiography

Analysis code for a retrospective cohort study of mitral annular calcification (MAC) in MIMIC-IV, from the BigQuery
extraction to every table and figure in the manuscript and its supplement. The pipeline is deterministic: two runs on the
same data produce identical results, and the final step checks the run against the values reported in the paper.

## What the study does

Adults with two or more transthoracic echocardiography episodes carrying the annular calcification field are followed from
their first (index) episode. Three questions are asked in order:

1. Which characteristics are associated with reaching moderate or severe MAC, from no reported calcification (onset) and
   from mild calcification (progression), under two outcome definitions (first-observed and confirmed moderate or severe
   read), with death as a competing event?
2. How reliable is the routinely reported grade? Transitions between consecutive reads, the fate of every first moderate or
   greater read, and a hidden Markov model that treats the underlying grade as unobserved.
3. What follows moderate or severe MAC: calcific mitral stenosis, moderate or greater mitral dysfunction, intervention, death.

The pre-specified analysis plan is `protocol/MAC_Study_Protocol_v2.4.md`.

## Data

The study uses only resources hosted on PhysioNet and mirrored on Google BigQuery in the `physionet-data` project:
MIMIC-IV v3.1 (hospital and ICU modules), MIMIC-IV-ECHO v1.0 (structured measurements), MIMIC-IV-ED v2.2 and MIMIC-IV-ECG
v1.0. All require credentialed access under the PhysioNet data use agreement. No patient-level data, extracted or derived,
are included in this repository, and `data/` and `outputs/` are ignored by git; the pipeline recreates them.

## Requirements

* Google Cloud SDK (`bq`; tested with BigQuery CLI 2.1.27), authenticated to an account with credentialed PhysioNet
  BigQuery access, and a Google Cloud project to bill the queries (about 8.5 GB scanned).
* R 4.5 (tested with 4.5.2) with `survival` (3.8-3), `mice` (3.18.0), `cmprsk` (2.2-12), `msm` (1.8-2), `sandwich` (3.1-1), `ggplot2` (4.0.0),
  `gridExtra` (2.3), `ragg` (1.5.0) and `scales` (1.4.0). Figures use the graphics device's default sans-serif family
  (Helvetica on macOS, Arial or DejaVu Sans elsewhere), so no additional fonts are needed; the R scripts that write text switch to a UTF-8 locale
  themselves, so they run correctly even from a shell that uses the C locale.
* `renv.lock` records the exact versions of these packages and their dependencies (82 packages, R 4.5.2) as used for the
  reported results. To reproduce that library, run `renv::restore(lockfile = "renv.lock")` in R before `run_all.sh`;
  a current CRAN library with the versions listed above also reproduces the results.

## How to run

```bash
export BQ_BILLING_PROJECT=<your-gcp-project-id>
bash run_all.sh
```

`run_all.sh` runs the ten steps below in order and stops at the first error. Total wall time on a laptop is about
two and a half hours, most of it in steps 4 (about 55 minutes) and 5 (about 70 minutes on five cores).

| Step | Script | Produces |
| --- | --- | --- |
| 1 | `sql/extract.sh` (runs `sql/01_analysis_frame.sql`, `sql/02_episode_panel.sql`) | `data/frame.csv` (one row per patient), `data/panel.csv` (one row per episode: defining-study date, first and last study, number of studies, inpatient flag, grade) |
| 2 | `R/01_prepare_frame.R` | `data/frame.rds` (derived variables: body mass index, CKD-EPI 2021 eGFR, categories, standardised markers, age at the qualifying echocardiogram; source-population and index-episode diagnostics printed) |
| 3 | `R/02_impute.R` | `data/mice.rds` (20 imputations of index covariates; seed 20260816) |
| 4 | `R/03_models_main.R` | `outputs/results/results_part1.rds` (co-primary piecewise-exponential models on person-intervals split at the time bands with patient-clustered variances, grouped-time complementary log-log (unsplit intervals) and Cox companions, confirmed-event timing and diagnosis-timing sensitivities, imputation bridge, onset-versus-progression interaction, stricter onset cohort, landmark at the second episode, complete-case models, rheumatic-report sensitivity, visit-process and censoring weights with diagnostics, total versus direct associations, Fine and Gray pooled across imputations, laboratory-definition and care-setting sensitivities) |
| 5 | `R/04_hmm.R` | `outputs/results/hmm.RData` (hidden Markov model from five dispersed starting sets, run in parallel with `parallel::mclapply`; set `MAC_CORES` to limit the cores, and note that on Windows `mclapply` runs the starts sequentially; about 70 minutes on five cores, longer sequentially) |
| 6 | `R/05_reliability_intervention.R` | `outputs/results/results_part2.rds` (baseline table, transitions, fate of first reads, model-derived quantities, stenosis incidence with the death-window sensitivity and cumulative incidence, landmarked model-independent check, dysfunction and intervention cascade, cumulative incidence of intervention, mortality at fixed horizons) |
| 7 | `R/06_secondary_mi.R` | `outputs/results/secondary_mi.rds` (stenosis time-varying-exposure model with the death-window sensitivity, cause-specific intervention model and survival model, all with multiply imputed covariates) |
| 8 | `R/07_figures.R`, `R/08_tables.R` | `outputs/figures/` (Figures 1–4 as PNG, TIFF and PDF at 600 dpi), `outputs/tables/` (Tables 1–4 and Supplementary Tables S1–S16 as CSV, with notes) |
| 9 | `R/09_verify.R` | regression check of this run against `expected/expected_results.csv`, a snapshot of the values reported in the manuscript taken from the archived run; exit status 1 if any check fails. It confirms that the code reproduces the reported numbers; it is not evidence that the numbers are right |
| 10 | `R/10_mi_stability.R` | `outputs/results/mi_stability.rds`, `outputs/tables/tableS17_imputation_stability.csv` (the imputation repeated with a second seed and the fully adjusted main models refitted, after checking that the archived imputations reproduce `results_part1.rds` exactly; about 15 minutes) |

`R/00_common.R` holds the shared conventions (death window, time bands, interval splitter) and is sourced by the scripts.

Steps 3–7 depend only on `data/`; step 5 is independent of steps 3, 4 and can run in parallel with them; step 10 needs steps 2, 3 and 4.

## Where each manuscript element comes from

| Manuscript | Source in `outputs/` |
| --- | --- |
| Figure 1 (cohort flow) | `figures/Figure1_flow.*`; counts from `results_part2.rds$flow` and `results_part1.rds$main` |
| Table 1 | `tables/table1_baseline.csv` |
| Table 2A, 2B | `tables/table2a_transitions.csv`, `tables/table2a_sequences.csv`, `tables/table2b_fate.csv` |
| Figure 2 | `figures/Figure2_reliability.*` (emission matrix and fate by age from `results_part2.rds`) |
| Table 3, Figure 3 | `tables/table3_main_models.csv`, `figures/Figure3_forest.*` |
| Figure 4, Table 4A | `figures/Figure4_calcific_stenosis_CIF.*`, `tables/table4a_stenosis_incidence.csv` and note |
| Table 4B, 4C | `tables/table4b_cascade.csv` and notes, `tables/table4c_intervention_or.csv` and notes |
| Supplementary Tables S1–S17 | `tables/tableS1_*` to `tables/tableS17*` (S7c weight diagnostics, S11 complementary log-log, S12 complete case, S13 landmark and rheumatic sensitivities, S14 death window, S15 hidden Markov starts, S16 fixed-horizon intervention and mortality, S17 imputation stability) |

## Reproducibility notes

* Every step is deterministic. One physical study defines each episode (highest grade, then earliest, then lowest
  identifier) and supplies its grade and, at index, the baseline covariates; time zero is the last study of the index episode;
  runs of studies whose calendar days overlap form one episode; the SQL resolves every other tie explicitly (most-recent and
  median laboratory values, first intervention, type of dysfunction); imputation and model fitting are seeded, and the hidden
  Markov model is a BFGS fit from five fixed starting sets (the fit with the lowest -2 log-likelihood among converged starts
  with a positive-definite Hessian is kept; the table of starts is Supplementary Table S15).
* Conventions: person-intervals are split at 0.5, 1, 2, 4 and 7 years; death within 30 days after the last echocardiogram is
  a competing event and later death is censoring at the last echocardiogram (`GRACE_DAYS` in `R/00_common.R`; varied in
  Supplementary Table S14); weights are truncated at 0.1 and 10.
* Patient identifiers in `data/` are hashes of the MIMIC-IV subject identifier and all times are day offsets from time
  zero (the last study of the index episode); these files still fall under the data use agreement and must not be shared.
* Random-number streams differ across R versions only if the default generator changes; the versions above reproduce
  the reported values exactly. Small last-digit differences in the hidden Markov model can arise from a different BLAS.

## Repository layout

```
sql/        BigQuery extraction (two queries and the wrapper script)
R/          shared conventions (00_common.R) and analysis scripts 01–10, run in order by run_all.sh
protocol/   the analysis plan as archived (v2.4) and DEVIATIONS.md, the dated note on how the final analysis differs from it
expected/   snapshot of the reported results used by the regression check in step 9
renv.lock   package versions used for the reported results
data/       created by step 1–3 (ignored by git)
outputs/    results, figures, tables and logs (ignored by git)
```

## Sharing the code

Share the repository only through git (`git push`) or a git export:

```bash
git archive --format=zip -o mac-serial-echo-mimic.zip HEAD
```

Both send exactly the tracked files. Do not zip the working folder itself: after a run it contains `data/` and
`outputs/`, which hold patient-level derived data covered by the PhysioNet data use agreement. If your working copy sits
in a cloud-synced folder, keep those two directories elsewhere and put symlinks named `data` and `outputs` in the
repository; the scripts only use the relative paths and `.gitignore` covers both forms.

## Citation and licence

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21983972.svg)](https://doi.org/10.5281/zenodo.21983972)

Code is released under the MIT licence (`LICENSE`). The repository lives at
https://github.com/nuel-otabormd/mac-serial-echo-mimic and every tagged release is archived on Zenodo. The concept DOI
10.5281/zenodo.21983972 resolves to the latest archived version; each release has its own version DOI, minted by Zenodo
when the release is published and shown on the Zenodo record and in the release notes on GitHub. The manuscript cites the
version DOI of the release it reports. Please cite the accompanying article and, for the code, the Zenodo record;
`CITATION.cff` carries the metadata.
