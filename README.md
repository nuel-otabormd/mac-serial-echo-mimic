# Onset, progression and consequences of mitral annular calcification on serial clinical echocardiography

Analysis code for a retrospective cohort study of mitral annular calcification (MAC) in MIMIC-IV, from the BigQuery
extraction to every table and figure in the manuscript and its supplement. The pipeline is deterministic: two runs on the
same data produce identical results, and the final step checks the run against the values reported in the paper.

## What the study does

Adults with two or more transthoracic echocardiography episodes carrying the annular calcification field are followed from
their first (index) episode. Three questions are asked in order:

1. Which characteristics are associated with reaching moderate or greater MAC, from no reported calcification (onset) and
   from mild calcification (progression), under two outcome definitions (first-observed and confirmed moderate or greater
   read), with death as a competing event?
2. How reliable is the routinely reported grade? Transitions between consecutive reads, the fate of every first moderate or
   greater read, and a hidden Markov model that treats the underlying grade as unobserved.
3. What follows moderate or greater MAC: calcific mitral stenosis, moderate or greater mitral dysfunction, intervention, death.

The pre-specified analysis plan is `protocol/MAC_Study_Protocol_v2.4.md`.

## Data

The study uses only resources hosted on PhysioNet and mirrored on Google BigQuery in the `physionet-data` project:
MIMIC-IV v3.1 (hospital and ICU modules), MIMIC-IV-ECHO v1.0 (structured measurements), MIMIC-IV-ED v2.2 and MIMIC-IV-ECG
v1.0. All require credentialed access under the PhysioNet data use agreement. No patient-level data, extracted or derived,
are included in this repository, and `data/` and `outputs/` are ignored by git; the pipeline recreates them.

## Requirements

* Google Cloud SDK (`bq`; tested with BigQuery CLI 2.1.27), authenticated to an account with credentialed PhysioNet
  BigQuery access, and a Google Cloud project to bill the queries (about 8.5 GB scanned).
* R 4.5 (tested with 4.5.2) with `survival` (3.8-3), `mice` (3.18.0), `cmprsk` (2.2-12), `msm` (1.8-2), `ggplot2` (4.0.0),
  `gridExtra` (2.3), `ragg` (1.5.0) and `scales` (1.4.0). Figures use the graphics device's default sans-serif family
  (Helvetica on macOS, Arial or DejaVu Sans elsewhere), so no additional fonts are needed; the R scripts that write text switch to a UTF-8 locale
  themselves, so they run correctly even from a shell that uses the C locale.

## How to run

```bash
export BQ_BILLING_PROJECT=<your-gcp-project-id>
bash run_all.sh
```

`run_all.sh` runs the nine steps below in order and stops at the first error. Total wall time on a laptop is about
one and a half hours, most of it in steps 4 and 5.

| Step | Script | Produces |
| --- | --- | --- |
| 1 | `sql/extract.sh` (runs `sql/01_analysis_frame.sql`, `sql/02_episode_panel.sql`) | `data/frame.csv` (one row per patient), `data/panel.csv` (one row per episode) |
| 2 | `R/01_prepare_frame.R` | `data/frame.rds` (derived variables: body mass index, CKD-EPI 2021 eGFR, categories, standardised markers) |
| 3 | `R/02_impute.R` | `data/mice.rds` (20 imputations of index covariates; seed 20260816) |
| 4 | `R/03_models_main.R` | `outputs/results/results_part1.rds` (co-primary piecewise-exponential models, imputation bridge, onset-versus-progression interaction, stricter onset cohort, informative-observation and censoring sensitivities, total versus direct associations, Fine and Gray companions, laboratory-definition and care-setting sensitivities) |
| 5 | `R/04_hmm.R` | `outputs/results/hmm.RData` (hidden Markov model; about 20 minutes) |
| 6 | `R/05_reliability_intervention.R` | `outputs/results/results_part2.rds` (baseline table, transitions, fate of first reads, model-derived quantities, stenosis incidence, model-independent check, intervention cascade) |
| 7 | `R/06_secondary_mi.R` | `outputs/results/secondary_mi.rds` (stenosis time-varying-exposure model, intervention logistic model and survival model, all with multiply imputed covariates) |
| 8 | `R/07_figures.R`, `R/08_tables.R` | `outputs/figures/` (Figures 1–4 as PNG, TIFF and PDF at 600 dpi), `outputs/tables/` (Tables 1–4 and Supplementary Tables S1–S10 as CSV, with notes) |
| 9 | `R/09_verify.R` | comparison of this run with `expected/expected_results.csv`, the values reported in the manuscript; exit status 1 if any check fails |

Steps 3–7 depend only on `data/`; step 5 is independent of steps 3, 4 and can run in parallel with them.

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
| Supplementary Tables S1–S10 | `tables/tableS1_*` to `tables/tableS10*` |

## Reproducibility notes

* Every step is deterministic. The SQL resolves ties explicitly (index study within an episode, ordering of same-day
  episodes, most-recent and median laboratory values, first intervention, type of dysfunction), imputation and model
  fitting are seeded, and the hidden Markov model is a BFGS fit from fixed starting values.
* Patient identifiers in `data/` are hashes of the MIMIC-IV subject identifier and all times are day offsets from the
  index episode; these files still fall under the data use agreement and must not be shared.
* Random-number streams differ across R versions only if the default generator changes; the versions above reproduce
  the reported values exactly. Small last-digit differences in the hidden Markov model can arise from a different BLAS.

## Repository layout

```
sql/        BigQuery extraction (two queries and the wrapper script)
R/          analysis scripts 01–09, run in order by run_all.sh
protocol/   the pre-specified analysis plan
expected/   manifest of the reported results used by the verification step
data/       created by step 1–3 (ignored by git)
outputs/    results, figures, tables and logs (ignored by git)
```

## Citation and licence

Code is released under the MIT licence (`LICENSE`). Please cite the accompanying article; a DOI for this repository is
provided in `CITATION.cff` once archived.
