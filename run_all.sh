#!/usr/bin/env bash
# Reproduce the study end to end. Requirements: Google Cloud SDK (bq) with credentialed access to the PhysioNet
# BigQuery project, R >= 4.3 with the packages listed in README.md, and BQ_BILLING_PROJECT set in the environment.
# Steps 1-6 write to data/ and outputs/results/; 7-8 write figures and tables; 9 verifies the run against the manifest.
set -euo pipefail
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8   # R reads the scripts and writes text as UTF-8 regardless of the shell locale
cd "$(dirname "$0")"; mkdir -p data outputs/results outputs/figures outputs/tables outputs/logs
step () { echo "[$(date +%H:%M:%S)] $1"; }
step "1/9 extracting analysis frame and episode panel from BigQuery"; bash sql/extract.sh
step "2/9 preparing analysis frame";            Rscript R/01_prepare_frame.R          2>&1 | tee outputs/logs/01_prepare_frame.log
step "3/9 multiple imputation (m = 20)";        Rscript R/02_impute.R                 2>&1 | tee outputs/logs/02_impute.log
step "4/9 main models and sensitivity analyses"; Rscript R/03_models_main.R           2>&1 | tee outputs/logs/03_models_main.log
step "5/9 hidden Markov model of the reported grade"; Rscript R/04_hmm.R              2>&1 | tee outputs/logs/04_hmm.log
step "6/9 reliability, descriptives, intervention cascade"; Rscript R/05_reliability_intervention.R 2>&1 | tee outputs/logs/05_reliability_intervention.log
step "7/9 secondary models with multiple imputation"; Rscript R/06_secondary_mi.R     2>&1 | tee outputs/logs/06_secondary_mi.log
step "8/9 figures and tables";                  Rscript R/07_figures.R                2>&1 | tee outputs/logs/07_figures.log
                                                Rscript R/08_tables.R                 2>&1 | tee outputs/logs/08_tables.log
step "9/9 verification against the reported results"; Rscript R/09_verify.R           2>&1 | tee outputs/logs/09_verify.log
step "done"
