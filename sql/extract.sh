#!/usr/bin/env bash
# Extract the analysis frame and the episode panel from BigQuery into data/ (CSV).
# Requires: Google Cloud SDK (bq) authenticated to an account with credentialed PhysioNet access to
# MIMIC-IV v3.1, MIMIC-IV-ECHO v1.0, MIMIC-IV-ED and MIMIC-IV-ECG on the `physionet-data` project.
# Set the billing project in the environment: export BQ_BILLING_PROJECT=<your-gcp-project-id>
set -euo pipefail
: "${BQ_BILLING_PROJECT:?Set BQ_BILLING_PROJECT to the Google Cloud project that will be billed for the queries}"
here="$(cd "$(dirname "$0")" && pwd)"; root="$(dirname "$here")"; mkdir -p "$root/data"
run () {  # $1 = sql file, $2 = output csv
  echo "[extract] $(basename "$1") -> $(basename "$2")"
  bq --project_id="$BQ_BILLING_PROJECT" query --use_legacy_sql=false --format=csv --max_rows=100000000 --quiet < "$1" > "$2"
  echo "[extract]   rows: $(( $(wc -l < "$2") - 1 ))"
}
run "$here/01_analysis_frame.sql" "$root/data/frame.csv"
run "$here/02_episode_panel.sql"  "$root/data/panel.csv"
