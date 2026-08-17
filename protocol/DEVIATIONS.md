# Deviations from, and clarifications to, the archived analysis plan

`MAC_Study_Protocol_v2.4.md` is reproduced as written on 16 August 2026 and has not been edited. The final analysis in
this repository differs from it in the respects below. None of them changes an outcome definition, an exposure, a
model, or a reporting rule.

1. **Data sources.** The plan lists MIMIC-IV-ECHO v0.1 and MIMIC-IV-Note. The structured echocardiography measurements
   used come from MIMIC-IV-ECHO version 1.0 (the version identity of the PhysioNet resource was verified after the plan
   was written; v0.1 is a smaller subset that was never used). MIMIC-IV-Note was not used: mitral interventions were
   identified from ICD-9-CM and ICD-10-PCS procedure codes only. MIMIC-IV-ED supplied the care setting at index.
2. **Sections 8.6 and 8.7 (penalised predictor selection, bounded machine learning) and TRIPOD reporting** were
   optional extensions for a prediction score. No score was developed, so they were not carried out and TRIPOD does not
   apply.
3. **Section 8.10 (uncertainty by patient-level bootstrap or repeated split resampling).** Confidence intervals are Wald
   intervals from the fitted models, combined across the twenty imputed datasets by Rubin's rules, and delta-method
   intervals for quantities derived from the hidden Markov model; the bootstrap was not needed once the models were
   fitted on person-interval data with multiply imputed covariates.
4. **Specification-stage scripts and logs.** The plan says these would be archived beside it. They were superseded by
   the deterministic pipeline in this repository, which regenerates every reported value from the raw tables and is
   verified against the manuscript (`R/09_verify.R`, `expected/expected_results.csv`). They are held by the
   corresponding author and are not published because they were run against interim extractions and their outputs are
   not separable from patient-level intermediate files; they can be shown to editors or reviewers under the terms of the
   data use agreement.
5. **Laboratory-definition and care-setting sensitivity analyses** are reported as Supplementary Tables S9 and S10.
6. **Extraction determinism.** During construction of this pipeline the choice of the index study within an episode was
   found to depend on an unordered aggregate in the original queries; ties are now broken explicitly (highest grade,
   earliest study, lowest measurement identifier), and the same-day episode order, median laboratory values, most-recent
   laboratory values, first intervention and type of dysfunction are likewise deterministic. Cohort membership, event
   counts, transitions and the hidden Markov inputs are unchanged; covariate-based estimates moved in the last reported
   digit.
7. **Time-varying stenosis model.** The model with reaching moderate or severe calcification as a time-varying exposure
   is fitted with multiply imputed covariates on all 28,515 at-risk patients (an earlier complete-case fit is
   superseded).
