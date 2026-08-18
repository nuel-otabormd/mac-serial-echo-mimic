# Deviations from, and clarifications to, the archived analysis plan

`MAC_Study_Protocol_v2.4.md` is reproduced as written on 16 August 2026 and has not been edited. It was developed
iteratively against the data between 12 and 16 August 2026, and that specification work produced preliminary hazard ratios
under superseded specifications, as the protocol itself discloses. The analysis in this repository (versions 1.1.0 to 1.2.0, 17 to 18 August
2026) differs from the plan, and from the first archived implementation (version 1.0.x, same day), in the respects below.
Items 1 to 6 are clarifications; items 7 to 19 are changes made after an independent methodological review of the first
implementation, before submission. None changes an outcome definition, an exposure or a reporting rule; several change the
timing conventions and the estimator's implementation, and every reported number was regenerated after them.

## Clarifications of the plan

1. **Data sources.** The plan lists MIMIC-IV-ECHO v0.1 and MIMIC-IV-Note. The structured echocardiography measurements
   used come from MIMIC-IV-ECHO version 1.0 (the version identity of the PhysioNet resource was verified after the plan was
   written; v0.1 is a smaller subset that was never used). MIMIC-IV-Note was not used: mitral interventions were identified
   from ICD-9-CM and ICD-10-PCS procedure codes only. MIMIC-IV-ED supplied the care setting at index and emergency
   department stays for the prior-contact definition.
2. **Sections 8.6 and 8.7 (penalised predictor selection, bounded machine learning) and TRIPOD reporting** were optional
   extensions for a prediction score. No score was developed, so they were not carried out and TRIPOD does not apply.
3. **Section 8.10 (uncertainty by patient-level bootstrap or repeated split resampling).** Confidence intervals are Wald
   intervals with patient-clustered sandwich variances (item 9), combined across the twenty imputed datasets by Rubin's
   rules, and delta-method intervals for quantities derived from the hidden Markov model; the bootstrap was not used.
4. **Specification-stage scripts and logs.** The plan says these would be archived beside it. They were superseded by the
   deterministic pipeline in this repository, which regenerates every reported value from the raw tables and is verified
   against the manuscript (`R/09_verify.R`, `expected/expected_results.csv`). They are held by the corresponding author and
   are not published because they were run against interim extractions and their outputs are not separable from
   patient-level intermediate files; they can be shown to editors or reviewers under the terms of the data use agreement.
5. **Blank-as-censored hidden Markov sensitivity (plan section 8.4).** The run that treated blank reports as censored
   between the calcification-free and mild states did not converge to an identifiable solution (Hessian not positive
   definite, degenerate initial-state distribution) and is not reported.
6. **Laboratory-definition and care-setting sensitivity analyses** are reported as Supplementary Tables S9 and S10.

## Changes made after the independent methodological review

7. **One physical study defines each episode, and time zero is the end of the index episode.** The first implementation
   dated an episode at its earliest study but took the grade and the baseline covariates from the highest-grade study of the
   episode, so covariates could post-date time zero. Now the study with the highest grade (earliest when tied) supplies the
   episode's grade and, at index, the covariates, and time zero is the last study of the index episode, so that every
   baseline quantity is measured at or before time zero and follow-up starts after the last index study. Later episodes are
   dated by their defining study. Episodes are hospital stays (overlapping admissions merged) or calendar weeks outside a
   stay, and runs whose calendar days overlap are one episode (the first implementation could produce two same-day episodes
   when an emergency department study preceded an admission). Diagnostics: 6.0% of index episodes (1,854 patients) contain
   more than one study; among them the defining study is the first study in 87.4%, and the covariates come from a study
   before time zero in 1,589 patients (5.1% of the cohort; median 6 days earlier, 90th percentile 27, maximum 135).
8. **Exposure and outcomes on one timeline.** Secondary outcomes (calcific stenosis, regurgitation, mitral dysfunction) are
   ascertained only on studies after the last study of the index episode; dysfunction under the confirmed definition is
   counted from the echocardiogram that confirms the pair; dysfunction already documented on the qualifying episode is
   reported separately from dysfunction first documented later.
9. **Estimator implementation.** Person-intervals are split at the time-band cut points (0.5, 1, 2, 4, 7 years) so that
   the piecewise baseline hazard is what the plan describes; standard errors are patient-clustered; a grouped-time
   complementary log-log model fitted to the unsplit intervals between examinations is reported as a companion (Supplementary
   Table S11; version 1.1.0 fitted it to the split intervals with the event kept in the last segment, which is not the
   interval-censored likelihood when an interval crosses a band cut point, as about 60% of event intervals do; corrected in
   version 1.1.2); model n is taken from the fitted data (the mineral-marker models are fitted in patients with all three
   markers measured).
10. **Death.** Observation ends at the last echocardiogram; in analyses that need a death time, death within 30 days after
    the last echocardiogram is a competing event and later death is censoring at the last echocardiogram, because MIMIC-IV
    records deaths for one year after the last hospital discharge. The rule, previously undisclosed, is now stated and varied
    (0, 90, 365 days; Supplementary Table S14). Fine and Gray models are pooled across the twenty imputed datasets.
11. **Baseline covariates strictly at or before index.** Atrial fibrillation from the electrocardiogram now uses ECGs from
    30 days before index to the end of the index episode (previously 30 days either side of index). "Prior contact of one
    year" counts admissions, emergency department stays and echocardiograms (previously admissions only).
12. **Cohort definition uses only information available by the end of the index episode.** Prosthesis or ring and rheumatic
    morphology are judged on studies up to the end of the index episode. Patients whose only rheumatic report came later
    (n = 56) are kept in the primary cohort and excluded in a sensitivity analysis (Supplementary Table S13b); the first
    implementation excluded them.
13. **Observation-process weights.** The constant patient-level inverse-intensity weight and the echo-rate covariate of the
    first implementation, both of which used each patient's eventual echo count, are replaced by a visit-process model
    (Andersen and Gill Cox model for the next echocardiogram with the baseline covariates and time-updated history) that
    supplies stabilised weights per interval; weight diagnostics are reported (Supplementary Table S7c). Truncation of the
    censoring weights at 0.1 and 10 is now stated.
14. **Confirmed-versus-refuted comparison** is landmarked at the classification echocardiogram (follow-up starts when the
    classification becomes known); the earlier construction started the clock at the first read.
15. **Intervention** is analysed as a time-to-event outcome from the dysfunction echocardiogram with death as a competing
    event (cumulative incidence at fixed horizons and cause-specific hazard ratios), replacing the logistic regression;
    mortality after dysfunction is reported at fixed horizons.
16. **Additional sensitivity analyses added:** landmark at the second episode for both cohorts (Supplementary Table S13a);
    complete-case main models (Supplementary Table S12); hidden Markov model refitted from five dispersed starting sets with
    the table of starts and observed-versus-expected prevalences (Supplementary Table S15); the imputation repeated with a
    second random seed to show the simulation error of the pooled estimates (Supplementary Table S17; version 1.1.1).
17. **Source population.** The cohort flow now begins with all adults who had a qualifying study, showing those with a
    single episode; the estimand is stated as conditional on serial echocardiography.
18. **Event and diagnosis timing (version 1.1.2).** The confirmed event, dated at the first study of the pair in the plan and
    the main analysis, is also analysed dated at the confirming echocardiogram; the multiply imputed confirmation status is
    described as addressing missing confirmation under a missing-at-random assumption, not as removing selection.
    Coded diagnoses were taken from admissions begun by the end of the index episode, but MIMIC-IV assigns codes per admission
    at discharge without a date within the stay, so the models were repeated with codes only from admissions completed by time
    zero (electrocardiogram and dialysis records unchanged). Both are in Supplementary Table S11, which now holds the model and
    timing companions in one table.
19. **Episode-level valve assessment and vital-status horizon (version 1.2.0).** For the secondary outcomes, baseline mitral
    stenosis or regurgitation is now judged on every study of the index episode, not only its defining study, and mitral
    dysfunction on any study of the qualifying (first moderate or severe) episode counts as already present, whether that
    study preceded or followed the episode's defining study; both follow the episode-level logic used elsewhere. The confirmed
    definition is stated precisely: the first moderate or severe report that was reproduced on the next episode (a first
    report followed by a lower grade did not qualify but a later reproduced report did; 214 of 250 and 464 of 528 confirmed
    events were the first report itself), with a sensitivity analysis restricted to the first report. Survivors in the
    descriptive mortality summaries are censored one year after the last hospital discharge, the period for which MIMIC-IV
    records deaths, or at the last echocardiogram if later (previously one year after the last echocardiogram). Missingness
    of the imputed covariates is stated. The supplement shows the principal correlates and the complete tables ship as
    machine-readable files (`release_tables/`).
    The manuscript now calls the confirmed endpoint "reproduced" and the first-report dispositions "reproduced", "lower grade
    next" and "not re-examined"; the code and CSV files keep the internal names (confirmed, refuted, unconfirmable).
