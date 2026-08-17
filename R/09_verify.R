# 09_verify.R: compare the numbers produced by this run with the reference values reported in the manuscript
# (expected/expected_results.csv). Usage:
#   Rscript R/09_verify.R            compare current outputs with the manifest; exits with status 1 if any check fails
#   Rscript R/09_verify.R --write    write the manifest from the current outputs (used once, from the reference run)
args <- commandArgs(trailingOnly=TRUE); write_mode <- "--write" %in% args
r1 <- readRDS("outputs/results/results_part1.rds"); r2 <- readRDS("outputs/results/results_part2.rds"); sm <- readRDS("outputs/results/secondary_mi.rds")
m3 <- r1$main[r1$main$domain=="D3",]; it <- r1$interaction[r1$interaction$def=="first",]; fl <- r2$flow; fa <- r2$first_read_fate[r2$first_read_fate$group=="all",]
hr <- function(st, de, tm) m3$hr[m3$stratum==st & m3$def==de & m3$term==tm]; itr <- function(tm) it$hr[it$term==paste0("prog:",tm)]
E <- r2$hmm$emission; iv <- r2$interv; inc <- r2$ms_incidence
vals <- c(
  cohort_analysed=unname(fl["analysed"]), cohort_onset=unname(fl["blank"]), cohort_progression=unname(fl["mild"]), cohort_advanced=unname(fl["modsev"]),
  excluded_prosthesis=unname(fl["excl_pros"]), excluded_rheumatic=unname(fl["excl_rheum"]), eligible=unname(fl["eligible"]),
  events_onset_first=unique(m3$events[m3$stratum=="onset"&m3$def=="first"]), events_onset_confirmed=unique(m3$events[m3$stratum=="onset"&m3$def=="confirmed"]),
  events_progression_first=unique(m3$events[m3$stratum=="progression"&m3$def=="first"]), events_progression_confirmed=unique(m3$events[m3$stratum=="progression"&m3$def=="confirmed"]),
  stricter_cohort_n=unique(r1$confirmed_baseline$n_pt), stricter_events_first=unique(r1$confirmed_baseline$events[r1$confirmed_baseline$def=="first"]),
  hr_age_onset_first=hr("onset","first","age10"), hr_age_progression_first=hr("progression","first","age10"), hr_male_onset_first=hr("onset","first","male"),
  hr_avsev_onset_first=hr("onset","first","av_catsev"), hr_avsev_progression_first=hr("progression","first","av_catsev"), hr_avsev_onset_confirmed=hr("onset","confirmed","av_catsev"),
  hr_zE_onset_first=hr("onset","first","zE"), hr_eglt30_onset_first=hr("onset","first","eglt30"), hr_esrd_onset_first=hr("onset","first","esrd"), hr_eflt40_onset_confirmed=hr("onset","confirmed","ef_catlt40"),
  ratio_age=itr("age10"), ratio_av_per_grade=itr("av_lin"), ratio_eg30to60=itr("eg30to60"), ratio_zE=itr("zE"), ratio_zbmi=itr("zbmi"),
  hmm_n_pt=r2$hmm$n_pt, hmm_mild_reported_none=E[2,1], hmm_moderate_reported_modplus=E[3,3]+E[3,4], hmm_none_reported_none=E[1,1], hmm_ppv_modplus=unname(r2$hmm$derived["PPV_modplus_read.State 3"]),
  hmm_latent_mild_given_blank=unname(r2$hmm$derived["P_true_mild_given_blank.State 2"]), hmm_free_given_mild=unname(r2$hmm$derived["P_true_free_given_mild.State 1"]), hmm_sojourn_none_years=unname(r2$hmm$sojourn[1,"estimates"]),
  fate_onset_confirmed_pct=100*fa$confirmed[fa$stratum=="onset"], fate_onset_refuted_pct=100*fa$refuted[fa$stratum=="onset"], fate_onset_never_pct=100*fa$last_echo[fa$stratum=="onset"],
  fate_progression_confirmed_pct=100*fa$confirmed[fa$stratum=="progression"], fate_progression_never_pct=100*fa$last_echo[fa$stratum=="progression"],
  transition_none_to_none_pct=r2$transitions_pct[1,1], transition_mild_to_none_pct=r2$transitions_pct[2,1], repeat_modplus_next_modplus_pct=100*unname(r2$repeat_mod_next["modplus"]),
  stenosis_rate_none=inc$rate100[1], stenosis_rate_mild=inc$rate100[2], stenosis_rate_modsev=inc$rate100[3], stenosis_atrisk_none=inc$n[1], stenosis_atrisk_mild=inc$n[2], stenosis_atrisk_modsev=inc$n[3],
  stenosis_tvc_hr_mi=sm$ms_tvc$hr[sm$ms_tvc$term=="mac"], internal_check_hr=unname(r2$internal_check["hr"]),
  interv_modplus=iv$n_modplus, interv_dysfunction=iv$n_pheno, interv_analysed=iv$n_analysed, interv_n=iv$n_interv, interv_died=iv$died, interv_median_days=iv$med_days,
  interv_or_age_mi=sm$pooled$or_mi[sm$pooled$term=="age10p"], interv_surv_hr_mi=unname(sm$surv_hr["hr"]),
  confirmed_def_modplus=unname(r2$interv_confirmed["n_modplus"]), confirmed_def_interv=unname(r2$interv_confirmed["n_interv"]))
tol <- ifelse(grepl("^(cohort|excluded|eligible|events|stricter|hmm_n_pt|stenosis_atrisk|interv_(modplus|dysfunction|analysed|n|died|median)|confirmed_def)", names(vals)), 0,
              ifelse(grepl("pct$|rate|sojourn", names(vals)), 0.05, 0.005))
if (write_mode) { dir.create("expected", showWarnings=FALSE); write.csv(data.frame(key=names(vals), value=unname(vals), tolerance=tol), "expected/expected_results.csv", row.names=FALSE); cat("manifest written: expected/expected_results.csv (", length(vals), "values )\n"); quit(status=0) }
exp <- read.csv("expected/expected_results.csv"); ok <- TRUE
for (i in seq_len(nrow(exp))) { k <- exp$key[i]; got <- vals[[k]]; pass <- !is.null(got) && abs(got - exp$value[i]) <= exp$tolerance[i]
  cat(sprintf("%-38s expected %12.4f  got %12.4f  %s\n", k, exp$value[i], ifelse(is.null(got), NA, got), ifelse(pass, "PASS", "FAIL"))); ok <- ok && pass }
cat(ifelse(ok, "\nALL CHECKS PASSED\n", "\nSOME CHECKS FAILED\n")); quit(status=ifelse(ok, 0, 1))
