# 08_tables.R: main and supplementary tables as CSV files in outputs/tables/, built only from the results objects.
# Formatting follows the manuscript: hazard ratios as "HR (lower–upper)", percentages as in the tables.
source("R/00_common.R")
r1 <- readRDS("outputs/results/results_part1.rds"); r2 <- readRDS("outputs/results/results_part2.rds"); sm <- readRDS("outputs/results/secondary_mi.rds")
TB <- "outputs/tables/"; dir.create(TB, showWarnings=FALSE, recursive=TRUE)
f  <- function(h, l, u, d=2) sprintf(paste0("%.",d,"f (%.",d,"f–%.",d,"f)"), h, l, u)
fp <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
lab <- c(age10="Age, per 10 years", male="Male sex", av_catmild="Aortic valve calcification, mild vs none", av_catmod="Aortic valve calcification, moderate vs none",
         av_catsev="Aortic valve calcification, severe vs none", ef_catlt40="LVEF <40% vs ≥50%", ef_cat40to49="LVEF 40–49% vs ≥50%", zE="E/e', per SD",
         zla="Left atrial dimension, per SD", zivs="Septal wall thickness, per SD", zbmi="Body mass index, per SD", af="Atrial fibrillation or flutter",
         eg30to60="eGFR 30–59 vs ≥60", eglt30="eGFR <30 vs ≥60", esrd="Dialysis dependence", zphos="Serum phosphate, per SD", zca="Serum calcium, per SD",
         zalp="Alkaline phosphatase, per SD (log)", dm="Diabetes", htn="Hypertension", cad="Coronary artery disease")
ord <- names(lab)
get <- function(df, st, de, tm, cols=c("hr","lo","hi")) { z <- df[df$stratum==st & df$def==de & df$term==tm, ]; if (nrow(z)) f(z[[cols[1]]], z[[cols[2]]], z[[cols[3]]]) else "" }
four_col <- function(df, terms=ord, cols=c("hr","lo","hi")) do.call(rbind, lapply(intersect(terms, unique(df$term)), function(t) data.frame(characteristic=unname(lab[t]),
   onset_first=get(df,"onset","first",t,cols), onset_confirmed=get(df,"onset","confirmed",t,cols), progression_first=get(df,"progression","first",t,cols), progression_confirmed=get(df,"progression","confirmed",t,cols))))
counts_of <- function(df) unique(df[, intersect(c("stratum","def","n_pt","n","events","deaths"), names(df))])
wcsv <- function(x, file) write.csv(x, paste0(TB, file), row.names=FALSE)

# ---------- Table 1 ----------
wcsv(r2$table1, "table1_baseline.csv"); writeLines(r2$table1_notes, paste0(TB, "table1_notes.txt"))
# ---------- Table 2A: transitions between consecutive episodes; 2B: fate of the first moderate or severe read ----------
tr <- r2$transitions; trp <- r2$transitions_pct
t2a <- data.frame(grade=c("None (blank)","Mild","Moderate","Severe"), consecutive_pairs=rowSums(tr), next_none=sprintf("%.1f",trp[,1]), next_mild=sprintf("%.1f",trp[,2]), next_moderate=sprintf("%.1f",trp[,3]), next_severe=sprintf("%.1f",trp[,4]))
wcsv(t2a, "table2a_transitions.csv")
extra <- data.frame(sequence=c("mild newly appeared (previous episode blank): next episode","two successive moderate or severe reads: next episode"),
                    n=c(r2$new_mild_next["n"], r2$repeat_mod_next["n"]), next_blank_pct=100*c(r2$new_mild_next["blank"], r2$repeat_mod_next["blank"]),
                    next_mild_pct=100*c(r2$new_mild_next["mild"], r2$repeat_mod_next["mild"]), next_modplus_pct=100*c(r2$new_mild_next["modplus"], r2$repeat_mod_next["modplus"]))
wcsv(extra, "table2a_sequences.csv")
fa <- r2$first_read_fate; a <- fa[fa$group=="all",]
t2b <- data.frame(cohort=ifelse(a$stratum=="onset","No reported calcification at index","Mild calcification at index"), reads=a$n, confirmed_pct=round(100*a$confirmed),
                  refuted_pct=round(100*a$refuted), never_reexamined_pct=round(100*a$last_echo), died_1y_if_never_reexamined_pct=round(100*a$died1y_last), died_1y_if_reexamined_pct=round(100*a$died1y_other))
wcsv(t2b, "table2b_fate.csv")
# ---------- Table 3: main models (D3) with the onset-versus-progression interaction ----------
m3 <- r1$main[r1$main$domain=="D3",]; it <- r1$interaction; it$term2 <- sub("^prog:","", sub(":prog$","", it$term)); itf <- it[it$def=="first",]
imap <- c(age10="age10", male="male", av_catmild="av_lin", av_catmod=NA, av_catsev=NA, ef_catlt40="ef_catlt40", ef_cat40to49="ef_cat40to49", zE="zE", zla="zla", zivs="zivs", zbmi="zbmi", af="af", eg30to60="eg30to60", eglt30="eglt30", esrd="esrd")
t3 <- do.call(rbind, lapply(names(imap), function(tm) { z <- if (is.na(imap[tm])) itf[0,] else itf[itf$term2==imap[tm],]   # moderate/severe AVC rows carry no ratio (the interaction is per grade)
  data.frame(characteristic=unname(lab[tm]), onset_first=get(m3,"onset","first",tm), onset_confirmed=get(m3,"onset","confirmed",tm), progression_first=get(m3,"progression","first",tm),
             progression_confirmed=get(m3,"progression","confirmed",tm), hr_ratio_progression_vs_onset=if (nrow(z)) f(z$hr,z$lo,z$hi) else "", p_interaction=if (nrow(z)) fp(z$p) else "") }))
t3$note <- ifelse(t3$characteristic==lab["av_catmild"], "HR ratio for aortic valve calcification is per grade", "")
wcsv(t3, "table3_main_models.csv"); wcsv(unique(m3[,c("stratum","def","n_pt","events")]), "table3_counts.csv")
# ---------- Table 4 ----------
inc <- r2$ms_incidence; tvc <- sm$ms_tvc[sm$ms_tvc$term=="mac",]; cif <- r2$ms_cif
cif_at <- function(g, yr) { rn <- grep(paste0("^", g, " 1$"), rownames(cif$est)); e <- cif$est[rn, as.character(yr)]; v <- cif$var[rn, as.character(yr)]; sprintf("%.1f (%.1f–%.1f)", 100*e, 100*pmax(e-1.96*sqrt(v),0), 100*(e+1.96*sqrt(v))) }
gname <- c("No reported calcification","Mild","Moderate or severe")
t4a <- data.frame(baseline_grade=c("None reported","Mild","Moderate or severe"), at_risk=inc$n, stenosis_events=inc$events, person_years=round(inc$py), rate_per_100py=sprintf("%.2f", inc$rate100), deaths_before_stenosis=inc$deaths,
                  cumulative_incidence_2y_pct=sapply(gname, cif_at, yr=2), cumulative_incidence_5y_pct=sapply(gname, cif_at, yr=5))
wcsv(t4a, "table4a_stenosis_incidence.csv")
writeLines(c(sprintf("Reaching moderate or severe MAC during follow-up (time-varying exposure, patients with no or mild MAC at index): HR %s; n = %d, events = %d, multiply imputed covariates.", f(tvc$hr,tvc$lo,tvc$hi), sm$ms_tvc_n["n"], sm$ms_tvc_n["events"]),
             sprintf("Follow-up ends at the last echocardiogram; death within %d days after it is a competing event and later death is censoring at the last echocardiogram (Supplementary Table S14 varies the window).", GRACE_DAYS)),
           paste0(TB, "table4a_time_varying_note.txt"))
iv <- r2$interv; bt <- iv$by_type; ic <- r2$interv_confirmed; ci_iv <- r2$interv_cif; mort <- r2$mortality_after_dysfunction
t4b <- data.frame(step=c("Patients with moderate or severe MAC at index or follow-up","Moderate or greater mitral dysfunction at or after the qualifying read","  Already documented at the qualifying episode","  First documented on a later study",
                          "  Regurgitation-led","  Stenosis-led or raised gradient","Without a prior mitral procedure","Mitral intervention after the dysfunction was documented","  After stenosis-led dysfunction","  After regurgitation-led dysfunction",
                          "  Isolated (no concomitant bypass or aortic valve surgery)","Died after the echocardiogram showing dysfunction (within the ascertainment window)"),
                  n=c(iv$n_modplus, iv$n_pheno, iv$n_pheno_same_episode, iv$n_pheno_later, iv$pheno_mr, iv$pheno_sten, iv$n_analysed, iv$n_interv, bt["sten_iv"], bt["mr_iv"], iv$isolated, iv$died),
                  denominator=c(NA, iv$n_modplus, iv$n_pheno, iv$n_pheno, iv$n_pheno, iv$n_pheno, iv$n_pheno, iv$n_analysed, bt["sten_n"], bt["mr_n"], iv$n_interv, iv$n_analysed))
t4b$pct <- ifelse(is.na(t4b$denominator), NA, round(100*t4b$n/t4b$denominator, 1))
wcsv(t4b, "table4b_cascade.csv")
cif_iv <- function(yr, code=1) { e <- ci_iv$est[code, as.character(yr)]; v <- ci_iv$var[code, as.character(yr)]; sprintf("%.1f%% (%.1f–%.1f)", 100*e, 100*pmax(e-1.96*sqrt(v),0), 100*(e+1.96*sqrt(v))) }
writeLines(c(sprintf("Procedure types among the %d interventions: %s.", iv$n_interv, paste(sprintf("%s %d", names(iv$modality), as.integer(iv$modality)), collapse="; ")),
             sprintf("Median days from documented dysfunction to intervention %s; median days to death among those who died %s.", format(iv$med_days), format(iv$med_days_death)),
             sprintf("Cumulative incidence of intervention with death as a competing event: %s at 1 year, %s at 2 years, %s at 5 years (n = %d; %d interventions; %d died before any intervention).", cif_iv(1), cif_iv(2), cif_iv(5), ci_iv$n, ci_iv$interventions, ci_iv$deaths_before_intervention),
             sprintf("Cumulative mortality after the echocardiogram showing dysfunction: %s at 1 year, %s at 3 years, %s at 5 years (Kaplan-Meier; median follow-up %.1f years by reverse Kaplan-Meier; vital status ascertained for one year after the last echocardiogram).", paste(sprintf("%.1f%% (%.1f–%.1f)", 100*mort$cum_mortality, 100*mort$lo, 100*mort$hi)[1]), paste(sprintf("%.1f%% (%.1f–%.1f)", 100*mort$cum_mortality, 100*mort$lo, 100*mort$hi)[2]), paste(sprintf("%.1f%% (%.1f–%.1f)", 100*mort$cum_mortality, 100*mort$lo, 100*mort$hi)[3]), r2$followup_after_dysfunction_median_years),
             sprintf("Confirmed definition, counted from the confirming echocardiogram: %d patients (%d moderate or severe at index plus confirmed events); %d with dysfunction at or after it; %d without a prior procedure; %d interventions; %d died.", ic["n_modplus"], ic["n_index"], ic["n_pheno"], ic["n_analysed"], ic["n_interv"], ic["died"])),
           paste0(TB, "table4b_notes.txt"))
po <- sm$pooled; olab <- c(age10p="Age, per 10 years", male="Male sex", ef_catlt40="LVEF <40% vs ≥50%", ef_cat40to49="LVEF 40–49% vs ≥50%", eg30to60="eGFR 30–59 vs ≥60", eglt30="eGFR <30 vs ≥60", esrd="Dialysis dependence", sten="Stenosis-led rather than regurgitation-led dysfunction", af="Atrial fibrillation or flutter")
t4c <- data.frame(characteristic=unname(olab), hr_mi=sapply(names(olab), function(k){ z <- po[po$term==k,]; if (nrow(z)) f(z$hr_mi, z$lo_mi, z$hi_mi) else "" }),
                  hr_complete_case=sapply(names(olab), function(k){ z <- po[po$term==k,]; if (nrow(z)) f(z$hr_cc, z$lo_cc, z$hi_cc) else "" }))
wcsv(t4c, "table4c_intervention_hr.csv")
writeLines(c(sprintf("Cause-specific hazard ratios for receipt of mitral intervention from the echocardiogram documenting dysfunction (death a competing event, treated as censoring for the cause-specific hazard). Multiply imputed: n = %d, interventions = %d; complete case: n = %d, interventions = %d.", sm$n, sm$events, sm$n_cc, sm$events_cc),
             sprintf("Survival after the echocardiogram showing dysfunction, intervention as a time-varying exposure (MI-pooled): HR %s; confounded by indication and treatment selection.", f(sm$surv_hr["hr"], sm$surv_hr["lo"], sm$surv_hr["hi"]))), paste0(TB, "table4c_notes.txt"))
# ---------- Supplementary tables ----------
gl <- c(age_cat="Age at the qualifying echocardiogram, years", ef_cat="LVEF", esrd="Dialysis", male="Sex"); ll <- c(lt40="<40%", `40to49`="40–49%", ge50="≥50%", `0`="No", `1`="Yes")
s1 <- fa[fa$group!="all",]; s1$level2 <- ifelse(s1$group=="male", ifelse(s1$level=="1","Men","Women"), ifelse(s1$level %in% names(ll), ll[s1$level], s1$level))
S1 <- data.frame(baseline_group=ifelse(s1$stratum=="onset","No reported calcification","Mild"), subgroup=gl[s1$group], level=s1$level2, reads=s1$n, confirmed_pct=round(100*s1$confirmed), refuted_pct=round(100*s1$refuted), never_reexamined_pct=round(100*s1$last_echo))
wcsv(S1, "tableS1_fate_by_subgroup.csv")
E <- r2$hmm$emission; L <- r2$hmm$emission_L; U <- r2$hmm$emission_U; ip <- r2$hmm$initprobs; sj <- r2$hmm$sojourn; g4 <- c("None","Mild","Moderate","Severe")
S2 <- data.frame(underlying_grade=g4, reported_none=f(100*E[1:4,1],100*L[1:4,1],100*U[1:4,1],1), reported_mild=f(100*E[1:4,2],100*L[1:4,2],100*U[1:4,2],1), reported_moderate=f(100*E[1:4,3],100*L[1:4,3],100*U[1:4,3],1),
                 reported_severe=f(100*E[1:4,4],100*L[1:4,4],100*U[1:4,4],1), share_at_index_pct=sprintf("%.1f",100*ip[1:4]), mean_years_in_grade=f(sj[1:4,"estimates"], sj[1:4,"L"], sj[1:4,"U"],1))
wcsv(S2, "tableS2_hidden_markov.csv")
dv <- r2$hmm$derived; writeLines(c(sprintf("Fitted to %d patients and %d observations; -2 log-likelihood %.1f; best of %d starting-value sets (Supplementary Table S15).", r2$hmm$n_pt, r2$hmm$n_obs, r2$hmm$minus2LL, nrow(r2$hmm$starts)), paste(names(dv), sprintf("%.4f", dv), sep=" = ")), paste0(TB, "tableS2_derived.txt"))
br <- r1$bridge; S3 <- data.frame(characteristic=unname(lab[ord[ord %in% br$term]]), onset=sapply(ord[ord %in% br$term], function(t){ z <- br[br$stratum=="onset" & br$term==t,]; if (nrow(z)) f(z$hr,z$lo,z$hi) else "" }),
                                  progression=sapply(ord[ord %in% br$term], function(t){ z <- br[br$stratum=="progression" & br$term==t,]; if (nrow(z)) f(z$hr,z$lo,z$hi) else "" }))
wcsv(S3, "tableS3_confirmation_imputed.csv")
cb <- r1$confirmed_baseline; S4 <- data.frame(characteristic=unname(lab[ord[ord %in% cb$term]]), first_observed=sapply(ord[ord %in% cb$term], function(t){ z <- cb[cb$def=="first" & cb$term==t,]; if (nrow(z)) f(z$hr,z$lo,z$hi) else "" }),
                                             confirmed=sapply(ord[ord %in% cb$term], function(t){ z <- cb[cb$def=="confirmed" & cb$term==t,]; if (nrow(z)) f(z$hr,z$lo,z$hi) else "" }))
wcsv(S4, "tableS4_stricter_onset_cohort.csv"); wcsv(unique(cb[,c("def","n_pt","events")]), "tableS4_counts.csv")
td <- r1$total_direct[r1$total_direct$def=="first",]; S5 <- do.call(rbind, lapply(c("zbmi","af","dm","htn","cad"), function(t) data.frame(characteristic=unname(lab[t]),
   onset_demographic_adjusted=f(td$total_hr[td$stratum=="onset"&td$term==t], td$total_lo[td$stratum=="onset"&td$term==t], td$total_hi[td$stratum=="onset"&td$term==t]),
   onset_fully_adjusted=f(td$direct_hr[td$stratum=="onset"&td$term==t], td$direct_lo[td$stratum=="onset"&td$term==t], td$direct_hi[td$stratum=="onset"&td$term==t]),
   progression_demographic_adjusted=f(td$total_hr[td$stratum=="progression"&td$term==t], td$total_lo[td$stratum=="progression"&td$term==t], td$total_hi[td$stratum=="progression"&td$term==t]),
   progression_fully_adjusted=f(td$direct_hr[td$stratum=="progression"&td$term==t], td$direct_lo[td$stratum=="progression"&td$term==t], td$direct_hi[td$stratum=="progression"&td$term==t]))))
wcsv(S5, "tableS5_total_vs_direct.csv")
doms <- c(D1="Demographic", D2="+ Cardiac structure", D3="+ Kidney function (main model)", D4="+ Mineral markers", D5="+ Coded diagnoses", D3_cox_companion="Cox companion")
for (st in c("onset","progression")) for (de in c("first","confirmed")) {
  m <- r1$main[r1$main$stratum==st & r1$main$def==de,]; hdr <- sapply(names(doms), function(dm){ z <- m[m$domain==dm,]; sprintf("%s (n = %s; events %d)", doms[dm], format(z$n_pt[1], big.mark=","), z$events[1]) })
  body <- do.call(rbind, lapply(ord, function(t){ cells <- sapply(names(doms), function(dm){ z <- m[m$domain==dm & m$term==t,]; if (nrow(z)) f(z$hr,z$lo,z$hi) else "" }); if (any(cells!="")) data.frame(characteristic=unname(lab[t]), t(cells)) else NULL }))
  names(body) <- c("characteristic", hdr); wcsv(body, sprintf("tableS6%s_stagewise_%s_%s.csv", c(onset=c(first="a",confirmed="b"), progression=c(first="c",confirmed="d"))[[paste0(st,".",de)]], st, de)) }
sn <- c(IIW="Inverse-intensity weighted (visit-process model)", IPCW="Censoring weighted"); se <- r1$sens
for (st in c("onset","progression")) { body <- do.call(rbind, lapply(ord, function(t){ cells <- unlist(lapply(c("first","confirmed"), function(de) sapply(names(sn), function(s){ z <- se[se$stratum==st & se$def==de & se$sens==s & se$term==t,]; if (nrow(z)) f(z$hr,z$lo,z$hi) else "" })))
   if (any(cells!="")) data.frame(characteristic=unname(lab[t]), t(cells)) else NULL })); names(body) <- c("characteristic", paste(rep(sn,2), rep(c("first-observed","confirmed"), each=2), sep=", "))
   wcsv(body, sprintf("tableS7%s_informative_observation_%s.csv", ifelse(st=="onset","a","b"), st)) }
wd <- r1$weight_diag; wd[,c("min","p5","p25","median","p75","p95","max","pct_truncated")] <- lapply(wd[,c("min","p5","p25","median","p75","p95","max","pct_truncated")], function(v) signif(v,3)); wcsv(wd, "tableS7c_weight_diagnostics.csv")
fg <- r1$finegray; S8 <- four_col(fg, cols=c("shr","lo","hi")); wcsv(S8, "tableS8_fine_gray.csv"); wcsv(unique(fg[,c("stratum","def","n","events","deaths")]), "tableS8_counts.csv")
sl <- r1$sens_lab; S9 <- do.call(rbind, lapply(unique(sl$definition), function(dfn) do.call(rbind, lapply(c("eg30to60","eglt30","esrd"), function(t) data.frame(definition=dfn, characteristic=unname(lab[t]),
   onset_first=get(sl[sl$definition==dfn,],"onset","first",t), onset_confirmed=get(sl[sl$definition==dfn,],"onset","confirmed",t), progression_first=get(sl[sl$definition==dfn,],"progression","first",t), progression_confirmed=get(sl[sl$definition==dfn,],"progression","confirmed",t))))))
sp <- r1$sens_phos; S9p <- do.call(rbind, lapply(unique(sp$definition), function(dfn) data.frame(definition=paste("phosphate:", dfn), characteristic=unname(lab["zphos"]),
   onset_first=get(sp[sp$definition==dfn,],"onset","first","zphos"), onset_confirmed=get(sp[sp$definition==dfn,],"onset","confirmed","zphos"), progression_first=get(sp[sp$definition==dfn,],"progression","first","zphos"), progression_confirmed=get(sp[sp$definition==dfn,],"progression","confirmed","zphos"))))
wcsv(rbind(S9, S9p), "tableS9_laboratory_definitions.csv"); wcsv(unique(rbind(sl[,c("definition","stratum","def","n_pt","events")], sp[,c("definition","stratum","def","n_pt","events")])), "tableS9_counts.csv")
ss <- r1$sens_setting
for (st in c("onset","progression")) { body <- do.call(rbind, lapply(ord, function(t){ cells <- unlist(lapply(c("outpatient index","inpatient index"), function(g) sapply(c("first","confirmed"), function(de){ z <- ss[ss$stratum==st & ss$setting==g & ss$def==de & ss$term==t,]; if (nrow(z)) f(z$hr,z$lo,z$hi) else "" })))
   if (any(cells!="")) data.frame(characteristic=unname(lab[t]), t(cells)) else NULL })); names(body) <- c("characteristic", paste(rep(c("outpatient index","inpatient index"), each=2), rep(c("first-observed","confirmed"),2), sep=", "))
   wcsv(body, sprintf("tableS10%s_care_setting_%s.csv", ifelse(st=="onset","a","b"), st)) }
wcsv(unique(ss[,c("setting","stratum","def","n_pt","events")]), "tableS10_counts.csv")
# ---- new supplementary tables ----
wcsv(four_col(r1$cloglog), "tableS11_cloglog_companion.csv"); wcsv(counts_of(r1$cloglog), "tableS11_counts.csv")                        # S11: exact interval-censored (complementary log-log) companion of the main model
wcsv(four_col(r1$complete_case), "tableS12_complete_case.csv"); wcsv(counts_of(r1$complete_case), "tableS12_counts.csv")                # S12: complete-case main model
wcsv(four_col(r1$landmark), "tableS13a_landmark_second_episode.csv"); wcsv(counts_of(r1$landmark), "tableS13a_counts.csv")               # S13a: follow-up from the second episode
wcsv(four_col(r1$rheum_excluded), "tableS13b_rheumatic_after_index_excluded.csv"); wcsv(counts_of(r1$rheum_excluded), "tableS13b_counts.csv") # S13b: rheumatic-after-index patients excluded
gr <- r2$ms_incidence_grace; gr$group <- c(blank="None reported", mild="Mild", modsev="Moderate or severe")[gr$group]; gr$rate_per_100py <- sprintf("%.2f", gr$rate100); gr$person_years <- round(gr$py)
S14a <- gr[, c("grace_days","group","n","events","person_years","rate_per_100py","deaths")]; wcsv(S14a, "tableS14a_death_window_stenosis_rates.csv")
S14b <- transform(sm$ms_tvc_grace, hr_ci=f(hr,lo,hi))[, c("grace_days","hr_ci","n","events")]; wcsv(S14b, "tableS14b_death_window_time_varying_hr.csv")   # S14: death-window sensitivity
hs <- r2$hmm$starts; hs$minus2LL <- round(hs$minus2LL, 1); wcsv(hs, "tableS15_hmm_starts.csv")                                             # S15: HMM starting values and convergence
if (!is.null(r2$hmm$prevalence)) { pv <- r2$hmm$prevalence; obs <- as.data.frame(round(pv$"Observed percentages",1)); exp <- as.data.frame(round(pv$"Expected percentages",1)); obs$kind <- "observed"; exp$kind <- "expected"; obs$years <- rownames(obs); exp$years <- rownames(exp); wcsv(rbind(obs, exp), "tableS15b_hmm_prevalence.csv") }
mort_txt <- setNames(sprintf("%.1f%% (%.1f–%.1f)", 100*mort$cum_mortality, 100*mort$lo, 100*mort$hi), mort$years)
S16 <- data.frame(horizon_years=c(1,2,3,5), cumulative_incidence_intervention=c(cif_iv(1), cif_iv(2), "", cif_iv(5)), cumulative_mortality=c(mort_txt["1"], "", mort_txt["3"], mort_txt["5"]))
wcsv(S16, "tableS16_intervention_and_mortality_horizons.csv"); wcsv(mort, "tableS16_mortality_km.csv")   # S16: fixed-horizon estimates
wcsv(it, "interaction_all_terms.csv")
writeLines(c(sprintf("Index episode: %d patients (%.1f%%) had more than one study in the index episode; span median %d days, 90th percentile %d, maximum %d; the defining study was later than the episode's first study for %d patients.", r2$index_episode["multi_study"], r2$index_episode["pct_multi"], r2$index_episode["span_median"], r2$index_episode["span_p90"], r2$index_episode["span_max"], r2$index_episode["def_not_first"]),
             sprintf("Model-independent check landmarked at the classification echocardiogram: HR %s (n = %d, events = %d); the earlier construction with the clock at the first read gave %s (n = %d, events = %d).", f(r2$internal_check["hr"], r2$internal_check["lo"], r2$internal_check["hi"]), r2$internal_check["n"], r2$internal_check["events"], f(r2$internal_check_from_first_read["hr"], r2$internal_check_from_first_read["lo"], r2$internal_check_from_first_read["hi"]), r2$internal_check_from_first_read["n"], r2$internal_check_from_first_read["events"])),
           paste0(TB, "diagnostics_notes.txt"))
cat("tables written to", TB, "\n")
