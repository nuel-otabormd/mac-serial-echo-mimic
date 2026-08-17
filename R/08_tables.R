# 08_tables.R: main and supplementary tables as CSV files in outputs/tables/, built only from the results objects.
# Formatting follows the manuscript: hazard ratios as "HR (lower\u2013upper)", percentages as in the tables.
if (!grepl("UTF-8", Sys.getlocale("LC_CTYPE"), ignore.case=TRUE)) for (loc in c("en_US.UTF-8","C.UTF-8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break   # write en dashes and other symbols as UTF-8 even when the calling shell uses the C locale
r1 <- readRDS("outputs/results/results_part1.rds"); r2 <- readRDS("outputs/results/results_part2.rds"); sm <- readRDS("outputs/results/secondary_mi.rds")
TB <- "outputs/tables/"; dir.create(TB, showWarnings=FALSE, recursive=TRUE)
f  <- function(h, l, u, d=2) sprintf(paste0("%.",d,"f (%.",d,"f\u2013%.",d,"f)"), h, l, u)
fp <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
lab <- c(age10="Age, per 10 years", male="Male sex", av_catmild="Aortic valve calcification, mild vs none", av_catmod="Aortic valve calcification, moderate vs none",
         av_catsev="Aortic valve calcification, severe vs none", ef_catlt40="LVEF <40% vs \u226550%", ef_cat40to49="LVEF 40\u201349% vs \u226550%", zE="E/e', per SD",
         zla="Left atrial dimension, per SD", zivs="Septal wall thickness, per SD", zbmi="Body mass index, per SD", af="Atrial fibrillation or flutter",
         eg30to60="eGFR 30\u201359 vs \u226560", eglt30="eGFR <30 vs \u226560", esrd="Dialysis dependence", zphos="Serum phosphate, per SD", zca="Serum calcium, per SD",
         zalp="Alkaline phosphatase, per SD (log)", dm="Diabetes", htn="Hypertension", cad="Coronary artery disease")
ord <- names(lab)
get <- function(df, st, de, tm, cols=c("hr","lo","hi")) { z <- df[df$stratum==st & df$def==de & df$term==tm, ]; if (nrow(z)) f(z[[cols[1]]], z[[cols[2]]], z[[cols[3]]]) else "" }

# ---------- Table 1 ----------
write.csv(r2$table1, paste0(TB, "table1_baseline.csv"), row.names=FALSE); writeLines(r2$table1_notes, paste0(TB, "table1_notes.txt"))

# ---------- Table 2A: transitions between consecutive episodes; 2B: fate of the first moderate or greater read ----------
tr <- r2$transitions; trp <- r2$transitions_pct
t2a <- data.frame(grade=c("None (blank)","Mild","Moderate","Severe"), consecutive_pairs=rowSums(tr), next_none=sprintf("%.1f",trp[,1]), next_mild=sprintf("%.1f",trp[,2]),
                  next_moderate=sprintf("%.1f",trp[,3]), next_severe=sprintf("%.1f",trp[,4]))
write.csv(t2a, paste0(TB, "table2a_transitions.csv"), row.names=FALSE)
extra <- data.frame(sequence=c("mild newly appeared (previous episode blank): next episode","two successive moderate or greater reads: next episode"),
                    n=c(r2$new_mild_next["n"], r2$repeat_mod_next["n"]), next_blank_pct=100*c(r2$new_mild_next["blank"], r2$repeat_mod_next["blank"]),
                    next_mild_pct=100*c(r2$new_mild_next["mild"], r2$repeat_mod_next["mild"]), next_modplus_pct=100*c(r2$new_mild_next["modplus"], r2$repeat_mod_next["modplus"]))
write.csv(extra, paste0(TB, "table2a_sequences.csv"), row.names=FALSE)
fa <- r2$first_read_fate; a <- fa[fa$group=="all",]
t2b <- data.frame(cohort=ifelse(a$stratum=="onset","No reported calcification at index","Mild calcification at index"), reads=a$n, confirmed_pct=round(100*a$confirmed),
                  refuted_pct=round(100*a$refuted), never_reexamined_pct=round(100*a$last_echo), died_1y_if_never_reexamined_pct=round(100*a$died1y_last), died_1y_if_reexamined_pct=round(100*a$died1y_other))
write.csv(t2b, paste0(TB, "table2b_fate.csv"), row.names=FALSE)

# ---------- Table 3: main models (D3) with the onset-versus-progression interaction ----------
m3 <- r1$main[r1$main$domain=="D3",]; it <- r1$interaction; it$term2 <- sub("^prog:","", sub(":prog$","", it$term)); itf <- it[it$def=="first",]
imap <- c(age10="age10", male="male", av_catmild="av_lin", av_catmod=NA, av_catsev=NA, ef_catlt40="ef_catlt40", ef_cat40to49="ef_cat40to49", zE="zE", zla="zla", zivs="zivs", zbmi="zbmi", af="af", eg30to60="eg30to60", eglt30="eglt30", esrd="esrd")
t3 <- do.call(rbind, lapply(names(imap), function(tm) { z <- if (is.na(imap[tm])) itf[0,] else itf[itf$term2==imap[tm],]   # moderate/severe AVC rows carry no ratio (the interaction is per grade)
  data.frame(characteristic=unname(lab[tm]), onset_first=get(m3,"onset","first",tm), onset_confirmed=get(m3,"onset","confirmed",tm), progression_first=get(m3,"progression","first",tm),
             progression_confirmed=get(m3,"progression","confirmed",tm), hr_ratio_progression_vs_onset=if (nrow(z)) f(z$hr,z$lo,z$hi) else "", p_interaction=if (nrow(z)) fp(z$p) else "") }))
t3$note <- ifelse(t3$characteristic==lab["av_catmild"], "HR ratio for aortic valve calcification is per grade", "")
write.csv(t3, paste0(TB, "table3_main_models.csv"), row.names=FALSE)
write.csv(unique(m3[,c("stratum","def","n_pt","events")]), paste0(TB, "table3_counts.csv"), row.names=FALSE)

# ---------- Table 4 ----------
inc <- r2$ms_incidence; tvc <- sm$ms_tvc[sm$ms_tvc$term=="mac",]
t4a <- data.frame(baseline_grade=c("None reported","Mild","Moderate or severe"), at_risk=inc$n, stenosis_events=inc$events, person_years=round(inc$py), rate_per_100py=sprintf("%.2f", inc$rate100), deaths_before_stenosis=inc$deaths)
write.csv(t4a, paste0(TB, "table4a_stenosis_incidence.csv"), row.names=FALSE)
writeLines(sprintf("Reaching moderate or greater MAC during follow-up (time-varying exposure, patients with no or mild MAC at index): HR %s; n = %d, events = %d, multiply imputed covariates.", f(tvc$hr,tvc$lo,tvc$hi), sm$ms_tvc_n["n"], sm$ms_tvc_n["events"]), paste0(TB, "table4a_time_varying_note.txt"))
iv <- r2$interv; bt <- iv$by_type; ic <- r2$interv_confirmed
t4b <- data.frame(step=c("Patients with moderate or greater MAC at index or follow-up","Developed moderate or greater mitral dysfunction","Regurgitation-led","Stenosis-led or raised gradient",
                          "Without a prior mitral procedure","Mitral intervention after the dysfunction was documented","After stenosis-led dysfunction","After regurgitation-led dysfunction",
                          "Isolated (no concomitant bypass or aortic valve surgery)","Died after the echocardiogram showing dysfunction"),
                  n=c(iv$n_modplus, iv$n_pheno, iv$pheno_mr, iv$pheno_sten, iv$n_analysed, iv$n_interv, bt["sten_iv"], bt["mr_iv"], iv$isolated, iv$died),
                  denominator=c(NA, iv$n_modplus, iv$n_pheno, iv$n_pheno, iv$n_pheno, iv$n_analysed, bt["sten_n"], bt["mr_n"], iv$n_interv, iv$n_analysed))
t4b$pct <- ifelse(is.na(t4b$denominator), NA, round(100*t4b$n/t4b$denominator, 1))
write.csv(t4b, paste0(TB, "table4b_cascade.csv"), row.names=FALSE)
writeLines(c(sprintf("Procedure types among the %d interventions: %s.", iv$n_interv, paste(sprintf("%s %d", names(iv$modality), as.integer(iv$modality)), collapse="; ")),
             sprintf("Median days from documented dysfunction to intervention %d; median days to death %d.", iv$med_days, iv$med_days_death),
             sprintf("Confirmed definition: %d patients (%d at index plus confirmed events); %d developed dysfunction; %d without prior procedure; %d interventions; %d died.", ic["n_modplus"], ic["n_index"], ic["n_pheno"], ic["n_analysed"], ic["n_interv"], ic["died"])),
           paste0(TB, "table4b_notes.txt"))
po <- sm$pooled; olab <- c(age10p="Age, per 10 years", male="Male sex", ef_catlt40="LVEF <40% vs \u226550%", ef_cat40to49="LVEF 40\u201349% vs \u226550%", eg30to60="eGFR 30\u201359 vs \u226560", eglt30="eGFR <30 vs \u226560", esrd="Dialysis dependence", sten="Stenosis-led rather than regurgitation-led dysfunction", af="Atrial fibrillation or flutter")
po$key <- sub("^ef_cat", "ef_cat", sub("^eg", "eg", po$term)); po$key[po$term=="ef_catlt40"] <- "ef_catlt40"; po$key[po$term=="ef_cat40to49"] <- "ef_cat40to49"; po$key[po$term=="eg30to60"] <- "eg30to60"; po$key[po$term=="eglt30"] <- "eglt30"
t4c <- data.frame(characteristic=unname(olab[names(olab)]), or_mi=sapply(names(olab), function(k) { z <- po[po$term==k,]; if (nrow(z)) f(z$or_mi, z$lo_mi, z$hi_mi) else "" }),
                  or_complete_case=sapply(names(olab), function(k) { z <- po[po$term==k,]; if (nrow(z)) f(z$or_cc, z$lo_cc, z$hi_cc) else "" }))
write.csv(t4c, paste0(TB, "table4c_intervention_or.csv"), row.names=FALSE)
writeLines(c(sprintf("Multiply imputed: n = %d, interventions = %d; complete case: n = %d, interventions = %d.", sm$n, sm$events, sm$n_cc, sm$events_cc),
             sprintf("Survival after documented dysfunction, intervention as time-varying exposure (MI-pooled): HR %s.", f(sm$surv_hr["hr"], sm$surv_hr["lo"], sm$surv_hr["hi"]))), paste0(TB, "table4c_notes.txt"))

# ---------- Supplementary tables ----------
gl <- c(age_cat="Age, years", ef_cat="LVEF", esrd="Dialysis", male="Sex"); ll <- c(lt40="<40%", `40to49`="40\u201349%", ge50="\u226550%", `0`="No", `1`="Yes")
s1 <- fa[fa$group!="all",]; s1$level2 <- ifelse(s1$group=="male", ifelse(s1$level=="1","Men","Women"), ifelse(s1$level %in% names(ll), ll[s1$level], s1$level))
S1 <- data.frame(baseline_group=ifelse(s1$stratum=="onset","No reported calcification","Mild"), subgroup=gl[s1$group], level=s1$level2, reads=s1$n, confirmed_pct=round(100*s1$confirmed), refuted_pct=round(100*s1$refuted), never_reexamined_pct=round(100*s1$last_echo))
write.csv(S1, paste0(TB, "tableS1_fate_by_subgroup.csv"), row.names=FALSE)
E <- r2$hmm$emission; L <- r2$hmm$emission_L; U <- r2$hmm$emission_U; ip <- r2$hmm$initprobs; sj <- r2$hmm$sojourn; g4 <- c("None","Mild","Moderate","Severe")
S2 <- data.frame(underlying_grade=g4, reported_none=f(100*E[1:4,1],100*L[1:4,1],100*U[1:4,1],1), reported_mild=f(100*E[1:4,2],100*L[1:4,2],100*U[1:4,2],1), reported_moderate=f(100*E[1:4,3],100*L[1:4,3],100*U[1:4,3],1),
                 reported_severe=f(100*E[1:4,4],100*L[1:4,4],100*U[1:4,4],1), share_at_index_pct=sprintf("%.1f",100*ip[1:4]), mean_years_in_grade=f(sj[1:4,"estimates"], sj[1:4,"L"], sj[1:4,"U"],1))
write.csv(S2, paste0(TB, "tableS2_hidden_markov.csv"), row.names=FALSE)
dv <- r2$hmm$derived; writeLines(c(sprintf("Fitted to %d patients; -2 log-likelihood %.1f.", r2$hmm$n_pt, r2$hmm$minus2LL), paste(names(dv), sprintf("%.4f", dv), sep=" = ")), paste0(TB, "tableS2_derived.txt"))
br <- r1$bridge; S3 <- data.frame(characteristic=unname(lab[ord[ord %in% br$term]]), onset=sapply(ord[ord %in% br$term], function(t){ z <- br[br$stratum=="onset" & br$term==t,]; if (nrow(z)) f(z$hr,z$lo,z$hi) else "" }),
                                  progression=sapply(ord[ord %in% br$term], function(t){ z <- br[br$stratum=="progression" & br$term==t,]; if (nrow(z)) f(z$hr,z$lo,z$hi) else "" }))
write.csv(S3, paste0(TB, "tableS3_confirmation_imputed.csv"), row.names=FALSE)
cb <- r1$confirmed_baseline; S4 <- data.frame(characteristic=unname(lab[ord[ord %in% cb$term]]), first_observed=sapply(ord[ord %in% cb$term], function(t){ z <- cb[cb$def=="first" & cb$term==t,]; if (nrow(z)) f(z$hr,z$lo,z$hi) else "" }),
                                             confirmed=sapply(ord[ord %in% cb$term], function(t){ z <- cb[cb$def=="confirmed" & cb$term==t,]; if (nrow(z)) f(z$hr,z$lo,z$hi) else "" }))
write.csv(S4, paste0(TB, "tableS4_stricter_onset_cohort.csv"), row.names=FALSE); write.csv(unique(cb[,c("def","n_pt","events")]), paste0(TB, "tableS4_counts.csv"), row.names=FALSE)
td <- r1$total_direct[r1$total_direct$def=="first",]; S5 <- do.call(rbind, lapply(c("zbmi","af","dm","htn","cad"), function(t) data.frame(characteristic=unname(lab[t]),
   onset_demographic_adjusted=f(td$total_hr[td$stratum=="onset"&td$term==t], td$total_lo[td$stratum=="onset"&td$term==t], td$total_hi[td$stratum=="onset"&td$term==t]),
   onset_fully_adjusted=f(td$direct_hr[td$stratum=="onset"&td$term==t], td$direct_lo[td$stratum=="onset"&td$term==t], td$direct_hi[td$stratum=="onset"&td$term==t]),
   progression_demographic_adjusted=f(td$total_hr[td$stratum=="progression"&td$term==t], td$total_lo[td$stratum=="progression"&td$term==t], td$total_hi[td$stratum=="progression"&td$term==t]),
   progression_fully_adjusted=f(td$direct_hr[td$stratum=="progression"&td$term==t], td$direct_lo[td$stratum=="progression"&td$term==t], td$direct_hi[td$stratum=="progression"&td$term==t]))))
write.csv(S5, paste0(TB, "tableS5_total_vs_direct.csv"), row.names=FALSE)
doms <- c(D1="Demographic", D2="+ Cardiac structure", D3="+ Kidney function (main model)", D4="+ Mineral markers", D5="+ Coded diagnoses", D3_cox_companion="Cox companion")
for (st in c("onset","progression")) for (de in c("first","confirmed")) {
  m <- r1$main[r1$main$stratum==st & r1$main$def==de,]; hdr <- sapply(names(doms), function(dm){ z <- m[m$domain==dm,]; sprintf("%s (n = %s; events %d)", doms[dm], format(z$n_pt[1], big.mark=","), z$events[1]) })
  body <- do.call(rbind, lapply(ord, function(t){ cells <- sapply(names(doms), function(dm){ z <- m[m$domain==dm & m$term==t,]; if (nrow(z)) f(z$hr,z$lo,z$hi) else "" }); if (any(cells!="")) data.frame(characteristic=unname(lab[t]), t(cells)) else NULL }))
  names(body) <- c("characteristic", hdr); write.csv(body, paste0(TB, sprintf("tableS6%s_stagewise_%s_%s.csv", c(onset=c(first="a",confirmed="b"), progression=c(first="c",confirmed="d"))[[paste0(st,".",de)]], st, de)), row.names=FALSE) }
sn <- c(echo_rate="Echo-rate adjusted", IIW="Inverse-intensity weighted", IPCW="Censoring weighted"); se <- r1$sens
for (st in c("onset","progression")) { body <- do.call(rbind, lapply(ord, function(t){ cells <- unlist(lapply(c("first","confirmed"), function(de) sapply(names(sn), function(s){ z <- se[se$stratum==st & se$def==de & se$sens==s & se$term==t,]; if (nrow(z)) f(z$hr,z$lo,z$hi) else "" })))
   if (any(cells!="")) data.frame(characteristic=unname(lab[t]), t(cells)) else NULL })); names(body) <- c("characteristic", paste(rep(sn,2), rep(c("first-observed","confirmed"), each=3), sep=", "))
   write.csv(body, paste0(TB, sprintf("tableS7%s_informative_observation_%s.csv", ifelse(st=="onset","a","b"), st)), row.names=FALSE) }
fg <- r1$finegray; S8 <- do.call(rbind, lapply(intersect(ord, unique(fg$term)), function(t) data.frame(characteristic=unname(lab[t]), onset_first=get(fg,"onset","first",t,c("shr","lo","hi")), onset_confirmed=get(fg,"onset","confirmed",t,c("shr","lo","hi")),
   progression_first=get(fg,"progression","first",t,c("shr","lo","hi")), progression_confirmed=get(fg,"progression","confirmed",t,c("shr","lo","hi")))))
write.csv(S8, paste0(TB, "tableS8_fine_gray.csv"), row.names=FALSE); write.csv(unique(fg[,c("stratum","def","n","events","deaths")]), paste0(TB, "tableS8_counts.csv"), row.names=FALSE)
sl <- r1$sens_lab; S9 <- do.call(rbind, lapply(unique(sl$definition), function(dfn) do.call(rbind, lapply(c("eg30to60","eglt30","esrd"), function(t) data.frame(definition=dfn, characteristic=unname(lab[t]),
   onset_first=get(sl[sl$definition==dfn,],"onset","first",t), onset_confirmed=get(sl[sl$definition==dfn,],"onset","confirmed",t), progression_first=get(sl[sl$definition==dfn,],"progression","first",t), progression_confirmed=get(sl[sl$definition==dfn,],"progression","confirmed",t))))))
sp <- r1$sens_phos; S9p <- do.call(rbind, lapply(unique(sp$definition), function(dfn) data.frame(definition=paste("phosphate:", dfn), characteristic=unname(lab["zphos"]),
   onset_first=get(sp[sp$definition==dfn,],"onset","first","zphos"), onset_confirmed=get(sp[sp$definition==dfn,],"onset","confirmed","zphos"), progression_first=get(sp[sp$definition==dfn,],"progression","first","zphos"), progression_confirmed=get(sp[sp$definition==dfn,],"progression","confirmed","zphos"))))
write.csv(rbind(S9, S9p), paste0(TB, "tableS9_laboratory_definitions.csv"), row.names=FALSE); write.csv(unique(rbind(sl[,c("definition","stratum","def","n_pt","events")], sp[,c("definition","stratum","def","n_pt","events")])), paste0(TB, "tableS9_counts.csv"), row.names=FALSE)
ss <- r1$sens_setting
for (st in c("onset","progression")) { body <- do.call(rbind, lapply(ord, function(t){ cells <- unlist(lapply(c("outpatient index","inpatient index"), function(g) sapply(c("first","confirmed"), function(de){ z <- ss[ss$stratum==st & ss$setting==g & ss$def==de & ss$term==t,]; if (nrow(z)) f(z$hr,z$lo,z$hi) else "" })))
   if (any(cells!="")) data.frame(characteristic=unname(lab[t]), t(cells)) else NULL })); names(body) <- c("characteristic", paste(rep(c("outpatient index","inpatient index"), each=2), rep(c("first-observed","confirmed"),2), sep=", "))
   write.csv(body, paste0(TB, sprintf("tableS10%s_care_setting_%s.csv", ifelse(st=="onset","a","b"), st)), row.names=FALSE) }
write.csv(unique(ss[,c("setting","stratum","def","n_pt","events")]), paste0(TB, "tableS10_counts.csv"), row.names=FALSE)
write.csv(it, paste0(TB, "interaction_all_terms.csv"), row.names=FALSE)
cat("tables written to", TB, "\n")
