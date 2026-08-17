# 05_reliability_intervention.R: baseline table, grade transitions, fate of first moderate-or-greater reads, hidden Markov
# derived quantities, calcific stenosis incidence and time-varying exposure model, model-independent check, intervention cascade.
suppressPackageStartupMessages({library(survival); library(msm)})
d <- readRDS("data/frame.rds"); p <- read.csv("data/panel.csv"); res2 <- list()
mk <- function(tp, td, tl){ ev <- ifelse(!is.na(tp),1, ifelse(!is.na(td) & (is.na(tl) | td<=tl+30),2,0)); tt <- ifelse(!is.na(tp),tp, ifelse(ev==2, td, tl)); list(ev=ev, tt=pmax(tt,1)/365.25) }
# ---------- TABLE 1 source ----------
d$egfr_base_cat <- cut(d$egfr_base, c(-Inf,30,60,Inf), labels=c("lt30","30to60","ge60"))
grp <- ifelse(d$sev1==0,"blank", ifelse(d$sev1==1,"mild","modsev")); d$grp <- factor(grp, levels=c("blank","mild","modsev"))
t1 <- list()
add <- function(lab, f) t1[[length(t1)+1]] <<- data.frame(row=lab, blank=f(d[d$grp=="blank",]), mild=f(d[d$grp=="mild",]), modsev=f(d[d$grp=="modsev",]), stringsAsFactors=FALSE)
msd <- function(x) sprintf("%.1f (%.1f)", mean(x,na.rm=T), sd(x,na.rm=T)); miqr <- function(x, dg=1) sprintf(paste0("%.",dg,"f (%.",dg,"f to %.",dg,"f)"), median(x,na.rm=T), quantile(x,.25,na.rm=T), quantile(x,.75,na.rm=T)); npc <- function(b) sprintf("%d (%.1f)", sum(b,na.rm=T), 100*mean(b,na.rm=T))
add("Patients, n", function(x) sprintf("%d", nrow(x)))
add("Age, years, mean (SD)", function(x) msd(x$age0)); add("Women, n (%)", function(x) npc(x$male==0))
add("Inpatient at index echocardiogram, n (%)", function(x) npc(x$idx_inpt==1)); add("  Intensive care, n (%)", function(x) npc(x$setting=="ICU"))
add("Aortic valve calcification, n (%)", function(x) ""); add("  None", function(x) npc(x$av==0)); add("  Mild", function(x) npc(x$av==1)); add("  Moderate", function(x) npc(x$av==2)); add("  Severe", function(x) npc(x$av==3)); add("  Not recorded", function(x) npc(is.na(x$av)))
add("Left ventricular ejection fraction, n (%)", function(x) ""); add("  <40%", function(x) npc(x$lvef<40)); add("  40 to 49%", function(x) npc(x$lvef>=40 & x$lvef<50)); add("  >=50%", function(x) npc(x$lvef>=50)); add("  Not recorded", function(x) npc(is.na(x$lvef)))
add("E/e' (septal), median (IQR)", function(x) miqr(x$E_sept)); add("Left atrial dimension, cm, mean (SD)", function(x) msd(x$la)); add("Septal wall thickness, cm, mean (SD)", function(x) msd(x$ivs)); add("Body mass index, kg/m2, mean (SD)", function(x) msd(x$bmi))
add("Baseline eGFR, mL/min/1.73 m2, n (%)", function(x) ""); add("  >=60", function(x) npc(x$egfr_base_cat=="ge60")); add("  30 to 59", function(x) npc(x$egfr_base_cat=="30to60")); add("  <30", function(x) npc(x$egfr_base_cat=="lt30")); add("  No creatinine in prior year", function(x) npc(is.na(x$egfr_base)))
add("Dialysis dependence, n (%)", function(x) npc(x$esrd==1)); add("Atrial fibrillation or flutter, n (%)", function(x) npc(x$af==1))
add("Serum phosphate, mg/dL, mean (SD)*", function(x) msd(x$phos_median)); add("Serum calcium, mg/dL, mean (SD)*", function(x) msd(x$ca_median)); add("Alkaline phosphatase, U/L, median (IQR)*", function(x) miqr(x$alp_median,0))
add("Diabetes, n (%)+", function(x) npc(x$dm[x$yr1==1]==1)); add("Hypertension, n (%)+", function(x) npc(x$htn[x$yr1==1]==1)); add("Coronary artery disease, n (%)+", function(x) npc(x$cad[x$yr1==1]==1))
add("Echocardiography episodes per patient, median (IQR)", function(x) miqr(x$n_ep,0)); add("Follow-up to last echocardiogram, years, median (IQR)", function(x) miqr(x$t_last/365.25))
add("Died during follow-up, n (%)", function(x) npc(!is.na(x$t_death)))
res2$table1 <- do.call(rbind, t1); res2$table1_notes <- sprintf("*Among patients with a value: phosphate %d, calcium %d, alkaline phosphatase %d. +Among the %d patients with at least one year of record before index.", sum(!is.na(d$phos_median)), sum(!is.na(d$ca_median)), sum(!is.na(d$alp_median)), sum(d$yr1==1))
# ---------- RELIABILITY ----------
tr <- table(factor(p$sev,levels=0:3), factor(p$nxt,levels=0:3)); res2$transitions <- tr; res2$transitions_pct <- round(100*prop.table(tr,1),1)
p <- p[order(p$pid,p$rn),]; p$prev <- ave(p$sev, p$pid, FUN=function(v) c(NA, head(v,-1)))
nm <- p[!is.na(p$prev) & p$prev==0 & p$sev==1 & !is.na(p$nxt),]; mm <- p[!is.na(p$prev) & p$prev>=2 & p$sev>=2 & !is.na(p$nxt),]
res2$new_mild_next <- c(n=nrow(nm), blank=mean(nm$nxt==0), mild=mean(nm$nxt==1), modplus=mean(nm$nxt>=2)); res2$repeat_mod_next <- c(n=nrow(mm), blank=mean(mm$nxt==0), mild=mean(mm$nxt==1), modplus=mean(mm$nxt>=2))
f <- d[d$sev1<2 & d$first_status %in% c("confirmed","refuted","unconfirmable"),]; f$died1y <- !is.na(f$t_death) & (f$t_death-f$t_first)<=365; f$age_cat <- cut(f$age0,c(0,65,75,85,200),labels=c("<65","65-74","75-84","85+")); f$ef_cat <- cut(f$lvef,c(-Inf,39.99,49.99,Inf),labels=c("lt40","40to49","ge50"))
fate <- list()
for (s in 0:1) { x <- f[f$sev1==s,]; fate[[length(fate)+1]] <- data.frame(stratum=ifelse(s==0,"onset","progression"), group="all", level="all", n=nrow(x), confirmed=mean(x$first_status=="confirmed"), refuted=mean(x$first_status=="refuted"), last_echo=mean(x$first_status=="unconfirmable"), died1y_last=mean(x$died1y[x$first_status=="unconfirmable"]), died1y_other=mean(x$died1y[x$first_status!="unconfirmable"]), med_days_death_last=median((x$t_death-x$t_first)[x$first_status=="unconfirmable" & !is.na(x$t_death)]))
  for (v in c("age_cat","ef_cat","esrd","male")) for (lv in levels(factor(x[[v]]))) { y <- x[!is.na(x[[v]]) & x[[v]]==lv,]; fate[[length(fate)+1]] <- data.frame(stratum=ifelse(s==0,"onset","progression"), group=v, level=lv, n=nrow(y), confirmed=mean(y$first_status=="confirmed"), refuted=mean(y$first_status=="refuted"), last_echo=mean(y$first_status=="unconfirmable"), died1y_last=NA, died1y_other=NA, med_days_death_last=NA) } }
res2$first_read_fate <- do.call(rbind, fate)
load("outputs/results/hmm.RData"); E <- ematrix.msm(mF); ip <- mF$hmodel$initprobs; ip <- if (is.matrix(ip)) ip[,1] else ip
res2$hmm <- list(emission=E$estimates[1:4,1:4], emission_L=E$L[1:4,1:4], emission_U=E$U[1:4,1:4], initprobs=ip[1:4], sojourn=sojourn.msm(mF), q=qmatrix.msm(mF)$estimates[1:4,], minus2LL=mF$minus2loglik, n_pt=length(unique(mF$data$mf$`(subject)`)))
Em <- E$estimates; pb <- ip[1:4]*Em[1:4,1]; pmi <- ip[1:4]*Em[1:4,2]; pm <- ip[1:4]*(Em[1:4,3]+Em[1:4,4])
res2$hmm$derived <- c(P_true_mild_given_blank=pb[2]/sum(pb), P_true_modsev_given_blank=(pb[3]+pb[4])/sum(pb), P_true_free_given_mild=pmi[1]/sum(pmi), P_true_modsev_given_mild=(pmi[3]+pmi[4])/sum(pmi), PPV_modplus_read=(pm[3]+pm[4])/sum(pm), sens_modplus_true_mod=Em[3,3]+Em[3,4], sens_modplus_true_sev=Em[4,3]+Em[4,4])
# ---------- CALCIFIC MS ENDPOINT + INTERNAL CHECK ----------
py <- function(tt, ev) c(events=sum(ev==1), py=sum(tt), rate100=100*sum(ev==1)/sum(tt))
ms <- list(); for (s in 0:2) { x <- d[if (s<2) d$sev1==s else d$sev1>=2,]; x <- x[x$ms0==0 & !is.na(x$ms0),]; z <- mk(x$t_ms,x$t_death,x$t_last); ms[[length(ms)+1]] <- data.frame(group=c("blank","mild","modsev")[s+1], n=nrow(x), t(py(z$tt,z$ev)), deaths=sum(z$ev==2)) }
res2$ms_incidence <- do.call(rbind, ms)
mr <- list(); for (s in 0:2) { x <- d[if (s<2) d$sev1==s else d$sev1>=2,]; x <- x[x$mr0==0 & !is.na(x$mr0),]; z <- mk(x$t_mr,x$t_death,x$t_last); mr[[length(mr)+1]] <- data.frame(group=c("blank","mild","modsev")[s+1], n=nrow(x), t(py(z$tt,z$ev))) }
res2$mr_incidence <- do.call(rbind, mr)
d$age10 <- (d$age0-70)/10; d$era_n <- as.numeric(factor(d$era))-3; d$ef_cat <- relevel(cut(d$lvef,c(-Inf,39.99,49.99,Inf),labels=c("lt40","40to49","ge50")),"ge50"); d$eg <- relevel(cut(d$egfr_base,c(-Inf,30,60,Inf),labels=c("lt30","30to60","ge60")),"ge60")
z <- function(v) as.numeric(scale(v)); d$zE <- z(d$E_sept); d$zla <- z(d$la); d$zivs <- z(d$ivs); d$zbmi <- z(d$bmi)
x <- d[d$sev1<2 & d$ms0==0 & !is.na(d$ms0),]; zz <- mk(x$t_ms, x$t_death, x$t_last); x$tt <- zz$tt; x$ev <- as.integer(zz$ev==1); x$id <- seq_len(nrow(x)); x$t_mac <- ifelse(!is.na(x$t_first) & x$t_first/365.25 < x$tt, x$t_first/365.25, NA)
base <- x[,c("id","tt","ev","age10","male","era_n","setting","av","ef_cat","zE","zla","zivs","zbmi","af","eg","esrd","sev1")]; tm <- tmerge(base, base, id=id, endpt=event(tt, ev)); tm <- tmerge(tm, x[!is.na(x$t_mac),c("id","t_mac")], id=id, mac=tdc(t_mac)); tm$mac[is.na(tm$mac)] <- 0
m <- coxph(Surv(tstart, tstop, endpt) ~ mac + age10 + male + era_n + setting + av + ef_cat + zE + zla + zivs + zbmi + af + eg + esrd + strata(sev1), data=tm); s_ <- summary(m)$conf.int
res2$ms_tvc <- data.frame(term=rownames(s_), hr=s_[,1], lo=s_[,3], hi=s_[,4]); res2$ms_tvc_n <- c(n=nrow(x), events=sum(x$ev), n_complete=m$n)
g <- d[d$sev1<2 & d$first_status %in% c("confirmed","refuted") & (is.na(d$t_ms) | d$t_ms > d$t_first),]
g$tt2 <- pmax((ifelse(!is.na(g$t_ms), g$t_ms, ifelse(!is.na(g$t_death), pmin(g$t_death, g$t_last+30), g$t_last)) - g$t_first)/365.25, 1/365.25); g$ev2 <- as.integer(!is.na(g$t_ms))
m2 <- coxph(Surv(tt2, ev2) ~ I(first_status=="confirmed") + age10 + male + strata(sev1), data=g); s2 <- summary(m2)$conf.int
res2$internal_check <- c(hr=s2[1,1], lo=s2[1,3], hi=s2[1,4], n=nrow(g), events=sum(g$ev2), conf_events=sum(g$ev2[g$first_status=="confirmed"]), conf_n=sum(g$first_status=="confirmed"), ref_events=sum(g$ev2[g$first_status=="refuted"]), ref_n=sum(g$first_status=="refuted"))
# ---------- INTERVENTION ----------
a <- d[!is.na(d$t_modany),]; ph <- a[!is.na(a$t_pheno),]; ph <- ph[!( !is.na(ph$t_interv) & ph$t_interv < ph$t_pheno),]; ph$iv <- !is.na(ph$t_interv) & ph$t_interv >= ph$t_pheno
res2$interv <- list(n_modplus=nrow(a), n_pheno=sum(!is.na(a$t_pheno)), pheno_mr=sum(a$pheno_type=="mr",na.rm=T), pheno_sten=sum(a$pheno_type=="sten",na.rm=T), n_analysed=nrow(ph), n_interv=sum(ph$iv), med_days=median((ph$t_interv-ph$t_pheno)[ph$iv]), isolated=sum(ph$interv_isolated[ph$iv]==1,na.rm=T), modality=table(ph$interv_modality[ph$iv]), died=sum(!is.na(ph$t_death)), med_days_death=median((ph$t_death-ph$t_pheno)[!is.na(ph$t_death)]))
ph$sten <- as.integer(ph$pheno_type=="sten"); ph$age10p <- (ph$age0 + ph$t_pheno/365.25 - 75)/10
# intervention by type of dysfunction, and the same cascade under the confirmed outcome definition (index moderate or severe, or a confirmed event)
res2$interv$by_type <- c(sten_n=sum(ph$sten==1), sten_iv=sum(ph$iv[ph$sten==1]), mr_n=sum(ph$sten==0), mr_iv=sum(ph$iv[ph$sten==0]))
confsel <- d$sev1>=2 | !is.na(d$t_conf); ac <- d[confsel,]; phc <- ac[!is.na(ac$t_pheno),]; phc <- phc[!( !is.na(phc$t_interv) & phc$t_interv < phc$t_pheno),]; ivc <- !is.na(phc$t_interv) & phc$t_interv >= phc$t_pheno
res2$interv_confirmed <- c(n_modplus=nrow(ac), n_index=sum(ac$sev1>=2), n_pheno=sum(!is.na(ac$t_pheno)), n_analysed=nrow(phc), n_interv=sum(ivc), died=sum(!is.na(phc$t_death)))
lg <- glm(iv ~ age10p + male + ef_cat + eg + esrd + sten + af, data=ph, family=binomial); s3 <- exp(cbind(coef(lg), confint.default(lg))); res2$interv_logit <- data.frame(term=rownames(s3), or=s3[,1], lo=s3[,2], hi=s3[,3]); res2$interv_logit_n <- c(n=nobs(lg), events=sum(ph$iv[!is.na(ph$lvef)&!is.na(ph$eg)]))
ph$tt3 <- pmax((ifelse(!is.na(ph$t_death), ph$t_death, pmax(ph$t_last, ph$t_pheno)) - ph$t_pheno)/365.25, 1/365.25); ph$dead <- as.integer(!is.na(ph$t_death)); ph$id <- seq_len(nrow(ph)); ph$t_iv <- ifelse(ph$iv & (ph$t_interv-ph$t_pheno)/365.25 < ph$tt3, (ph$t_interv-ph$t_pheno)/365.25, NA)
b <- ph[,c("id","tt3","dead","age10p","male","ef_cat","eg","esrd","sten","af")]; tm <- tmerge(b,b,id=id,death=event(tt3,dead)); tm <- tmerge(tm, ph[!is.na(ph$t_iv),c("id","t_iv")], id=id, ivt=tdc(t_iv)); tm$ivt[is.na(tm$ivt)] <- 0
m4 <- coxph(Surv(tstart,tstop,death) ~ ivt + age10p + male + ef_cat + eg + esrd + sten + af, data=tm); s4 <- summary(m4)$conf.int; res2$interv_surv <- data.frame(term=rownames(s4), hr=s4[,1], lo=s4[,3], hi=s4[,4])
res2$flow <- c(eligible=d$n_adult[1], excl_pros=d$n_pros[1], excl_rheum=d$n_rheum[1], analysed=nrow(d), blank=sum(d$sev1==0), mild=sum(d$sev1==1), modsev=sum(d$sev1>=2))
saveRDS(res2, "outputs/results/results_part2.rds"); cat("PART 2 saved\n"); print(res2$hmm$derived); print(res2$internal_check)
