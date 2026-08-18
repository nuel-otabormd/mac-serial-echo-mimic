# 06_secondary_mi.R: secondary models refitted across the 20 imputed datasets and pooled by Rubin's rules:
#   (a) calcific mitral stenosis with reaching moderate or severe MAC as a time-varying exposure (Cox, stratified by index grade),
#       with the death-window sensitivity (0, 30, 90, 365 days),
#   (b) receipt of intervention among patients with moderate or greater mitral dysfunction: cause-specific Cox model of the
#       hazard of intervention from the echocardiogram documenting dysfunction, death treated as a competing event,
#   (c) survival after documented dysfunction with intervention as a time-varying exposure (Cox; confounded by indication).
# Cohorts as in 05_reliability_intervention.R; only ef_cat / eg / af / esrd / echo covariates come from the imputed data.
suppressPackageStartupMessages({library(mice); library(survival)}); source("R/00_common.R")
d <- readRDS("data/frame.rds"); mi <- readRDS("data/mice.rds"); imp <- mi$imp; pid <- mi$pid
d$t_vs_end <- vital_horizon(d$t_last_disch, d$t_last_study); d$dead_vital <- as.integer(!is.na(d$t_death) & d$t_death <= d$t_vs_end); d$t_end_vital <- ifelse(d$dead_vital==1, d$t_death, d$t_vs_end)
mk_ef <- function(v) factor(ifelse(v < 40, "lt40", ifelse(v < 50, "40to49", "ge50")), levels=c("ge50","lt40","40to49"))
mk_eg <- function(v) factor(ifelse(v >= 60, "ge60", ifelse(v >= 30, "30to60", "lt30")), levels=c("ge60","30to60","lt30"))
MU <- c(E_sept=mean(d$E_sept,na.rm=T), la=mean(d$la,na.rm=T), ivs=mean(d$ivs,na.rm=T), bmi=mean(d$bmi,na.rm=T)); SD <- c(E_sept=sd(d$E_sept,na.rm=T), la=sd(d$la,na.rm=T), ivs=sd(d$ivs,na.rm=T), bmi=sd(d$bmi,na.rm=T))
# ---- (a) calcific stenosis: reaching moderate or severe MAC as a time-varying exposure, MI-pooled, with the death-window sensitivity ----
sel <- d$sev1<2 & d$ms0==0 & !is.na(d$ms0); x0 <- d[sel,]; x0$id <- seq_len(nrow(x0)); x0$age10 <- (x0$age0-70)/10; x0$era_n <- as.numeric(factor(d$era))[sel]-3
idx0 <- match(x0$pid, pid); stopifnot(!any(is.na(idx0)))
tvc_fit <- function(grace){ zz <- mk(x0$t_ms, x0$t_death, x0$t_last_study, grace); x0$tt <- zz$tt; x0$ev <- as.integer(zz$ev==1); x0$t_mac <- ifelse(!is.na(x0$t_first) & x0$t_first/365.25 < x0$tt, x0$t_first/365.25, NA)
  fits <- vector("list", imp$m)
  for (i in seq_len(imp$m)) { cd <- complete(imp, i)[idx0,]; x <- x0
    x$av <- as.numeric(cd$av); x$ef_cat <- mk_ef(cd$lvef); x$eg <- mk_eg(cd$egfr_base); x$af <- as.integer(as.logical(cd$af)); x$esrd <- as.integer(as.logical(cd$esrd))
    x$zE <- (cd$E_sept-MU["E_sept"])/SD["E_sept"]; x$zla <- (cd$la-MU["la"])/SD["la"]; x$zivs <- (cd$ivs-MU["ivs"])/SD["ivs"]; x$zbmi <- (cd$bmi-MU["bmi"])/SD["bmi"]
    base <- x[,c("id","tt","ev","age10","male","era_n","setting","av","ef_cat","zE","zla","zivs","zbmi","af","eg","esrd","sev1")]; tm <- tmerge(base, base, id=id, endpt=event(tt, ev)); tm <- tmerge(tm, x[!is.na(x$t_mac),c("id","t_mac")], id=id, mac=tdc(t_mac)); tm$mac[is.na(tm$mac)] <- 0
    fits[[i]] <- coxph(Surv(tstart, tstop, endpt) ~ mac + age10 + male + era_n + setting + av + ef_cat + zE + zla + zivs + zbmi + af + eg + esrd + strata(sev1), data=tm) }
  tp <- summary(pool(as.mira(fits)), conf.int=TRUE, exponentiate=TRUE)
  list(tab=data.frame(term=tp$term, hr=tp$estimate, lo=tp$`2.5 %`, hi=tp$`97.5 %`, p=tp$p.value), n=c(n=nrow(x0), events=sum(x0$ev))) }
main_tvc <- tvc_fit(GRACE_DAYS); ms_tvc_mi <- main_tvc$tab; ms_tvc_mi_n <- main_tvc$n
cat(sprintf("stenosis, reaching moderate or severe MAC (time-varying), MI-pooled: HR %.2f (%.2f-%.2f); n = %d, events = %d\n", ms_tvc_mi$hr[ms_tvc_mi$term=="mac"], ms_tvc_mi$lo[ms_tvc_mi$term=="mac"], ms_tvc_mi$hi[ms_tvc_mi$term=="mac"], ms_tvc_mi_n["n"], ms_tvc_mi_n["events"]))
ms_tvc_grace <- do.call(rbind, lapply(c(0,30,90,365), function(g){ r <- tvc_fit(g); z <- r$tab[r$tab$term=="mac",]; data.frame(grace_days=g, hr=z$hr, lo=z$lo, hi=z$hi, n=r$n["n"], events=r$n["events"]) }))
print(ms_tvc_grace, row.names=FALSE)
# ---- (b) intervention: cause-specific hazard of intervention from the dysfunction echocardiogram, death competing, MI-pooled ----
a <- d[!is.na(d$t_modany),]; ph <- a[!is.na(a$t_pheno),]; ph <- ph[!( !is.na(ph$t_interv) & ph$t_interv < ph$t_pheno),]
ph$iv <- !is.na(ph$t_interv) & ph$t_interv >= ph$t_pheno; ph$sten <- as.integer(ph$pheno_type=="sten"); ph$age10p <- (ph$age0 + ph$t_pheno/365.25 - 75)/10
ph$fu_end <- ifelse(ph$dead_vital==1, ph$t_death, ph$t_end_vital); ph$tt3 <- pmax(ph$fu_end - ph$t_pheno, 1)/365.25
ph$t_iv_y <- ifelse(ph$iv, (ph$t_interv - ph$t_pheno)/365.25, NA); ph$ev_iv <- as.integer(ph$iv & ph$t_iv_y <= ph$tt3); ph$tt_iv <- ifelse(ph$ev_iv==1, pmax(ph$t_iv_y, 1/365.25), ph$tt3)
idx <- match(ph$pid, pid); stopifnot(!any(is.na(idx)))
ph$ef_cat <- factor(as.character(ph$ef_cat), levels=c("ge50","lt40","40to49")); ph$eg <- factor(as.character(ph$eg), levels=c("ge60","30to60","lt30"))
cc <- coxph(Surv(tt_iv, ev_iv) ~ age10p + male + ef_cat + eg + esrd + sten + af, data=ph)
fits <- vector("list", imp$m)
for (i in seq_len(imp$m)) { cd <- complete(imp, i)[idx,]; x <- ph
  x$ef_cat <- mk_ef(cd$lvef); x$eg <- mk_eg(cd$egfr_base); x$af <- as.integer(as.logical(cd$af)); x$esrd <- as.integer(as.logical(cd$esrd))
  fits[[i]] <- coxph(Surv(tt_iv, ev_iv) ~ age10p + male + ef_cat + eg + esrd + sten + af, data=x) }
pooled <- summary(pool(as.mira(fits)), conf.int=TRUE, exponentiate=TRUE)
out <- data.frame(term=pooled$term, hr_mi=pooled$estimate, lo_mi=pooled$`2.5 %`, hi_mi=pooled$`97.5 %`, p_mi=pooled$p.value)
ccs <- summary(cc)$conf.int; m <- match(out$term, rownames(ccs)); out$hr_cc <- ccs[m,1]; out$lo_cc <- ccs[m,3]; out$hi_cc <- ccs[m,4]
print(out, row.names=FALSE, digits=3)
cat("MI: n =", nrow(ph), " interventions =", sum(ph$ev_iv), " | complete-case: n =", cc$n, " interventions =", cc$nevent, "\n")
# ---- (c) survival after documented dysfunction, intervention as a time-varying exposure, MI-pooled ----
ph$id <- seq_len(nrow(ph)); ph$t_iv <- ifelse(ph$iv & ph$t_iv_y < ph$tt3, ph$t_iv_y, NA)
sfits <- vector("list", imp$m)
for (i in seq_len(imp$m)) { cd <- complete(imp, i)[idx,]; x <- ph
  x$ef_cat <- mk_ef(cd$lvef); x$eg <- mk_eg(cd$egfr_base); x$af <- as.integer(as.logical(cd$af)); x$esrd <- as.integer(as.logical(cd$esrd))
  b <- x[,c("id","tt3","dead_vital","age10p","male","ef_cat","eg","esrd","sten","af")]; tm <- tmerge(b,b,id=id,death=event(tt3,dead_vital)); tm <- tmerge(tm, x[!is.na(x$t_iv),c("id","t_iv")], id=id, ivt=tdc(t_iv)); tm$ivt[is.na(tm$ivt)] <- 0
  sfits[[i]] <- coxph(Surv(tstart,tstop,death) ~ ivt + age10p + male + ef_cat + eg + esrd + sten + af, data=tm) }
sp <- summary(pool(as.mira(sfits)), conf.int=TRUE, exponentiate=TRUE)
cat("Survival, intervention (time-varying) HR pooled:", round(sp$estimate[sp$term=="ivt"],2), round(sp$`2.5 %`[sp$term=="ivt"],2), round(sp$`97.5 %`[sp$term=="ivt"],2), "\n")
saveRDS(list(pooled=out, n=nrow(ph), events=sum(ph$ev_iv), n_cc=cc$n, events_cc=cc$nevent, ms_tvc=ms_tvc_mi, ms_tvc_n=ms_tvc_mi_n, ms_tvc_grace=ms_tvc_grace,
             surv_hr=c(hr=sp$estimate[sp$term=="ivt"], lo=sp$`2.5 %`[sp$term=="ivt"], hi=sp$`97.5 %`[sp$term=="ivt"])), "outputs/results/secondary_mi.rds")
cat("secondary_mi.rds saved\n")
