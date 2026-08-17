# 06_secondary_mi.R: secondary models refitted across the 20 imputed datasets and pooled by Rubin's rules:
#   (a) calcific mitral stenosis with reaching moderate or greater MAC as a time-varying exposure (Cox, stratified by index grade),
#   (b) receipt of intervention among patients with moderate or greater mitral dysfunction (logistic),
#   (c) survival after documented dysfunction with intervention as a time-varying exposure (Cox).
# Intervention logistic model (Table 4C) refitted across the 20 imputed datasets (data/mice.rds) and pooled by Rubin's rules.
# Cohort and model as in 05_reliability_intervention.R; only ef_cat / eg / af / esrd come from the imputed data.
suppressPackageStartupMessages(library(mice))
d <- readRDS("data/frame.rds"); mi <- readRDS("data/mice.rds"); imp <- mi$imp; pid <- mi$pid
a <- d[!is.na(d$t_modany),]; ph <- a[!is.na(a$t_pheno),]; ph <- ph[!( !is.na(ph$t_interv) & ph$t_interv < ph$t_pheno),]
ph$iv <- !is.na(ph$t_interv) & ph$t_interv >= ph$t_pheno; ph$sten <- as.integer(ph$pheno_type=="sten"); ph$age10p <- (ph$age0 + ph$t_pheno/365.25 - 75)/10
cat("classes: ef_cat", class(d$ef_cat), "levels", levels(d$ef_cat), "| eg", class(d$eg), levels(d$eg), "| esrd", class(d$esrd), "| af", class(d$af), "| imp esrd", class(imp$data$esrd), "| imp af", class(imp$data$af), "\n")
cat("frame esrd table:", table(d$esrd, useNA="a"), "| imp esrd table:", table(imp$data$esrd, useNA="a"), "\n")
idx <- match(ph$pid, pid); stopifnot(!any(is.na(idx)))
efl <- levels(factor(d$ef_cat)); egl <- levels(factor(d$eg))
mk_ef <- function(v) factor(ifelse(v < 40, "lt40", ifelse(v < 50, "40to49", "ge50")), levels=c("ge50","lt40","40to49"))
mk_eg <- function(v) factor(ifelse(v >= 60, "ge60", ifelse(v >= 30, "30to60", "lt30")), levels=c("ge60","30to60","lt30"))
# sanity: derived categories reproduce the frame's own categories where observed
chk <- !is.na(ph$lvef); cat("ef_cat reproduces frame:", mean(as.character(mk_ef(ph$lvef[chk]))==as.character(ph$ef_cat[chk])), "\n")
chk2 <- !is.na(ph$egfr_base); cat("eg reproduces frame:", mean(as.character(mk_eg(ph$egfr_base[chk2]))==as.character(ph$eg[chk2])), "\n")
ph$ef_cat <- factor(as.character(ph$ef_cat), levels=c("ge50","lt40","40to49")); ph$eg <- factor(as.character(ph$eg), levels=c("ge60","30to60","lt30"))
# ---- (a) calcific stenosis: reaching moderate or greater MAC as a time-varying exposure, MI-pooled ----
suppressPackageStartupMessages(library(survival))
mk <- function(tp, td, tl){ ev <- ifelse(!is.na(tp),1, ifelse(!is.na(td) & (is.na(tl) | td<=tl+30),2,0)); tt <- ifelse(!is.na(tp),tp, ifelse(ev==2, td, tl)); list(ev=ev, tt=pmax(tt,1)/365.25) }
MU <- c(E_sept=mean(d$E_sept,na.rm=T), la=mean(d$la,na.rm=T), ivs=mean(d$ivs,na.rm=T), bmi=mean(d$bmi,na.rm=T)); SD <- c(E_sept=sd(d$E_sept,na.rm=T), la=sd(d$la,na.rm=T), ivs=sd(d$ivs,na.rm=T), bmi=sd(d$bmi,na.rm=T))
sel <- d$sev1<2 & d$ms0==0 & !is.na(d$ms0); x0 <- d[sel,]; zz <- mk(x0$t_ms, x0$t_death, x0$t_last); x0$tt <- zz$tt; x0$ev <- as.integer(zz$ev==1); x0$id <- seq_len(nrow(x0))
x0$t_mac <- ifelse(!is.na(x0$t_first) & x0$t_first/365.25 < x0$tt, x0$t_first/365.25, NA); x0$age10 <- (x0$age0-70)/10; x0$era_n <- as.numeric(factor(d$era))[sel]-3
idx0 <- match(x0$pid, pid); tvc <- vector("list", imp$m)
for (i in seq_len(imp$m)) {
  cd <- complete(imp, i)[idx0,]; x <- x0
  x$av <- as.numeric(cd$av); x$ef_cat <- mk_ef(cd$lvef); x$eg <- mk_eg(cd$egfr_base); x$af <- as.integer(as.logical(cd$af)); x$esrd <- as.integer(as.logical(cd$esrd))
  x$zE <- (cd$E_sept-MU["E_sept"])/SD["E_sept"]; x$zla <- (cd$la-MU["la"])/SD["la"]; x$zivs <- (cd$ivs-MU["ivs"])/SD["ivs"]; x$zbmi <- (cd$bmi-MU["bmi"])/SD["bmi"]
  base <- x[,c("id","tt","ev","age10","male","era_n","setting","av","ef_cat","zE","zla","zivs","zbmi","af","eg","esrd","sev1")]; tm <- tmerge(base, base, id=id, endpt=event(tt, ev)); tm <- tmerge(tm, x[!is.na(x$t_mac),c("id","t_mac")], id=id, mac=tdc(t_mac)); tm$mac[is.na(tm$mac)] <- 0
  tvc[[i]] <- coxph(Surv(tstart, tstop, endpt) ~ mac + age10 + male + era_n + setting + av + ef_cat + zE + zla + zivs + zbmi + af + eg + esrd + strata(sev1), data=tm) }
tp <- summary(pool(as.mira(tvc)), conf.int=TRUE, exponentiate=TRUE)
ms_tvc_mi <- data.frame(term=tp$term, hr=tp$estimate, lo=tp$`2.5 %`, hi=tp$`97.5 %`, p=tp$p.value); ms_tvc_mi_n <- c(n=nrow(x0), events=sum(x0$ev))
cat(sprintf("stenosis, reaching moderate or greater MAC (time-varying), MI-pooled: HR %.2f (%.2f-%.2f); n = %d, events = %d\n", ms_tvc_mi$hr[ms_tvc_mi$term=="mac"], ms_tvc_mi$lo[ms_tvc_mi$term=="mac"], ms_tvc_mi$hi[ms_tvc_mi$term=="mac"], ms_tvc_mi_n["n"], ms_tvc_mi_n["events"]))
# ---- (b) intervention logistic ----
cc <- glm(iv ~ age10p + male + ef_cat + eg + esrd + sten + af, data=ph, family=binomial)
fits <- vector("list", imp$m)
for (i in seq_len(imp$m)) {
  cd <- complete(imp, i)[idx,]; x <- ph
  x$ef_cat <- mk_ef(cd$lvef); x$eg <- mk_eg(cd$egfr_base); x$af <- as.integer(as.logical(cd$af)); x$esrd <- as.integer(as.logical(cd$esrd))
  fits[[i]] <- glm(iv ~ age10p + male + ef_cat + eg + esrd + sten + af, data=x, family=binomial)
}
pooled <- summary(pool(as.mira(fits)), conf.int=TRUE, exponentiate=TRUE)
out <- data.frame(term=pooled$term, or_mi=round(pooled$estimate,2), lo_mi=round(pooled$`2.5 %`,2), hi_mi=round(pooled$`97.5 %`,2), p_mi=signif(pooled$p.value,2))
ccs <- exp(cbind(coef(cc), confint.default(cc))); m <- match(out$term, rownames(ccs)); out$or_cc <- round(ccs[m,1],2); out$lo_cc <- round(ccs[m,2],2); out$hi_cc <- round(ccs[m,3],2)
print(out, row.names=FALSE)
cat("MI: n =", nrow(ph), " events =", sum(ph$iv), " | complete-case: n =", nobs(cc), " events =", sum(ph$iv[!is.na(ph$lvef)&!is.na(ph$egfr_base)]), "\n")
saveRDS(list(pooled=out, n=nrow(ph), events=sum(ph$iv), n_cc=nobs(cc), events_cc=sum(ph$iv[!is.na(ph$lvef)&!is.na(ph$egfr_base)]), ms_tvc=ms_tvc_mi, ms_tvc_n=ms_tvc_mi_n), "outputs/results/secondary_mi.rds")
# ---- companion survival model: intervention as time-varying exposure, same imputations ----
suppressPackageStartupMessages(library(survival))
ph$tt3 <- pmax((ifelse(!is.na(ph$t_death), ph$t_death, pmax(ph$t_last, ph$t_pheno)) - ph$t_pheno)/365.25, 1/365.25); ph$dead <- as.integer(!is.na(ph$t_death)); ph$id <- seq_len(nrow(ph))
ph$t_iv <- ifelse(ph$iv & (ph$t_interv-ph$t_pheno)/365.25 < ph$tt3, (ph$t_interv-ph$t_pheno)/365.25, NA)
sfits <- vector("list", imp$m)
for (i in seq_len(imp$m)) {
  cd <- complete(imp, i)[idx,]; x <- ph
  x$ef_cat <- mk_ef(cd$lvef); x$eg <- mk_eg(cd$egfr_base); x$af <- as.integer(as.logical(cd$af)); x$esrd <- as.integer(as.logical(cd$esrd))
  b <- x[,c("id","tt3","dead","age10p","male","ef_cat","eg","esrd","sten","af")]; tm <- tmerge(b,b,id=id,death=event(tt3,dead)); tm <- tmerge(tm, x[!is.na(x$t_iv),c("id","t_iv")], id=id, ivt=tdc(t_iv)); tm$ivt[is.na(tm$ivt)] <- 0
  sfits[[i]] <- coxph(Surv(tstart,tstop,death) ~ ivt + age10p + male + ef_cat + eg + esrd + sten + af, data=tm)
}
sp <- summary(pool(as.mira(sfits)), conf.int=TRUE, exponentiate=TRUE)
cat("Survival, intervention (time-varying) HR pooled:", round(sp$estimate[sp$term=="ivt"],2), round(sp$`2.5 %`[sp$term=="ivt"],2), round(sp$`97.5 %`[sp$term=="ivt"],2), "\n")
b0 <- ph[,c("id","tt3","dead","age10p","male","ef_cat","eg","esrd","sten","af")]; tm0 <- tmerge(b0,b0,id=id,death=event(tt3,dead)); tm0 <- tmerge(tm0, ph[!is.na(ph$t_iv),c("id","t_iv")], id=id, ivt=tdc(t_iv)); tm0$ivt[is.na(tm0$ivt)] <- 0
m0 <- coxph(Surv(tstart,tstop,death) ~ ivt + age10p + male + ef_cat + eg + esrd + sten + af, data=tm0); cat("complete-case HR:", round(summary(m0)$conf.int["ivt",c(1,3,4)],2), " n used:", m0$n, "\n")
r <- readRDS("outputs/results/secondary_mi.rds"); r$surv_hr <- c(hr=sp$estimate[sp$term=="ivt"], lo=sp$`2.5 %`[sp$term=="ivt"], hi=sp$`97.5 %`[sp$term=="ivt"]); saveRDS(r, "outputs/results/secondary_mi.rds")
