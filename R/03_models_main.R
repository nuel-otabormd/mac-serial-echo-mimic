# 03_models_main.R: co-primary models (piecewise-exponential on person-intervals split at the time-band cut points, patient-
# clustered variances, pooled across imputations), complementary log-log companion on the unsplit intervals (exact
# interval-censored likelihood, hazard constant within an interval), confirmed-event timing sensitivity,
# Cox companion, imputation bridge, stratum interaction, stricter onset cohort, landmark at the second episode, complete-case
# models, sensitivity excluding patients whose only rheumatic report came after index, visit-process (inverse-intensity) and
# censoring weights with weight diagnostics, total versus direct associations, Fine and Gray pooled across imputations, and
# the laboratory-definition and care-setting sensitivities.
suppressPackageStartupMessages({library(survival); library(mice); library(sandwich)}); source("R/00_common.R")
set.seed(20260816)
M_USE <- as.integer(Sys.getenv("MAC_M", "20"))          # number of imputed datasets to use (20 for the reported results; smaller for smoke tests)
dA <- readRDS("data/frame.rds"); d <- dA
p <- read.csv("data/panel.csv"); p <- p[order(p$pid,p$rn),]
p$t_prev <- ave(p$t, p$pid, FUN=function(v) c(NA, head(v,-1))); p$sev_prev <- ave(p$sev, p$pid, FUN=function(v) c(NA, head(v,-1))); p$inpt_prev <- ave(p$inpt, p$pid, FUN=function(v) c(NA, head(v,-1)))
M <- readRDS("data/mice.rds"); imp <- M$imp; stopifnot(all(M$pid==dA$pid))
res <- list(settings=list(grace_days=GRACE_DAYS, band_cuts=BAND_CUTS, weight_cap=c(0.1,10), m=M_USE, variance="patient-clustered sandwich (HC0), pooled by Rubin's rules"))
msg <- function(...) cat(sprintf(...), "\n")
# ---------- standardisation constants from observed data (primary cohort) ----------
SD <- c(E_sept=sd(d$E_sept,na.rm=T), la=sd(d$la,na.rm=T), ivs=sd(d$ivs,na.rm=T), bmi=sd(d$bmi,na.rm=T), phos=sd(d$phos_median,na.rm=T), ca=sd(d$ca_median,na.rm=T), alp=sd(log(d$alp_median),na.rm=T))
MU <- c(E_sept=mean(d$E_sept,na.rm=T), la=mean(d$la,na.rm=T), ivs=mean(d$ivs,na.rm=T), bmi=mean(d$bmi,na.rm=T), phos=mean(d$phos_median,na.rm=T), ca=mean(d$ca_median,na.rm=T), alp=mean(log(d$alp_median),na.rm=T))
derive_cov <- function(x){   # x: one imputed dataset (or the observed data) with the raw index variables; returns the analysis covariates
  x$age10 <- (x$age0-70)/10; x$setting <- factor(as.character(x$setting), levels=c("outpt","ED","ward","ICU"))
  x$av_cat <- factor(pmin(pmax(round(as.numeric(x$av)),0),3), levels=0:3, labels=c("nl","mild","mod","sev")); x$av_lin <- as.numeric(x$av_cat)-1
  x$ef_cat <- relevel(cut(x$lvef, c(-Inf,39.99,49.99,Inf), labels=c("lt40","40to49","ge50")), "ge50")
  x$eg <- relevel(cut(x$egfr_base, c(-Inf,30,60,Inf), labels=c("lt30","30to60","ge60")), "ge60")
  x$zE <- (x$E_sept-MU["E_sept"])/SD["E_sept"]; x$zla <- (x$la-MU["la"])/SD["la"]; x$zivs <- (x$ivs-MU["ivs"])/SD["ivs"]; x$zbmi <- (x$bmi-MU["bmi"])/SD["bmi"]
  x$zphos <- (x$phos_median-MU["phos"])/SD["phos"]; x$zca <- (x$ca_median-MU["ca"])/SD["ca"]; x$zalp <- (log(x$alp_median)-MU["alp"])/SD["alp"]
  x }
FIXED <- c("pid","yr1","t_first","t_conf","t_last","t_last_study","t_death","first_status","n_ep","idx_inpt","rheum_post","dm_prior","htn_prior","cad_prior","af_prior","esrd_prior")   # *_prior: diagnosis codes only from admissions completed by time zero (sensitivity)
covs_m <- function(m){ x <- complete(imp, m); x$pid <- dA$pid; x$sev1n <- dA$sev1; x$era_n <- x$era_n
  for (v in FIXED[-1]) x[[v]] <- dA[[v]]
  x <- derive_cov(x)
  x[, c(FIXED,"age10","male","era_n","setting","sev1n","av_cat","av_lin","ef_cat","eg","esrd","zE","zla","zivs","zbmi","af","zphos","zca","zalp","dm","htn","cad")] }
covs <- lapply(seq_len(M_USE), covs_m); covsX <- lapply(covs, function(x) x[x$rheum_post==0,])   # covs = primary cohort; covsX excludes the patients whose only rheumatic report came after index
# ---------- person-interval builder (intervals between consecutive episodes, split at the band cut points) ----------
build_iv <- function(cov, stratum, def, from_rn=1, t_event_override=NULL){
  ids <- cov$pid[cov$sev1n==stratum]; pp <- p[p$pid %in% ids & p$rn>from_rn, c("pid","rn","t","t_prev","sev_prev","inpt_prev")]
  if (from_rn>1) { t0 <- p$t[p$rn==from_rn]; names(t0) <- p$pid[p$rn==from_rn]; pp$t_prev <- pmax(pp$t_prev, t0[pp$pid]) }
  te <- if (!is.null(t_event_override)) t_event_override else if (def=="first") cov$t_first else cov$t_conf; names(te) <- cov$pid
  pp$te <- te[pp$pid]; pp$ev <- as.integer(!is.na(pp$te) & pp$t==pp$te)
  pp$after <- ave(pp$ev, pp$pid, FUN=function(v) c(0, head(cumsum(v),-1))); pp <- pp[pp$after==0,]   # drop intervals after the event
  if (from_rn>1) { pp <- pp[pp$t > pp$t_prev,]; pp <- pp[is.na(pp$te) | pp$te > (t0[pp$pid]),] }      # landmark: follow-up starts at the landmark episode; events at or before it excluded (they define it)
  pp <- split_intervals(pp[,c("pid","t_prev","t","ev","sev_prev","inpt_prev")])
  merge(pp, cov, by="pid") }
build_iv_unsplit <- function(cov, stratum, def){
  # one row per observation interval, NOT split at the band cut points: with a complementary log-log link this is the exact
  # likelihood for an event known only to lie within the interval, the hazard being taken as constant over each interval and
  # indexed by the band in which the interval starts (splitting and keeping the event in the last segment would instead
  # assert that the event fell after the last cut point crossed, which is not known)
  ids <- cov$pid[cov$sev1n==stratum]; pp <- p[p$pid %in% ids & p$rn>1, c("pid","rn","t","t_prev","sev_prev","inpt_prev")]
  te <- if (def=="first") cov$t_first else cov$t_conf; names(te) <- cov$pid
  pp$te <- te[pp$pid]; pp$ev <- as.integer(!is.na(pp$te) & pp$t==pp$te)
  pp$after <- ave(pp$ev, pp$pid, FUN=function(v) c(0, head(cumsum(v),-1))); pp <- pp[pp$after==0,]
  pp$dt <- (pmax(pp$t, pp$t_prev+1)-pp$t_prev)/365.25; pp$band <- band_of(pp$t_prev/365.25)
  merge(pp[,c("pid","t_prev","t","ev","dt","band","sev_prev","inpt_prev")], cov, by="pid") }
fit_glm <- function(iv, form, wts=NULL, family=poisson()){
  f <- suppressWarnings(glm(as.formula(paste("ev ~ band + offset(log(dt)) +", form)), data=iv, family=family, weights=wts))
  stopifnot(nobs(f)==nrow(iv))                                            # no silent row loss: every fitted model's n is the n of its data
  V <- sandwich::vcovCL(f, cluster=iv$pid)
  list(coef=coef(f), var=diag(V), n_pt=length(unique(iv$pid)), n_iv=nrow(iv), events=sum(iv$ev)) }
pool_hr <- function(fits, terms){ terms <- terms[terms %in% names(fits[[1]]$coef)]
  Q <- sapply(fits, function(f) f$coef[terms]); U <- sapply(fits, function(f) f$var[terms]); if (is.null(dim(Q))) { Q <- matrix(Q,nrow=1,dimnames=list(terms,NULL)); U <- matrix(U,nrow=1) }
  qb <- rowMeans(Q); ub <- rowMeans(U); b <- if (ncol(Q)>1) apply(Q,1,var) else 0; T <- ub + (1+1/ncol(Q))*b; z <- qb/sqrt(T)
  data.frame(term=terms, hr=exp(qb), lo=exp(qb-1.96*sqrt(T)), hi=exp(qb+1.96*sqrt(T)), p=2*pnorm(-abs(z)), row.names=NULL) }
D1 <- "age10 + male + era_n + setting"; D2 <- paste(D1, "+ av_cat + ef_cat + zE + zla + zivs + zbmi + af"); D3 <- paste(D2, "+ eg + esrd"); D4 <- paste(D3, "+ zphos + zca + zalp"); D5 <- paste(D3, "+ dm + htn + cad")
KEY <- c("age10","male","av_catmild","av_catmod","av_catsev","ef_catlt40","ef_cat40to49","zE","zla","zivs","zbmi","af","eglt30","eg30to60","esrd","zphos","zca","zalp","dm","htn","cad")
lab_st <- function(s) ifelse(s==0,"onset","progression")
# ---------- 1. MAIN MODELS: piecewise-exponential (primary), MI-pooled, both strata x definitions x domains; Cox and cloglog companions ----------
main <- list(); clog <- list()
for (s in 0:1) for (def in c("first","confirmed")) {
  ivs <- lapply(covs, function(cv) build_iv(cv, s, def))
  for (dom in c("D1","D2","D3","D4","D5")) { form <- get(dom)
    ivd <- if (dom=="D4") lapply(ivs, function(x) x[!is.na(x$zphos) & !is.na(x$zca) & !is.na(x$zalp),]) else if (dom=="D5") lapply(ivs, function(x) x[x$yr1==1,]) else ivs   # D4: complete on all three markers, so the stated n is the fitted n
    fits <- lapply(ivd, function(x) fit_glm(x, form))
    out <- pool_hr(fits, KEY); out$stratum <- lab_st(s); out$def <- def; out$domain <- dom; out$n_pt <- fits[[1]]$n_pt; out$n_iv <- fits[[1]]$n_iv; out$events <- fits[[1]]$events
    main[[length(main)+1]] <- out
    msg("PE %-11s %-9s %s: patients %d intervals %d events %d", lab_st(s), def, dom, out$n_pt[1], out$n_iv[1], out$events[1]) }
  # complementary log-log companion on the UNSPLIT observation intervals: exact likelihood for an event known only to lie
  # within the interval (hazard constant over the interval, baseline indexed by the band at the interval start)
  fits <- lapply(covs, function(cv) fit_glm(build_iv_unsplit(cv, s, def), D3, family=binomial(link="cloglog")))
  out <- pool_hr(fits, KEY); out$stratum <- lab_st(s); out$def <- def; out$n_pt <- fits[[1]]$n_pt; out$events <- fits[[1]]$events; clog[[length(clog)+1]] <- out
  # Cox companion (event at the qualifying date; death within the grace window a competing event, otherwise censoring at the last study), MI-pooled, D3
  cx <- lapply(covs, function(cv){ x <- cv[cv$sev1n==s,]; tp <- if (def=="first") x$t_first else x$t_conf; z <- mk(tp, x$t_death, x$t_last_study); x$tt <- z$tt; x$ev <- z$ev
    f <- coxph(as.formula(paste("Surv(tt, ev==1) ~", D3)), data=x, cluster=pid); list(coef=coef(f), var=diag(vcov(f)), n_pt=f$n, events=f$nevent) })
  out <- pool_hr(cx, KEY); out$stratum <- lab_st(s); out$def <- def; out$domain <- "D3_cox_companion"; out$n_pt <- cx[[1]]$n_pt; out$n_iv <- NA; out$events <- cx[[1]]$events; main[[length(main)+1]] <- out
}
res$main <- do.call(rbind, main); res$cloglog <- do.call(rbind, clog)
# ---------- 2. MI BRIDGE (confirmation imputed for never-re-examined reads), joint with covariate imputation, PE D3 ----------
bridge <- list()
for (s in 0:1) { fits <- list()
  for (m in seq_len(M_USE)) { cv <- covs[[m]]; x <- cv[cv$sev1n==s,]; x$t_end <- ifelse(!is.na(x$t_death), x$t_death, x$t_last_study); x$t_after <- (x$t_end-x$t_first)/365.25; x$died <- !is.na(x$t_death)
    re <- x[x$first_status %in% c("confirmed","refuted"),]; un <- x[x$first_status=="unconfirmable",]
    im <- glm(I(first_status=="confirmed") ~ age10 + male + setting + era_n + ef_cat + esrd + av_lin + t_after + died, data=re, family=binomial)
    pun <- predict(im, newdata=un, type="response"); hit <- un$pid[rbinom(nrow(un),1,pun)==1]
    tconf <- x$t_conf; tconf[x$pid %in% hit] <- x$t_first[x$pid %in% hit]; names(tconf) <- x$pid
    full <- rep(NA_real_, nrow(cv)); names(full) <- cv$pid; full[names(tconf)] <- tconf
    fits[[m]] <- fit_glm(build_iv(cv, s, "confirmed", t_event_override=full), D3) }
  out <- pool_hr(fits, KEY); out$stratum <- lab_st(s); out$events_imp1 <- fits[[1]]$events; bridge[[s+1]] <- out; msg("BRIDGE %s events(imp1) %d", lab_st(s), fits[[1]]$events) }
res$bridge <- do.call(rbind, bridge)
# ---------- 3. INTERACTION (first primary; confirmed check), PE D3 with stratum interaction, MI-pooled ----------
inter <- list()
for (def in c("first","confirmed")) { fits <- lapply(covs, function(cv){ iv <- rbind(build_iv(cv,0,def), build_iv(cv,1,def)); iv$prog <- as.integer(iv$sev1n==1)
    f <- suppressWarnings(glm(ev ~ band*prog + offset(log(dt)) + (age10 + male + av_lin + ef_cat + zE + zla + zivs + eg + esrd + af + zbmi)*prog + era_n + setting, data=iv, family=poisson))
    list(coef=coef(f), var=diag(sandwich::vcovCL(f, cluster=iv$pid))) })
  it <- grep(":prog$|^prog:", names(fits[[1]]$coef), value=TRUE); it <- it[!grepl("band", it)]
  out <- pool_hr(fits, it); out$def <- def; inter[[def]] <- out }
res$interaction <- do.call(rbind, inter)
# ---------- 4. CONFIRMED MAC-FREE BASELINE (index blank & 2nd blank; from 2nd episode), PE D3, MI-pooled ----------
sec <- p[p$rn==2, c("pid","sev")]; ok2 <- sec$pid[sec$sev==0]
cb <- list()
for (def in c("first","confirmed")) { fits <- lapply(covs, function(cv) fit_glm(build_iv(cv[cv$pid %in% ok2,], 0, def, from_rn=2), D3))
  out <- pool_hr(fits, KEY); out$def <- def; out$n_pt <- fits[[1]]$n_pt; out$events <- fits[[1]]$events; cb[[def]] <- out; msg("CONFIRMED-MAC-FREE baseline %s: patients %d events %d", def, out$n_pt[1], out$events[1]) }
res$confirmed_baseline <- do.call(rbind, cb)
# ---------- 4b. LANDMARK AT THE SECOND EPISODE (both cohorts, no grade requirement at the second episode), PE D3, MI-pooled ----------
lm_ <- list()
for (s in 0:1) for (def in c("first","confirmed")) { fits <- lapply(covs, function(cv) fit_glm(build_iv(cv, s, def, from_rn=2), D3))
  out <- pool_hr(fits, KEY); out$stratum <- lab_st(s); out$def <- def; out$n_pt <- fits[[1]]$n_pt; out$events <- fits[[1]]$events; lm_[[length(lm_)+1]] <- out; msg("LANDMARK 2nd episode %s %s: patients %d events %d", lab_st(s), def, out$n_pt[1], out$events[1]) }
res$landmark <- do.call(rbind, lm_)
# ---------- 4c. COMPLETE-CASE D3 (observed covariates only, no imputation) ----------
cc0 <- derive_cov(transform(d, av=as.numeric(av))); cc0$sev1n <- d$sev1; cc0$era_n <- d$era_n
ccv <- c("age10","male","era_n","setting","av_cat","ef_cat","zE","zla","zivs","zbmi","af","eg","esrd"); cc <- cc0[complete.cases(cc0[, ccv]), c(FIXED, "sev1n", ccv)]
ccl <- list()
for (s in 0:1) for (def in c("first","confirmed")) { f1 <- fit_glm(build_iv(cc, s, def), D3); out <- pool_hr(list(f1), KEY); out$stratum <- lab_st(s); out$def <- def; out$n_pt <- f1$n_pt; out$events <- f1$events; ccl[[length(ccl)+1]] <- out
  msg("COMPLETE-CASE %s %s: patients %d events %d", lab_st(s), def, f1$n_pt, f1$events) }
res$complete_case <- do.call(rbind, ccl)
# ---------- 4d. SENSITIVITY: patients whose only rheumatic report came after index excluded, PE D3, MI-pooled ----------
rh <- list()
for (s in 0:1) for (def in c("first","confirmed")) { fits <- lapply(covsX, function(cv) fit_glm(build_iv(cv, s, def), D3))
  out <- pool_hr(fits, KEY); out$stratum <- lab_st(s); out$def <- def; out$n_pt <- fits[[1]]$n_pt; out$events <- fits[[1]]$events; rh[[length(rh)+1]] <- out }
res$rheum_excluded <- do.call(rbind, rh)
# ---------- 4e. SENSITIVITY: confirmed event dated at the CONFIRMING echocardiogram (the second of the pair) rather than the first, PE D3, MI-pooled ----------
ct <- list()
for (s in 0:1) { fits <- lapply(covs, function(cv) fit_glm(build_iv(cv, s, "confirmed", t_event_override=dA$t_conf_next[match(cv$pid, dA$pid)]), D3))
  out <- pool_hr(fits, KEY); out$stratum <- lab_st(s); out$def <- "confirmed, dated at the confirming echocardiogram"; out$n_pt <- fits[[1]]$n_pt; out$events <- fits[[1]]$events; ct[[length(ct)+1]] <- out
  msg("CONFIRMED TIMING (event at confirming echo) %s: patients %d events %d", lab_st(s), out$n_pt[1], out$events[1]) }
res$confirmed_timing <- do.call(rbind, ct)
# ---------- 4f. SENSITIVITY: diagnosis codes only from admissions COMPLETED by time zero (MIMIC-IV codes are assigned per admission at
# discharge and carry no date within the stay, so the main definition, admissions begun by time zero, can include conditions documented
# later in an index admission). D3 with prior-only atrial fibrillation and dialysis dependence; D5 with prior-only coded diagnoses. ----------
D3p <- sub("\\+ af \\+ eg \\+ esrd", "+ af_prior + eg + esrd_prior", D3); D5p <- paste(D3p, "+ dm_prior + htn_prior + cad_prior"); stopifnot(D3p != D3)
KEYp <- c(KEY, "af_prior","esrd_prior","dm_prior","htn_prior","cad_prior")
dxl <- list()
for (s in 0:1) for (def in c("first","confirmed")) { ivs <- lapply(covs, function(cv) build_iv(cv, s, def))
  fits <- lapply(ivs, function(x) fit_glm(x, D3p)); out <- pool_hr(fits, KEYp); out$stratum <- lab_st(s); out$def <- def; out$domain <- "D3 prior-only AF and dialysis"; out$n_pt <- fits[[1]]$n_pt; out$events <- fits[[1]]$events; dxl[[length(dxl)+1]] <- out
  fits <- lapply(ivs, function(x) fit_glm(x[x$yr1==1,], D5p)); out <- pool_hr(fits, KEYp); out$stratum <- lab_st(s); out$def <- def; out$domain <- "D5 prior-only coded diagnoses"; out$n_pt <- fits[[1]]$n_pt; out$events <- fits[[1]]$events; dxl[[length(dxl)+1]] <- out
  msg("DIAGNOSIS TIMING (completed admissions only) %s %s done", lab_st(s), def) }
res$dx_timing <- do.call(rbind, dxl)
# ---------- 5. SENSITIVITIES on PE D3, MI-pooled: visit-process inverse-intensity weights (IIW); inverse probability of censoring weights (IPCW) ----------
qw <- function(w) c(min=min(w), p5=unname(quantile(w,.05)), p25=unname(quantile(w,.25)), median=median(w), p75=unname(quantile(w,.75)), p95=unname(quantile(w,.95)), max=max(w))
sens <- list(); wdiag <- list()
for (s in 0:1) for (def in c("first","confirmed")) {
  fits_iiw <- list(); fits_ipcw <- list()
  for (m in seq_len(M_USE)) { cv <- covs[[m]]; iv <- build_iv(cv, s, def)
    # IIW: intensity of the next echocardiogram modelled on all observation intervals of the stratum (Andersen-Gill Cox model on time
    # since index) with the baseline covariates and the history known at the start of each interval (previous reported grade, whether the
    # previous episode was inpatient); stabilised weight for each outcome interval = exp(mean linear predictor - linear predictor), truncated
    ids <- cv$pid[cv$sev1n==s]; vp <- p[p$pid %in% ids & p$rn>1 & p$t>p$t_prev, c("pid","t_prev","t","sev_prev","inpt_prev")]
    vp <- merge(vp, cv, by="pid"); vp$one <- 1L
    vm <- coxph(as.formula(paste("Surv(t_prev/365.25, t/365.25, one) ~", D3, "+ sev_prev + inpt_prev")), data=vp)
    lp_iv <- predict(vm, newdata=iv, type="lp"); lp_ref <- mean(predict(vm, newdata=vp, type="lp"))
    w_raw <- exp(lp_ref - lp_iv); iv$w_iiw <- pmin(pmax(w_raw, 0.1), 10); fits_iiw[[m]] <- fit_glm(iv, D3, wts=iv$w_iiw)
    # IPCW: censoring = alive at the last study without the event; Cox model for censoring on the baseline covariates; weight = S_c(marginal)/S_c(individual) at the interval end, truncated
    x <- cv[cv$sev1n==s,]; tp <- if (def=="first") x$t_first else x$t_conf; x$cens <- as.integer(is.na(tp) & is.na(x$t_death)); x$tc <- pmax(ifelse(!is.na(tp), tp, x$t_last_study),1)/365.25
    cm <- coxph(as.formula(paste("Surv(tc, cens) ~", D3)), data=x); bh <- basehaz(cm, centered=FALSE); lp <- predict(cm, newdata=x, type="lp"); names(lp) <- x$pid
    iv$H0 <- approx(bh$time, bh$hazard, xout=iv$t/365.25, rule=2)$y; iv$Sc <- exp(-iv$H0*exp(lp[iv$pid])); iv$Sm <- exp(-iv$H0*exp(mean(lp))); w2_raw <- iv$Sm/iv$Sc; iv$w_ipcw <- pmin(pmax(w2_raw, 0.1), 10)
    fits_ipcw[[m]] <- fit_glm(iv, D3, wts=iv$w_ipcw)
    if (m==1) wdiag[[length(wdiag)+1]] <- rbind(data.frame(stratum=lab_st(s), def=def, weight="IIW", t(qw(w_raw)), pct_truncated=100*mean(w_raw<0.1 | w_raw>10), n_intervals=nrow(iv)),
                                               data.frame(stratum=lab_st(s), def=def, weight="IPCW", t(qw(w2_raw)), pct_truncated=100*mean(w2_raw<0.1 | w2_raw>10), n_intervals=nrow(iv))) }
  for (nm in c("IIW","IPCW")) { f <- if (nm=="IIW") fits_iiw else fits_ipcw; out <- pool_hr(f, KEY); out$stratum <- lab_st(s); out$def <- def; out$sens <- nm; sens[[length(sens)+1]] <- out }
  msg("SENS %s %s done", lab_st(s), def) }
res$sens <- do.call(rbind, sens); res$weight_diag <- do.call(rbind, wdiag)
# ---------- 6. TOTAL vs DIRECT (D1+x vs D3), PE, MI-pooled ----------
tvd <- list()
for (s in 0:1) for (def in c("first","confirmed")) { ivs <- lapply(covs, function(cv) build_iv(cv,s,def)); ivy <- lapply(ivs, function(iv) iv[iv$yr1==1,])
  for (v in c("zbmi","af","dm","htn","cad")) { sub <- v %in% c("dm","htn","cad"); use <- if (sub) ivy else ivs
    ft <- lapply(use, function(iv) fit_glm(iv, paste(D1,"+",v))); fd <- lapply(use, function(iv) fit_glm(iv, if (v %in% c("zbmi","af")) D3 else paste(D3,"+",v)))
    a <- pool_hr(ft, v); b <- pool_hr(fd, v); tvd[[length(tvd)+1]] <- data.frame(stratum=lab_st(s), def=def, term=v, total_hr=a$hr, total_lo=a$lo, total_hi=a$hi, direct_hr=b$hr, direct_lo=b$lo, direct_hi=b$hi) } }
res$total_direct <- do.call(rbind, tvd)
# ---------- 7. FINE-GRAY subdistribution hazards, pooled across imputations (finegray expansion + weighted Cox, patient-clustered), D3 ----------
fgv <- c("age10","male","era_n","setting","av_cat","ef_cat","zE","zla","zivs","zbmi","af","eg","esrd"); fg <- list()
for (s in 0:1) for (def in c("first","confirmed")) {
  fits <- lapply(covs, function(cv){ x <- cv[cv$sev1n==s,]; tp <- if (def=="first") x$t_first else x$t_conf; z <- mk(tp, x$t_death, x$t_last_study)
    xx <- x[, c("pid", fgv)]; xx$tt <- z$tt; xx$ev <- factor(z$ev, levels=0:2, labels=c("censor","event","death"))
    fgd <- finegray(Surv(tt, ev) ~ ., data=xx, etype="event")
    f <- coxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", D3)), data=fgd, weights=fgwt, cluster=pid)
    list(coef=coef(f), var=diag(vcov(f)), n_pt=nrow(x), events=sum(z$ev==1), deaths=sum(z$ev==2)) })
  out <- pool_hr(fits, KEY); names(out)[names(out)=="hr"] <- "shr"; out$stratum <- lab_st(s); out$def <- def; out$n <- fits[[1]]$n_pt; out$events <- fits[[1]]$events; out$deaths <- fits[[1]]$deaths; fg[[length(fg)+1]] <- out
  msg("FINE-GRAY %s %s pooled: n %d events %d deaths %d", lab_st(s), def, out$n[1], out$events[1], out$deaths[1]) }
res$finegray <- do.call(rbind, fg)
# ---------- 8. SENSITIVITIES: kidney-function definition, phosphate definition, care setting at index (PE D3/D4, MI-pooled) ----------
# eGFR variants use the median and the most recent creatinine of the year before index (protocol sensitivity definitions),
# fitted among patients with a measured creatinine so that the three definitions describe the same patients.
cutg <- function(x) relevel(cut(x, c(-Inf,30,60,Inf), labels=c("lt30","30to60","ge60")), "ge60")
d$eg_med <- cutg(d$egfr_med); d$eg_last <- cutg(d$egfr_last); d$zphos_last <- as.numeric(scale(d$phos_last)); measured_cr <- d$pid[!is.na(d$cr_min)]
lab_def <- c(eg="lowest creatinine (main)", eg_med="median creatinine", eg_last="most recent creatinine"); labsens <- list()
IVS <- list(); for (s in 0:1) for (def in c("first","confirmed")) IVS[[paste(s,def)]] <- lapply(covs, function(cv) build_iv(cv, s, def))
for (s in 0:1) for (def in c("first","confirmed")) for (v in names(lab_def)) {
  fits <- lapply(IVS[[paste(s,def)]], function(iv){ iv <- iv[iv$pid %in% measured_cr,]; if (v != "eg") iv$eg <- d[[v]][match(iv$pid, d$pid)]; fit_glm(iv, D3) })
  out <- pool_hr(fits, KEY); out$stratum <- lab_st(s); out$def <- def; out$definition <- lab_def[[v]]; out$n_pt <- fits[[1]]$n_pt; out$events <- fits[[1]]$events
  labsens[[length(labsens)+1]] <- out; msg("lab definition %-24s %-11s %-9s: patients %d events %d", lab_def[[v]], lab_st(s), def, out$n_pt[1], out$events[1]) }
res$sens_lab <- do.call(rbind, labsens)
phossens <- list()
for (s in 0:1) for (def in c("first","confirmed")) for (v in c("median of the year (main)","most recent value")) {
  fits <- lapply(IVS[[paste(s,def)]], function(iv){ if (v == "most recent value") iv$zphos <- d$zphos_last[match(iv$pid, d$pid)]; iv <- iv[!is.na(iv$zphos) & !is.na(iv$zca) & !is.na(iv$zalp),]; fit_glm(iv, D4) })
  out <- pool_hr(fits, c("zphos","zca","zalp")); out$stratum <- lab_st(s); out$def <- def; out$definition <- v; out$n_pt <- fits[[1]]$n_pt; out$events <- fits[[1]]$events
  phossens[[length(phossens)+1]] <- out }
res$sens_phos <- do.call(rbind, phossens)
# main model repeated by care setting at index (outpatient versus inpatient index study); the setting covariate is dropped within subgroups
D3_noset <- sub("\\+ setting ", "", D3); setsens <- list()
for (s in 0:1) for (def in c("first","confirmed")) for (g in 0:1) {
  fits <- lapply(IVS[[paste(s,def)]], function(iv){ iv <- iv[iv$idx_inpt==g,]; fit_glm(iv, D3_noset) })
  out <- pool_hr(fits, KEY); out$stratum <- lab_st(s); out$def <- def; out$setting <- ifelse(g==0,"outpatient index","inpatient index"); out$n_pt <- fits[[1]]$n_pt; out$events <- fits[[1]]$events
  setsens[[length(setsens)+1]] <- out; msg("care setting %-16s %-11s %-9s: patients %d events %d", out$setting[1], lab_st(s), def, out$n_pt[1], out$events[1]) }
res$sens_setting <- do.call(rbind, setsens)
saveRDS(res, "outputs/results/results_part1.rds"); msg("results_part1.rds saved")
