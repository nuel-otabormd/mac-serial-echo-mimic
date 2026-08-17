# 03_models_main.R: co-primary models (piecewise-exponential on person-intervals, MI-pooled), imputation bridge,
# stratum interaction, stricter onset cohort, informative-observation and censoring sensitivities, total versus direct
# associations, Fine and Gray companions, and the laboratory-definition and care-setting sensitivities.
suppressPackageStartupMessages({library(survival); library(mice); library(cmprsk)})
set.seed(20260816)
d <- readRDS("data/frame.rds"); p <- read.csv("data/panel.csv"); M <- readRDS("data/mice.rds"); imp <- M$imp; stopifnot(all(M$pid==d$pid))
p <- p[order(p$pid,p$rn),]; p$t_prev <- ave(p$t, p$pid, FUN=function(v) c(NA, head(v,-1)))
res <- list(); msg <- function(...) cat(sprintf(...), "\n")
# ---------- standardisation constants from observed data ----------
SD <- c(E_sept=sd(d$E_sept,na.rm=T), la=sd(d$la,na.rm=T), ivs=sd(d$ivs,na.rm=T), bmi=sd(d$bmi,na.rm=T), phos=sd(d$phos_median,na.rm=T), ca=sd(d$ca_median,na.rm=T), alp=sd(log(d$alp_median),na.rm=T))
MU <- c(E_sept=mean(d$E_sept,na.rm=T), la=mean(d$la,na.rm=T), ivs=mean(d$ivs,na.rm=T), bmi=mean(d$bmi,na.rm=T), phos=mean(d$phos_median,na.rm=T), ca=mean(d$ca_median,na.rm=T), alp=mean(log(d$alp_median),na.rm=T))
covs_m <- function(m){ x <- complete(imp, m); x$pid <- d$pid
  x$age10 <- (x$age0-70)/10; x$era_n <- x$era_n; x$setting <- factor(as.character(x$setting), levels=c("outpt","ED","ward","ICU"))
  x$av_cat <- factor(pmin(pmax(round(x$av),0),3), levels=0:3, labels=c("nl","mild","mod","sev")); x$av_lin <- as.numeric(x$av_cat)-1
  x$ef_cat <- relevel(cut(x$lvef, c(-Inf,39.99,49.99,Inf), labels=c("lt40","40to49","ge50")), "ge50")
  x$eg <- relevel(cut(x$egfr_base, c(-Inf,30,60,Inf), labels=c("lt30","30to60","ge60")), "ge60")
  x$zE <- (x$E_sept-MU["E_sept"])/SD["E_sept"]; x$zla <- (x$la-MU["la"])/SD["la"]; x$zivs <- (x$ivs-MU["ivs"])/SD["ivs"]; x$zbmi <- (x$bmi-MU["bmi"])/SD["bmi"]
  x$zphos <- (x$phos_median-MU["phos"])/SD["phos"]; x$zca <- (x$ca_median-MU["ca"])/SD["ca"]; x$zalp <- (log(x$alp_median)-MU["alp"])/SD["alp"]
  x$yr1 <- d$yr1; x$t_first <- d$t_first; x$t_conf <- d$t_conf; x$t_last <- d$t_last; x$t_death <- d$t_death; x$first_status <- d$first_status; x$sev1n <- d$sev1; x$n_ep <- d$n_ep
  x[, c("pid","age10","male","era_n","setting","sev1n","av_cat","av_lin","ef_cat","eg","esrd","zE","zla","zivs","zbmi","af","zphos","zca","zalp","dm","htn","cad","yr1","t_first","t_conf","t_last","t_death","first_status","n_ep","idx_inpt")] }
# ---------- person-interval builder ----------
build_iv <- function(cov, stratum, def, from_rn=1, t_event_override=NULL){
  ids <- cov$pid[cov$sev1n==stratum]; pp <- p[p$pid %in% ids & p$rn>from_rn,]
  if (from_rn>1) { t0 <- p$t[p$rn==from_rn]; names(t0) <- p$pid[p$rn==from_rn]; pp$t_prev <- pmax(pp$t_prev, t0[pp$pid]) }
  te <- if (!is.null(t_event_override)) t_event_override else if (def=="first") cov$t_first else cov$t_conf; names(te) <- cov$pid
  pp$te <- te[pp$pid]; pp$ev <- as.integer(!is.na(pp$te) & pp$t==pp$te)
  pp$after <- ave(pp$ev, pp$pid, FUN=function(v) c(0, head(cumsum(v),-1))); pp <- pp[pp$after==0,]   # drop intervals after the event
  if (from_rn>1) { pp <- pp[pp$t > pp$t_prev,]; pp <- pp[is.na(pp$te) | pp$te > (t0[pp$pid]),] }      # stricter baseline: events before/at 2nd episode not counted (they define it)
  pp$dt <- pmax(pp$t-pp$t_prev,1)/365.25; pp$band <- cut(pp$t_prev/365.25, c(-Inf,0.5,1,2,4,7,Inf), labels=c("a","b","c","d","e","f"))
  merge(pp[,c("pid","t_prev","t","dt","band","ev")], cov, by="pid") }
pe_fit <- function(iv, form, wts=NULL){ glm(as.formula(paste("ev ~ band + offset(log(dt)) +", form)), data=iv, family=poisson, weights=wts) }
pool_hr <- function(fits, terms){ Q <- sapply(fits, function(f) coef(f)[terms]); U <- sapply(fits, function(f) diag(vcov(f))[terms]); if (is.null(dim(Q))) { Q <- matrix(Q,nrow=1,dimnames=list(terms,NULL)); U <- matrix(U,nrow=1) }
  qb <- rowMeans(Q); ub <- rowMeans(U); b <- apply(Q,1,var); T <- ub + (1+1/length(fits))*b; z <- qb/sqrt(T)
  data.frame(term=terms, hr=exp(qb), lo=exp(qb-1.96*sqrt(T)), hi=exp(qb+1.96*sqrt(T)), p=2*pnorm(-abs(z)), row.names=NULL) }
D1 <- "age10 + male + era_n + setting"; D2 <- paste(D1, "+ av_cat + ef_cat + zE + zla + zivs + zbmi + af"); D3 <- paste(D2, "+ eg + esrd"); D4 <- paste(D3, "+ zphos + zca + zalp"); D5 <- paste(D3, "+ dm + htn + cad")
KEY <- c("age10","male","av_catmild","av_catmod","av_catsev","ef_catlt40","ef_cat40to49","zE","zla","zivs","zbmi","af","eglt30","eg30to60","esrd","zphos","zca","zalp","dm","htn","cad")
covs <- lapply(1:20, covs_m)
# ---------- 1. MAIN MODELS: piecewise-exponential (primary), MI-pooled, both strata x definitions x domains ----------
main <- list()
for (s in 0:1) for (def in c("first","confirmed")) {
  ivs <- lapply(covs, function(cv) build_iv(cv, s, def))
  n_pt <- length(unique(ivs[[1]]$pid)); n_iv <- nrow(ivs[[1]]); n_ev <- sum(ivs[[1]]$ev)
  for (dom in c("D1","D2","D3","D4","D5")) { form <- get(dom)
    ivd <- if (dom=="D4") lapply(ivs, function(x) x[!is.na(x$zphos),]) else if (dom=="D5") lapply(ivs, function(x) x[x$yr1==1,]) else ivs
    fits <- lapply(ivd, function(x) pe_fit(x, form)); tm <- KEY[KEY %in% names(coef(fits[[1]]))]
    out <- pool_hr(fits, tm); out$stratum <- ifelse(s==0,"onset","progression"); out$def <- def; out$domain <- dom; out$n_pt <- length(unique(ivd[[1]]$pid)); out$n_iv <- nrow(ivd[[1]]); out$events <- sum(ivd[[1]]$ev)
    main[[length(main)+1]] <- out
    msg("PE %-11s %-9s %s: patients %d intervals %d events %d", out$stratum[1], def, dom, out$n_pt[1], out$n_iv[1], out$events[1]) }
  # Cox companion (event at qualifying date), MI-pooled, D3
  cx <- lapply(covs, function(cv){ x <- cv[cv$sev1n==s,]; tp <- if (def=="first") x$t_first else x$t_conf; ev <- ifelse(!is.na(tp),1, ifelse(!is.na(x$t_death) & (is.na(x$t_last) | x$t_death<=x$t_last+30),2,0)); tt <- pmax(ifelse(!is.na(tp),tp, ifelse(ev==2, x$t_death, x$t_last)),1)/365.25
    coxph(as.formula(paste("Surv(tt, ev==1) ~", D3)), data=cbind(x,tt=tt,ev=ev)) })
  out <- pool_hr(cx, KEY[KEY %in% names(coef(cx[[1]]))]); out$stratum <- ifelse(s==0,"onset","progression"); out$def <- def; out$domain <- "D3_cox_companion"; out$n_pt <- cx[[1]]$n; out$n_iv <- NA; out$events <- cx[[1]]$nevent; main[[length(main)+1]] <- out
}
res$main <- do.call(rbind, main)
# ---------- 2. MI BRIDGE (confirmation imputed for never-re-examined reads), joint with covariate imputation, PE D3 ----------
bridge <- list()
for (s in 0:1) { fits <- list()
  for (m in 1:20) { cv <- covs[[m]]; x <- cv[cv$sev1n==s,]; x$t_end <- ifelse(!is.na(x$t_death), x$t_death, x$t_last); x$t_after <- (x$t_end-x$t_first)/365.25; x$died <- !is.na(x$t_death)
    re <- x[x$first_status %in% c("confirmed","refuted"),]; un <- x[x$first_status=="unconfirmable",]
    im <- glm(I(first_status=="confirmed") ~ age10 + male + setting + era_n + ef_cat + esrd + av_lin + t_after + died, data=re, family=binomial)
    pun <- predict(im, newdata=un, type="response"); hit <- un$pid[rbinom(nrow(un),1,pun)==1]
    tconf <- x$t_conf; tconf[x$pid %in% hit] <- x$t_first[x$pid %in% hit]; names(tconf) <- x$pid
    full <- rep(NA_real_, nrow(cv)); names(full) <- cv$pid; full[names(tconf)] <- tconf
    iv <- build_iv(cv, s, "confirmed", t_event_override=full); fits[[m]] <- pe_fit(iv, D3); if (m==1) nev <- sum(iv$ev) }
  out <- pool_hr(fits, KEY[KEY %in% names(coef(fits[[1]]))]); out$stratum <- ifelse(s==0,"onset","progression"); out$events_imp1 <- nev; bridge[[s+1]] <- out; msg("BRIDGE %s events(imp1) %d", out$stratum[1], nev) }
res$bridge <- do.call(rbind, bridge)
# ---------- 3. INTERACTION (first primary; confirmed check), PE D3 with stratum interaction, MI-pooled ----------
inter <- list()
for (def in c("first","confirmed")) { fits <- lapply(covs, function(cv){ iv <- rbind(build_iv(cv,0,def), build_iv(cv,1,def)); iv$prog <- as.integer(iv$sev1n==1)
    glm(as.formula(paste("ev ~ band*prog + offset(log(dt)) + (age10 + male + av_lin + ef_cat + zE + zla + zivs + eg + esrd + af + zbmi)*prog + era_n + setting")), data=iv, family=poisson) })
  it <- grep(":prog$|^prog:", names(coef(fits[[1]])), value=TRUE); it <- it[!grepl("band", it)]
  out <- pool_hr(fits, it); out$def <- def; inter[[def]] <- out }
res$interaction <- do.call(rbind, inter)
# ---------- 4. CONFIRMED MAC-FREE BASELINE (index blank & 2nd blank; from 2nd episode), PE D3, MI-pooled ----------
sec <- p[p$rn==2, c("pid","sev")]; ok2 <- sec$pid[sec$sev==0]
cb <- list()
for (def in c("first","confirmed")) { fits <- lapply(covs, function(cv){ cv2 <- cv[cv$pid %in% ok2,]; iv <- build_iv(cv2, 0, def, from_rn=2); pe_fit(iv, D3) })
  iv1 <- build_iv(covs[[1]][covs[[1]]$pid %in% ok2,], 0, def, from_rn=2)
  out <- pool_hr(fits, KEY[KEY %in% names(coef(fits[[1]]))]); out$def <- def; out$n_pt <- length(unique(iv1$pid)); out$events <- sum(iv1$ev); cb[[def]] <- out; msg("CONFIRMED-MAC-FREE baseline %s: patients %d events %d", def, out$n_pt[1], out$events[1]) }
res$confirmed_baseline <- do.call(rbind, cb)
# ---------- 5. SENSITIVITIES on PE D3, MI-pooled: echo-rate covariate; IIW; IPCW ----------
sens <- list()
for (s in 0:1) for (def in c("first","confirmed")) {
  fits_rate <- list(); fits_iiw <- list(); fits_ipcw <- list()
  for (m in 1:20) { cv <- covs[[m]]; iv <- build_iv(cv, s, def)
    # echo-rate covariate
    cv$lrate <- log((cv$n_ep-1)/pmax(cv$t_last,1)*365.25); iv$lrate <- cv$lrate[match(iv$pid, cv$pid)]; fits_rate[[m]] <- pe_fit(iv, paste(D3, "+ lrate"))
    # IIW: Poisson model for episode count ~ covariates, offset log follow-up; stabilised weight = marginal/individual predicted rate
    x <- cv[cv$sev1n==s,]; x$nint <- x$n_ep-1; x$fu <- pmax(x$t_last,1)/365.25
    pm <- glm(as.formula(paste("nint ~ offset(log(fu)) +", D3)), data=x, family=poisson); x$rate_i <- predict(pm, type="response")/x$fu; x$w_iiw <- pmin(pmax(mean(x$rate_i)/x$rate_i, 0.1), 10)
    iv$w_iiw <- x$w_iiw[match(iv$pid, x$pid)]; fits_iiw[[m]] <- pe_fit(iv, D3, wts=iv$w_iiw)
    # IPCW: censoring = alive at last echo without event; Cox for censoring on covariates; weight = 1/S_c(t_interval_end)
    tp <- if (def=="first") x$t_first else x$t_conf; x$cens <- as.integer(is.na(tp) & is.na(x$t_death)); x$tc <- pmax(ifelse(!is.na(tp), tp, x$t_last),1)/365.25
    cm <- coxph(as.formula(paste("Surv(tc, cens) ~", D3)), data=x); bh <- basehaz(cm, centered=FALSE); lp <- predict(cm, newdata=x, type="lp"); names(lp) <- x$pid
    iv$H0 <- approx(bh$time, bh$hazard, xout=iv$t/365.25, rule=2)$y; iv$Sc <- exp(-iv$H0*exp(lp[iv$pid])); iv$Sm <- exp(-iv$H0*exp(mean(lp))); iv$w_ipcw <- pmin(pmax(iv$Sm/iv$Sc, 0.1), 10)
    fits_ipcw[[m]] <- pe_fit(iv, D3, wts=iv$w_ipcw) }
  for (nm in c("echo_rate","IIW","IPCW")) { f <- get(paste0("fits_", tolower(gsub("echo_rate","rate",nm)))); out <- pool_hr(f, KEY[KEY %in% names(coef(f[[1]]))]); out$stratum <- ifelse(s==0,"onset","progression"); out$def <- def; out$sens <- nm; sens[[length(sens)+1]] <- out }
  msg("SENS %s %s done", ifelse(s==0,"onset","progression"), def) }
res$sens <- do.call(rbind, sens)
# ---------- 6. TOTAL vs DIRECT (D1+x vs D3), PE, MI-pooled ----------
tvd <- list()
for (s in 0:1) for (def in c("first","confirmed")) for (v in c("zbmi","af","dm","htn","cad")) { sub <- v %in% c("dm","htn","cad")
  ft <- lapply(covs, function(cv){ iv <- build_iv(cv,s,def); if (sub) iv <- iv[iv$yr1==1,]; pe_fit(iv, paste(D1,"+",v)) }); fd <- lapply(covs, function(cv){ iv <- build_iv(cv,s,def); if (sub) iv <- iv[iv$yr1==1,]; pe_fit(iv, if (v %in% c("zbmi","af")) D3 else paste(D3,"+",v)) })
  a <- pool_hr(ft, v); b <- pool_hr(fd, v); tvd[[length(tvd)+1]] <- data.frame(stratum=ifelse(s==0,"onset","progression"), def=def, term=v, total_hr=a$hr, total_lo=a$lo, total_hi=a$hi, direct_hr=b$hr, direct_lo=b$lo, direct_hi=b$hi) }
res$total_direct <- do.call(rbind, tvd)
# ---------- 7. FINE-GRAY (complete-case, cause-specific comparison), D3 ----------
fg <- list()
for (s in 0:1) for (def in c("first","confirmed")) { x <- covs[[1]][covs[[1]]$sev1n==s,]  # imputation 1 covariates (FG on all imputations is slow); note in paper
  tp <- if (def=="first") x$t_first else x$t_conf; ev <- ifelse(!is.na(tp),1, ifelse(!is.na(x$t_death) & (is.na(x$t_last) | x$t_death<=x$t_last+30),2,0)); tt <- pmax(ifelse(!is.na(tp),tp, ifelse(ev==2, x$t_death, x$t_last)),1)/365.25
  X <- model.matrix(as.formula(paste("~", D3)), data=x)[,-1]; f <- crr(tt, ev, X, failcode=1, cencode=0); s_ <- summary(f)$conf.int
  fg[[length(fg)+1]] <- data.frame(stratum=ifelse(s==0,"onset","progression"), def=def, term=rownames(s_), shr=s_[,1], lo=s_[,3], hi=s_[,4], n=nrow(x), events=sum(ev==1), deaths=sum(ev==2)) }
res$finegray <- do.call(rbind, fg)
# ---------- 6. SENSITIVITIES: kidney-function definition, phosphate definition, care setting at index (PE D3/D4, MI-pooled) ----------
# eGFR variants use the median and the most recent creatinine of the year before index (protocol sensitivity definitions),
# fitted among patients with a measured creatinine so that the three definitions describe the same patients.
cutg <- function(x) relevel(cut(x, c(-Inf,30,60,Inf), labels=c("lt30","30to60","ge60")), "ge60")
d$eg_med <- cutg(d$egfr_med); d$eg_last <- cutg(d$egfr_last); d$zphos_last <- as.numeric(scale(d$phos_last)); measured_cr <- d$pid[!is.na(d$cr_min)]
lab_def <- c(eg="lowest creatinine (main)", eg_med="median creatinine", eg_last="most recent creatinine"); labsens <- list()
for (s in 0:1) for (def in c("first","confirmed")) for (v in names(lab_def)) {
  fits <- lapply(covs, function(cv){ iv <- build_iv(cv, s, def); iv <- iv[iv$pid %in% measured_cr,]
    if (v != "eg") iv$eg <- d[[v]][match(iv$pid, d$pid)]
    pe_fit(iv, D3) })
  out <- pool_hr(fits, KEY[KEY %in% names(coef(fits[[1]]))]); out$stratum <- ifelse(s==0,"onset","progression"); out$def <- def; out$definition <- lab_def[[v]]
  iv1 <- build_iv(covs[[1]], s, def); iv1 <- iv1[iv1$pid %in% measured_cr,]; out$n_pt <- length(unique(iv1$pid)); out$events <- sum(iv1$ev)
  labsens[[length(labsens)+1]] <- out; msg("lab definition %-24s %-11s %-9s: patients %d events %d", lab_def[[v]], out$stratum[1], def, out$n_pt[1], out$events[1]) }
res$sens_lab <- do.call(rbind, labsens)
phossens <- list()
for (s in 0:1) for (def in c("first","confirmed")) for (v in c("median of the year (main)","most recent value")) {
  fits <- lapply(covs, function(cv){ iv <- build_iv(cv, s, def); iv <- iv[!is.na(iv$zphos),]
    if (v == "most recent value") iv$zphos <- d$zphos_last[match(iv$pid, d$pid)]
    pe_fit(iv, D4) })
  out <- pool_hr(fits, c("zphos","zca","zalp")); out$stratum <- ifelse(s==0,"onset","progression"); out$def <- def; out$definition <- v
  iv1 <- build_iv(covs[[1]], s, def); iv1 <- iv1[!is.na(iv1$zphos),]; out$n_pt <- length(unique(iv1$pid)); out$events <- sum(iv1$ev)
  phossens[[length(phossens)+1]] <- out }
res$sens_phos <- do.call(rbind, phossens)
# main model repeated by care setting at index (outpatient versus inpatient index study); the setting covariate is dropped within subgroups
D3_noset <- sub("\\+ setting ", "", D3); setsens <- list()
for (s in 0:1) for (def in c("first","confirmed")) for (g in 0:1) {
  fits <- lapply(covs, function(cv){ iv <- build_iv(cv, s, def); iv <- iv[iv$idx_inpt==g,]; pe_fit(iv, D3_noset) })
  out <- pool_hr(fits, KEY[KEY %in% names(coef(fits[[1]]))]); out$stratum <- ifelse(s==0,"onset","progression"); out$def <- def; out$setting <- ifelse(g==0,"outpatient index","inpatient index")
  iv1 <- build_iv(covs[[1]], s, def); iv1 <- iv1[iv1$idx_inpt==g,]; out$n_pt <- length(unique(iv1$pid)); out$events <- sum(iv1$ev)
  setsens[[length(setsens)+1]] <- out; msg("care setting %-16s %-11s %-9s: patients %d events %d", out$setting[1], out$stratum[1], def, out$n_pt[1], out$events[1]) }
res$sens_setting <- do.call(rbind, setsens)
saveRDS(res, "outputs/results/results_part1.rds"); msg("results_part1.rds saved")
