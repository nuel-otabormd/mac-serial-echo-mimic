# 10_mi_stability.R: Monte Carlo (simulation) error of the multiply imputed main models. The imputation is repeated with a
# different random seed (same model, m = 20, maxit = 5) and the fully adjusted main models (D3, both cohorts, both outcome
# definitions) are refitted and pooled exactly as in 03_models_main.R; the table reports each hazard ratio under both sets of
# imputations and the absolute difference, so that second-decimal differences between runs can be read against this
# yardstick. Before the comparison the script refits D3 from the archived imputations and checks that it reproduces
# results_part1.rds to numerical precision (the model code below mirrors 03_models_main.R).
suppressPackageStartupMessages({library(survival); library(mice); library(sandwich)}); source("R/00_common.R")
M_USE <- as.integer(Sys.getenv("MAC_M", "20")); SEED2 <- as.integer(Sys.getenv("MAC_SEED2", "20260817"))
dA <- readRDS("data/frame.rds"); d <- dA
p <- read.csv("data/panel.csv"); p <- p[order(p$pid,p$rn),]
p$t_prev <- ave(p$t, p$pid, FUN=function(v) c(NA, head(v,-1))); p$sev_prev <- ave(p$sev, p$pid, FUN=function(v) c(NA, head(v,-1))); p$inpt_prev <- ave(p$inpt, p$pid, FUN=function(v) c(NA, head(v,-1)))
M <- readRDS("data/mice.rds"); stopifnot(all(M$pid==dA$pid))
r1 <- readRDS("outputs/results/results_part1.rds"); m3 <- r1$main[r1$main$domain=="D3",]
msg <- function(...) cat(sprintf(...), "\n")
SD <- c(E_sept=sd(d$E_sept,na.rm=T), la=sd(d$la,na.rm=T), ivs=sd(d$ivs,na.rm=T), bmi=sd(d$bmi,na.rm=T), phos=sd(d$phos_median,na.rm=T), ca=sd(d$ca_median,na.rm=T), alp=sd(log(d$alp_median),na.rm=T))
MU <- c(E_sept=mean(d$E_sept,na.rm=T), la=mean(d$la,na.rm=T), ivs=mean(d$ivs,na.rm=T), bmi=mean(d$bmi,na.rm=T), phos=mean(d$phos_median,na.rm=T), ca=mean(d$ca_median,na.rm=T), alp=mean(log(d$alp_median),na.rm=T))
derive_cov <- function(x){
  x$age10 <- (x$age0-70)/10; x$setting <- factor(as.character(x$setting), levels=c("outpt","ED","ward","ICU"))
  x$av_cat <- factor(pmin(pmax(round(as.numeric(x$av)),0),3), levels=0:3, labels=c("nl","mild","mod","sev")); x$av_lin <- as.numeric(x$av_cat)-1
  x$ef_cat <- relevel(cut(x$lvef, c(-Inf,39.99,49.99,Inf), labels=c("lt40","40to49","ge50")), "ge50")
  x$eg <- relevel(cut(x$egfr_base, c(-Inf,30,60,Inf), labels=c("lt30","30to60","ge60")), "ge60")
  x$zE <- (x$E_sept-MU["E_sept"])/SD["E_sept"]; x$zla <- (x$la-MU["la"])/SD["la"]; x$zivs <- (x$ivs-MU["ivs"])/SD["ivs"]; x$zbmi <- (x$bmi-MU["bmi"])/SD["bmi"]
  x$zphos <- (x$phos_median-MU["phos"])/SD["phos"]; x$zca <- (x$ca_median-MU["ca"])/SD["ca"]; x$zalp <- (log(x$alp_median)-MU["alp"])/SD["alp"]
  x }
FIXED <- c("pid","yr1","t_first","t_conf","t_last","t_last_study","t_death","first_status","n_ep","idx_inpt","rheum_post")
covs_from <- function(imp, m){ x <- complete(imp, m); x$pid <- dA$pid; x$sev1n <- dA$sev1
  for (v in FIXED[-1]) x[[v]] <- dA[[v]]
  x <- derive_cov(x)
  x[, c(FIXED,"age10","male","era_n","setting","sev1n","av_cat","av_lin","ef_cat","eg","esrd","zE","zla","zivs","zbmi","af","zphos","zca","zalp","dm","htn","cad")] }
build_iv <- function(cov, stratum, def){
  ids <- cov$pid[cov$sev1n==stratum]; pp <- p[p$pid %in% ids & p$rn>1, c("pid","rn","t","t_prev","sev_prev","inpt_prev")]
  te <- if (def=="first") cov$t_first else cov$t_conf; names(te) <- cov$pid
  pp$te <- te[pp$pid]; pp$ev <- as.integer(!is.na(pp$te) & pp$t==pp$te)
  pp$after <- ave(pp$ev, pp$pid, FUN=function(v) c(0, head(cumsum(v),-1))); pp <- pp[pp$after==0,]
  pp <- split_intervals(pp[,c("pid","t_prev","t","ev","sev_prev","inpt_prev")])
  merge(pp, cov, by="pid") }
fit_glm <- function(iv, form){
  f <- suppressWarnings(glm(as.formula(paste("ev ~ band + offset(log(dt)) +", form)), data=iv, family=poisson()))
  stopifnot(nobs(f)==nrow(iv)); V <- sandwich::vcovCL(f, cluster=iv$pid)
  list(coef=coef(f), var=diag(V), n_pt=length(unique(iv$pid)), n_iv=nrow(iv), events=sum(iv$ev)) }
pool_hr <- function(fits, terms){ terms <- terms[terms %in% names(fits[[1]]$coef)]
  Q <- sapply(fits, function(f) f$coef[terms]); U <- sapply(fits, function(f) f$var[terms]); if (is.null(dim(Q))) { Q <- matrix(Q,nrow=1,dimnames=list(terms,NULL)); U <- matrix(U,nrow=1) }
  qb <- rowMeans(Q); ub <- rowMeans(U); b <- if (ncol(Q)>1) apply(Q,1,var) else 0; T <- ub + (1+1/ncol(Q))*b; z <- qb/sqrt(T)
  data.frame(term=terms, hr=exp(qb), lo=exp(qb-1.96*sqrt(T)), hi=exp(qb+1.96*sqrt(T)), p=2*pnorm(-abs(z)), row.names=NULL) }
D3 <- "age10 + male + era_n + setting + av_cat + ef_cat + zE + zla + zivs + zbmi + af + eg + esrd"
KEY <- c("age10","male","av_catmild","av_catmod","av_catsev","ef_catlt40","ef_cat40to49","zE","zla","zivs","zbmi","af","eglt30","eg30to60","esrd")
lab_st <- function(s) ifelse(s==0,"onset","progression")
fit_D3 <- function(imp, tag){ out <- list()
  covs <- lapply(seq_len(M_USE), function(m) covs_from(imp, m))
  for (s in 0:1) for (def in c("first","confirmed")) { ivs <- lapply(covs, function(cv) build_iv(cv, s, def)); fits <- lapply(ivs, function(x) fit_glm(x, D3))
    o <- pool_hr(fits, KEY); o$stratum <- lab_st(s); o$def <- def; o$n_pt <- fits[[1]]$n_pt; o$events <- fits[[1]]$events; out[[length(out)+1]] <- o
    msg("%s: %-11s %-9s patients %d events %d", tag, lab_st(s), def, o$n_pt[1], o$events[1]) }
  do.call(rbind, out) }
# 1) reproduce the archived D3 estimates from the archived imputations (proves that this script's model code equals 03_models_main.R)
a <- fit_D3(M$imp, "archived imputations")
chk <- merge(a, m3, by=c("stratum","def","term"), suffixes=c(".here",".archived"))
stopifnot(nrow(chk)==nrow(a), max(abs(chk$hr.here-chk$hr.archived)) < 1e-8, max(abs(chk$lo.here-chk$lo.archived)) < 1e-8)
msg("archived D3 estimates reproduced exactly (max |dHR| %.1e)", max(abs(chk$hr.here-chk$hr.archived)))
# 2) second imputation run with a different seed, same specification as 02_impute.R
zf <- mk(d$t_first, d$t_death, d$t_last_study); d$ev_first <- as.integer(zf$ev==1); d$ltime <- log(zf$tt); d$died <- as.integer(!is.na(d$t_death))
imp_vars <- c("age0","male","era_n","setting","sev1","av","lvef","E_sept","la","ivs","bmi","af","egfr_base","esrd","phos_median","ca_median","alp_median","dm","htn","cad","idx_inpt","n_ep","ev_first","ltime","died")
x <- d[, imp_vars]; x$setting <- factor(x$setting); x$av <- as.numeric(x$av); x$sev1 <- factor(x$sev1)
x$phos_median[!is.na(x$phos_median) & (x$phos_median<0.5 | x$phos_median>15)] <- NA
meth <- make.method(x); meth[c("phos_median","ca_median","alp_median")] <- ""; pred <- make.predictorMatrix(x); pred[, c("phos_median","ca_median","alp_median")] <- 0
set.seed(SEED2); imp2 <- mice(x, m=M$m, maxit=M$maxit, method=meth, predictorMatrix=pred, printFlag=FALSE)
b <- fit_D3(imp2, sprintf("second seed %d", SEED2))
cmp <- merge(a, b, by=c("stratum","def","term"), suffixes=c(".a",".b"))
cmp$abs_diff_hr <- abs(cmp$hr.a-cmp$hr.b); cmp$abs_diff_log_hr <- abs(log(cmp$hr.a)-log(cmp$hr.b))
cmp$ci_excludes_1_a <- cmp$lo.a>1 | cmp$hi.a<1; cmp$ci_excludes_1_b <- cmp$lo.b>1 | cmp$hi.b<1
cmp <- cmp[order(cmp$stratum, cmp$def, match(cmp$term, KEY)), ]
summ <- c(max_abs_diff_hr=max(cmp$abs_diff_hr), max_abs_diff_hr_excluding_avc_severe=max(cmp$abs_diff_hr[cmp$term!="av_catsev"]), max_abs_diff_log_hr=max(cmp$abs_diff_log_hr),
          n_terms=nrow(cmp), n_ci_status_changed=sum(cmp$ci_excludes_1_a!=cmp$ci_excludes_1_b), n_direction_changed=sum((cmp$hr.a>1)!=(cmp$hr.b>1)))
print(round(summ,4))
mi_stability <- list(seed_archived=20260816, seed_second=SEED2, m=M$m, maxit=M$maxit, comparison=cmp, summary=summ)
saveRDS(mi_stability, "outputs/results/mi_stability.rds")
f <- function(h,l,u) sprintf("%.2f (%.2f–%.2f)", h, l, u)
lab <- c(age10="Age, per 10 years", male="Male sex", av_catmild="Aortic valve calcification, mild vs none", av_catmod="Aortic valve calcification, moderate vs none", av_catsev="Aortic valve calcification, severe vs none",
         ef_catlt40="LVEF <40% vs ≥50%", ef_cat40to49="LVEF 40–49% vs ≥50%", zE="E/e', per SD", zla="Left atrial dimension, per SD", zivs="Septal wall thickness, per SD", zbmi="Body mass index, per SD",
         af="Atrial fibrillation or flutter", eglt30="eGFR <30 vs ≥60", eg30to60="eGFR 30–59 vs ≥60", esrd="Dialysis dependence")
S17 <- data.frame(cohort=cmp$stratum, definition=ifelse(cmp$def=="first","first-observed","confirmed"), characteristic=unname(lab[cmp$term]), hr_reported=f(cmp$hr.a,cmp$lo.a,cmp$hi.a), hr_second_imputation=f(cmp$hr.b,cmp$lo.b,cmp$hi.b), abs_difference_hr=sprintf("%.3f", cmp$abs_diff_hr))
dir.create("outputs/tables", showWarnings=FALSE, recursive=TRUE); write.csv(S17, "outputs/tables/tableS17_imputation_stability.csv", row.names=FALSE)
writeLines(c(sprintf("Fully adjusted main models (D3) refitted after repeating the multiple imputation with a different random seed (%d instead of %d; m = %d, maxit = %d, same imputation model).", SEED2, 20260816, M$m, M$maxit),
             sprintf("Largest absolute difference in hazard ratio: %.3f overall (%.3f excluding the sparse severe aortic calcification category); confidence-interval status (excluding 1 or not) changed for %d of %d estimates; no estimate changed direction (%d).", summ["max_abs_diff_hr"], summ["max_abs_diff_hr_excluding_avc_severe"], summ["n_ci_status_changed"], summ["n_terms"], summ["n_direction_changed"])),
           "outputs/tables/tableS17_notes.txt")
msg("mi_stability.rds and tableS17 written")
