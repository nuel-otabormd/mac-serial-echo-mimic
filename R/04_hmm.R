# 04_hmm.R: hidden Markov model of the reported grade (msm): underlying grade unobserved and non-regressing, death an exactly
# observed absorbing state, initial-state distribution estimated. Fitted by BFGS from five dispersed sets of starting values
# (in parallel); the fit with the lowest -2 log-likelihood among converged fits is kept, and the table of starts is saved so
# that convergence and the stability of the maximum can be reported. Diagnostics: observed versus model-expected prevalence
# of the reported grades over time.
suppressPackageStartupMessages({library(msm); library(parallel)}); source("R/00_common.R")
p <- read.csv("data/panel.csv"); p <- p[order(p$pid,p$rn),]; f <- readRDS("data/frame.rds")
p <- p[p$pid %in% f$pid,]                                            # primary cohort only
d <- p[,c("pid","t","sev")]; d$state <- d$sev+1L; d$t <- d$t/365.25
dth <- f[!is.na(f$t_death) & f$t_death > f$t_last, c("pid","t_death")]; dth <- data.frame(pid=dth$pid, t=dth$t_death/365.25, sev=NA, state=5L)
d <- rbind(d, dth); d <- d[order(d$pid,d$t),]; d <- d[!duplicated(d[c("pid","t")]),]   # two episodes on the same day contribute one observation time
n_obs <- table(d$pid); d <- d[d$pid %in% names(n_obs[n_obs>=2]),]
n_same_day <- length(unique(f$pid)) - length(unique(d$pid))
cat(sprintf("HMM data: %d patients, %d observations (%d patients whose only two episodes fell on the same date are excluded)\n", length(unique(d$pid)), nrow(d), n_same_day))
Q0 <- rbind(c(0,0.032,0,0,0.047), c(0,0,0.045,0,0.093), c(0,0,0,0.068,0.13), c(0,0,0,0,0.18), c(0,0,0,0,0))
E0 <- rbind(c(0.91,0.09,0.002,0.0005,0), c(0.39,0.57,0.04,0.001,0), c(0.07,0.33,0.59,0.01,0), c(0.06,0.02,0.48,0.43,0), c(0,0,0,0,1))
# five dispersed starting sets: the protocol values, faster/slower progression, less/more misclassification, and a shifted initial distribution
starts <- list(
  list(name="protocol values",           Q=Q0,       E=E0,                                                                          ip=c(0.67,0.23,0.08,0.02,0)),
  list(name="faster progression",        Q=Q0*2,     E=E0,                                                                          ip=c(0.67,0.23,0.08,0.02,0)),
  list(name="slower progression",        Q=Q0*0.5,   E=E0,                                                                          ip=c(0.67,0.23,0.08,0.02,0)),
  list(name="little misclassification",  Q=Q0,       E=rbind(c(0.97,0.03,0.001,0.0005,0), c(0.15,0.80,0.05,0.001,0), c(0.03,0.15,0.80,0.02,0), c(0.02,0.02,0.20,0.76,0), c(0,0,0,0,1)), ip=c(0.67,0.23,0.08,0.02,0)),
  list(name="more calcification at index", Q=Q0,     E=E0,                                                                          ip=c(0.55,0.30,0.12,0.03,0)))
fit_one <- function(st){ t0 <- Sys.time()
  m <- tryCatch(msm(state ~ t, subject=pid, data=d, qmatrix=st$Q, ematrix=st$E, deathexact=5, est.initprobs=TRUE, initprobs=st$ip, control=list(fnscale=90000,maxit=6000,reltol=1e-8), method="BFGS"), error=function(e) e)
  if (inherits(m, "error")) return(list(name=st$name, ok=FALSE, msg=conditionMessage(m), minutes=as.numeric(difftime(Sys.time(),t0,units="mins"))))
  list(name=st$name, ok=TRUE, model=m, converged=(m$opt$convergence==0), foundse=isTRUE(m$foundse), minus2LL=m$minus2loglik, minutes=as.numeric(difftime(Sys.time(),t0,units="mins"))) }
RNGkind("L'Ecuyer-CMRG"); set.seed(20260816)
nc <- suppressWarnings(as.integer(Sys.getenv("MAC_CORES", ""))); if (is.na(nc)) { nc <- detectCores(); nc <- if (is.na(nc)) 1L else max(1L, nc-1L) }   # portable: detectCores() can be NA in restricted sandboxes; MAC_CORES overrides
fits <- mclapply(starts, fit_one, mc.cores=min(5L, nc), mc.set.seed=TRUE)
tab <- do.call(rbind, lapply(fits, function(z) data.frame(start=z$name, fitted=z$ok, converged=if (z$ok) z$converged else NA, hessian_pd=if (z$ok) z$foundse else NA, minus2LL=if (z$ok) z$minus2LL else NA, minutes=round(z$minutes,1), note=if (z$ok) "" else z$msg)))
print(tab, row.names=FALSE)
ok <- which(tab$fitted & tab$converged %in% TRUE & tab$hessian_pd %in% TRUE)   # a fit is usable only if the optimiser converged AND the Hessian is positive definite
if (!length(ok)) stop("no starting set converged with a positive-definite Hessian; the hidden Markov model cannot be reported from this run")
best <- ok[which.min(tab$minus2LL[ok])]; mF <- fits[[best]]$model
tab$best <- seq_len(nrow(tab))==best
cat(sprintf("kept: start %d (%s), -2LL %.1f; spread of -2LL across converged starts: %.2f\n", best, tab$start[best], tab$minus2LL[best], diff(range(tab$minus2LL[ok]))))
cat("initial true-state probabilities:\n"); print(round(mF$hmodel$initprobs,3))
cat("emission (rows=true, cols=reported blank/mild/mod/sev):\n"); print(round(ematrix.msm(mF)$estimates[1:4,1:4],3))
cat("intensities /yr:\n"); print(round(qmatrix.msm(mF)$estimates[1:4,],4)); cat("sojourn yrs:\n"); print(sojourn.msm(mF))
# diagnostics: observed and model-expected prevalence of the reported grades and death at fixed times since index
prev <- tryCatch(prevalence.msm(mF, times=c(0,1,2,4,7,10)), error=function(e) NULL)
if (!is.null(prev)) { cat("observed vs expected prevalence (%):\n"); print(round(prev$"Observed percentages",1)); print(round(prev$"Expected percentages",1)) }
hmm_starts <- tab; hmm_prevalence <- prev; hmm_n <- c(patients=length(unique(d$pid)), observations=nrow(d), same_day_excluded=n_same_day)
save(mF, hmm_starts, hmm_prevalence, hmm_n, file="outputs/results/hmm.RData")
