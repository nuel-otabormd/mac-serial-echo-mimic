# 04_hmm.R: hidden Markov model of the reported grade (msm), underlying grade unobserved, death exact, initial
# state distribution estimated; deterministic BFGS fit from fixed starting values.
suppressPackageStartupMessages(library(msm))
p <- read.csv("data/panel.csv"); f <- read.csv("data/frame.csv")
d <- p[,c("pid","t","sev")]; d$state <- d$sev+1L; d$t <- d$t/365.25
dth <- f[!is.na(f$t_death) & f$t_death > f$t_last, c("pid","t_death")]; dth <- data.frame(pid=dth$pid, t=dth$t_death/365.25, sev=NA, state=5L)
d <- rbind(d, dth); d <- d[order(d$pid,d$t),]; d <- d[!duplicated(d[c("pid","t")]),]
n_obs <- table(d$pid); d <- d[d$pid %in% names(n_obs[n_obs>=2]),]
Q <- rbind(c(0,0.032,0,0,0.047), c(0,0,0.045,0,0.093), c(0,0,0,0.068,0.13), c(0,0,0,0,0.18), c(0,0,0,0,0))
E <- rbind(c(0.91,0.09,0.002,0.0005,0), c(0.39,0.57,0.04,0.001,0), c(0.07,0.33,0.59,0.01,0), c(0.06,0.02,0.48,0.43,0), c(0,0,0,0,1))
t0 <- Sys.time()
mF <- msm(state ~ t, subject=pid, data=d, qmatrix=Q, ematrix=E, deathexact=5, est.initprobs=TRUE, initprobs=c(0.67,0.23,0.08,0.02,0), control=list(fnscale=90000,maxit=6000,reltol=1e-8), method="BFGS")
cat(sprintf("FULL v4, initial state estimated | converged %s | -2LL %.1f | %.1f min\n", mF$opt$convergence==0, mF$minus2loglik, as.numeric(difftime(Sys.time(),t0,units="mins"))))
cat("initial true-state probabilities (est, lo, hi):\n"); print(round(mF$hmodel$initprobs,3))
cat("emission (rows=true, cols=reported blank/mild/mod/sev):\n"); print(round(ematrix.msm(mF)$estimates[1:4,1:4],3)); cat("SE:\n"); print(round(ematrix.msm(mF)$SE[1:4,1:4],3))
cat("intensities /yr:\n"); print(round(qmatrix.msm(mF)$estimates[1:4,],4)); cat("sojourn yrs:\n"); print(sojourn.msm(mF))
ip <- mF$hmodel$initprobs; if (is.matrix(ip)) ip <- ip[,1]; Em <- ematrix.msm(mF)$estimates; pb <- ip[1:4]*Em[1:4,1]
cat(sprintf("P(true state | reported BLANK at index): MAC-free %.3f, mild %.3f, moderate %.3f, severe %.3f\n", pb[1]/sum(pb), pb[2]/sum(pb), pb[3]/sum(pb), pb[4]/sum(pb)))
pm <- ip[1:4]*(Em[1:4,3]+Em[1:4,4]); cat(sprintf("model-implied PPV of a moderate-or-greater read at index (true mod/sev): %.3f\n", (pm[3]+pm[4])/sum(pm)))
save(mF, file="outputs/results/hmm.RData")
