# 02_impute.R: multiple imputation (m = 20) of index covariates by chained equations; mineral analytes are not imputed.
suppressPackageStartupMessages({library(mice)})
d <- readRDS("data/frame.rds")
mk <- function(tp, td, tl){ ev <- ifelse(!is.na(tp),1, ifelse(!is.na(td) & (is.na(tl) | td<=tl+30),2,0)); tt <- ifelse(!is.na(tp),tp, ifelse(ev==2, td, tl)); list(ev=ev, tt=pmax(tt,1)/365.25) }
zf <- mk(d$t_first, d$t_death, d$t_last); d$ev_first <- as.integer(zf$ev==1); d$ltime <- log(zf$tt); d$died <- as.integer(!is.na(d$t_death))
# imputation model: index-level variables, with outcome indicators as auxiliaries (standard); impute continuous, cut afterwards
imp_vars <- c("age0","male","era_n","setting","sev1","av","lvef","E_sept","la","ivs","bmi","af","egfr_base","esrd","phos_median","ca_median","alp_median","dm","htn","cad","idx_inpt","n_ep","ev_first","ltime","died")
x <- d[, imp_vars]; x$setting <- factor(x$setting); x$av <- as.numeric(x$av); x$sev1 <- factor(x$sev1)
x$phos_median[!is.na(x$phos_median) & (x$phos_median<0.5 | x$phos_median>15)] <- NA
meth <- make.method(x); meth[c("phos_median","ca_median","alp_median")] <- ""   # mineral analytes are NOT imputed: analysed only in the subcohort with values (protocol)
pred <- make.predictorMatrix(x); pred[, c("phos_median","ca_median","alp_median")] <- 0
cat("missing per variable:\n"); print(round(100*colMeans(is.na(x)),1))
t0 <- Sys.time(); set.seed(20260816)
imp <- mice(x, m=20, maxit=5, method=meth, predictorMatrix=pred, printFlag=FALSE)
cat("mice done in", round(as.numeric(difftime(Sys.time(), t0, units="mins")),1), "min\n")
saveRDS(list(imp=imp, pid=d$pid), "data/mice.rds")
