# 01_prepare_frame.R: derive analysis variables from the extracted frame (data/frame.csv) and save data/frame.rds (the primary
# cohort, defined only with information available by the end of the index episode; the 56 patients whose only rheumatic report
# came after index are kept and flagged rheum_post = 1 so that a sensitivity analysis can exclude them).
# Body mass index from reported height and weight (implausible values set missing); E/e' restricted to a plausible range;
# eGFR by CKD-EPI 2021 from the lowest (baseline), median and most recent creatinine of the year before index.
source("R/00_common.R")
d <- read.csv("data/frame.csv", stringsAsFactors=FALSE)
d$bmi <- d$wt/(d$ht/100)^2; d$bmi[d$bmi<12 | d$bmi>70 | !is.finite(d$bmi)] <- NA
d$E_sept[d$E_sept>60 | d$E_sept<2] <- NA
k <- ifelse(d$male==1,0.9,0.7); a <- ifelse(d$male==1,-0.302,-0.241); f <- ifelse(d$male==1,1,1.012)
egfr <- function(cr) 142*pmin(cr/k,1)^a*pmax(cr/k,1)^-1.2*0.9938^d$age0*f
d$egfr_base <- egfr(d$cr_min); d$egfr_med <- egfr(d$cr_median); d$egfr_last <- egfr(d$cr_last)
cutg <- function(x) relevel(cut(x, c(-Inf,30,60,Inf), labels=c("lt30","30to60","ge60")), "ge60")
d$eg <- cutg(d$egfr_base); d$ef_cat <- relevel(cut(d$lvef,c(-Inf,39.99,49.99,Inf),labels=c("lt40","40to49","ge50")),"ge50")
d$age10 <- (d$age0-70)/10; d$era_n <- as.numeric(factor(d$era))-3; d$setting <- factor(d$setting, levels=0:3, labels=c("outpt","ED","ward","ICU"))
z <- function(v) as.numeric(scale(v)); d$zE <- z(d$E_sept); d$zla <- z(d$la); d$zivs <- z(d$ivs); d$zbmi <- z(d$bmi); d$zphos <- z(d$phos_median); d$zca <- z(d$ca_median); d$zalp <- z(log(d$alp_median))
d$age_q <- d$age0 + d$t_first/365.25          # age at the qualifying (first moderate or severe) echocardiogram
saveRDS(d, "data/frame.rds"); dp <- d
cat(sprintf("frame prepared: %d patients in the primary cohort (of whom %d have a rheumatic report only after index; a sensitivity analysis excludes them)\n", nrow(d), sum(d$rheum_post==1)))
cat(sprintf("source population: %d adults with a qualifying study; %d with a single episode; %d with two or more; prosthesis/ring %d; rheumatic at/before index %d\n", d$n_adult_any[1], d$n_single_episode[1], d$n_adult[1], d$n_pros[1], d$n_rheum[1]))
cat(sprintf("index episode: %d (%.1f%%) with more than one study; defining study is the first study in %.1f%% of those; index date later than the episode's first study for %d patients\n", sum(dp$idx_n_stud>1), 100*mean(dp$idx_n_stud>1), 100*mean(dp$idx_def_is_first[dp$idx_n_stud>1]==1), sum(dp$idx_def_is_first==0)))
