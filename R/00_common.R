# 00_common.R: helpers shared by the analysis scripts (sourced by them, not run on its own).
if (!grepl("UTF-8", Sys.getlocale("LC_CTYPE"), ignore.case=TRUE)) for (loc in c("en_US.UTF-8","C.UTF-8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break

# Observation of the annulus ends at the last echocardiogram. In the time-to-event constructions that need a death time
# (Cox companion, Fine and Gray, cumulative incidence, time-varying-exposure models) death is a competing event when it
# occurs within GRACE_DAYS after the last study; a later death is treated as censoring at the last study, because MIMIC-IV
# records deaths only up to one year after the last hospital contact. Sensitivity analyses vary the window (0, 90, 365 days).
GRACE_DAYS <- 30
mk <- function(tp, td, tl, grace=GRACE_DAYS){
  ev <- ifelse(!is.na(tp), 1L, ifelse(!is.na(td) & (is.na(tl) | td <= tl + grace), 2L, 0L))
  tt <- ifelse(!is.na(tp), tp, ifelse(ev==2L, td, tl))
  list(ev=ev, tt=pmax(tt, 1)/365.25) }

# Piecewise-exponential time bands (years since index) and the interval splitter used by the person-interval models
BAND_CUTS <- c(0.5, 1, 2, 4, 7)
band_of <- function(t_years) cut(t_years, c(-Inf, BAND_CUTS, Inf), labels=c("a","b","c","d","e","f"), right=FALSE)
split_intervals <- function(iv){
  # iv: data frame with t_prev, t (days from index) and ev; every interval is split at each band cut point that falls
  # strictly inside it, the event (if any) is kept in the last segment, and each segment takes the band of its own start.
  iv$a <- iv$t_prev/365.25; iv$b <- pmax(iv$t, iv$t_prev + 1)/365.25
  for (cp in BAND_CUTS) {
    hit <- iv$a < cp & iv$b > cp
    if (any(hit)) { left <- iv[hit,]; left$b <- cp; left$ev <- 0L; iv$a[hit] <- cp; iv <- rbind(iv, left) } }
  iv$dt <- iv$b - iv$a; iv$band <- band_of(iv$a); iv$t_prev <- iv$a*365.25; iv$t <- iv$b*365.25
  iv[order(iv$pid, iv$a), setdiff(names(iv), c("a","b"))] }
