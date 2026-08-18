-- Analysis frame: one row per patient (adults with >=2 transthoracic echocardiography episodes carrying the annular
-- calcification field). Sources: MIMIC-IV v3.1 hospital and ICU modules, MIMIC-IV-ECHO v1.0 structured measurements,
-- MIMIC-IV-ED, MIMIC-IV-ECG (all `physionet-data` on BigQuery; credentialed access required).
-- Episodes: studies within one hospital stay (overlapping admissions merged), or within one calendar week outside an admission,
-- form one episode; runs whose calendar days overlap are merged. One physical study defines each episode (highest grade, then
-- earliest, then lowest identifier): it supplies the episode's grade and, for later episodes, the episode's date; at index it
-- also supplies the baseline covariates. Time zero is the last study of the index episode, so that everything used as baseline
-- information precedes time zero, and follow-up outcomes are ascertained only on studies after that date. Index = first
-- episode. Exclusions: mitral prosthesis or ring at or before the end of the index episode; rheumatic mitral valve morphology
-- reported by then (patients whose only rheumatic report is later are kept with rheum_post = 1 so that a sensitivity analysis
-- can exclude them). Times are day offsets from time zero; patient identifier is a hash of subject_id.
-- The shared cohort/episode chain is materialised once as temporary tables (a multi-statement script), because inlining it at each
-- of its many references makes BigQuery re-evaluate the admission windowing hundreds of times.
CREATE TEMP TABLE raw AS 
  SELECT subject_id, measurement_id, measurement_datetime dt,
    CASE TRIM(LOWER(result)) WHEN "severe" THEN 3 WHEN "mod mac" THEN 2 WHEN "mild" THEN 1 ELSE 0 END sev
  FROM `physionet-data.mimiciv_echo.structured_measurement` WHERE measurement="mac_severity" AND LOWER(test_type)="tte";
CREATE TEMP TABLE stays AS
WITH adm0 AS (SELECT subject_id, hadm_id, admittime, dischtime, IF(admittime > IFNULL(MAX(dischtime) OVER (PARTITION BY subject_id ORDER BY admittime, hadm_id ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING), DATETIME("1900-01-01")), 1, 0) newstay
         FROM `physionet-data.mimiciv_3_1_hosp.admissions`)
SELECT subject_id, stay, MIN(admittime) s0, MAX(dischtime) s1 FROM (SELECT *, SUM(newstay) OVER (PARTITION BY subject_id ORDER BY admittime, hadm_id) stay FROM adm0) GROUP BY 1,2;
CREATE TEMP TABLE seq AS
WITH adm AS (SELECT r.*, s.stay hadm_id FROM raw r LEFT JOIN stays s ON s.subject_id=r.subject_id AND r.dt BETWEEN s.s0 AND s.s1),
keyed AS (SELECT *, COALESCE(CONCAT("stay_", CAST(hadm_id AS STRING)), CONCAT("out_", CAST(DATE_TRUNC(DATE(dt), WEEK) AS STRING))) epi_key FROM adm),
-- consecutive studies (in time order) with the same key form a run, so an outpatient week cannot span an admission
runs0 AS (SELECT *, IF(epi_key = LAG(epi_key) OVER (PARTITION BY subject_id ORDER BY dt, measurement_id), 0, 1) newrun FROM keyed),
grp AS (SELECT *, SUM(newrun) OVER (PARTITION BY subject_id ORDER BY dt, measurement_id ROWS UNBOUNDED PRECEDING) epi FROM runs0),
ep0 AS (SELECT subject_id, epi, MIN(DATE(dt)) d_first, MAX(DATE(dt)) d_last FROM grp GROUP BY 1,2),
-- runs that overlap in calendar days (for example an emergency department study before an admission and an inpatient study later that
-- day) are one episode: gaps-and-islands over the runs' calendar-day spans
ep1 AS (SELECT *, IF(d_first > IFNULL(MAX(d_last) OVER (PARTITION BY subject_id ORDER BY d_first, epi ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING), DATE "1900-01-01"), 1, 0) newep FROM ep0),
ep2 AS (SELECT subject_id, epi, SUM(newep) OVER (PARTITION BY subject_id ORDER BY d_first, epi ROWS UNBOUNDED PRECEDING) island FROM ep1),
grp2 AS (SELECT g.*, e.island FROM grp g JOIN ep2 e USING(subject_id, epi)),
ep AS (SELECT subject_id, CAST(island AS STRING) epi, ARRAY_AGG(dt ORDER BY sev DESC, dt ASC, measurement_id ASC LIMIT 1)[OFFSET(0)] dt, MIN(dt) dt_first, MAX(dt) dt_last, COUNT(*) n_stud, MAX(sev) sev,
       ARRAY_AGG(measurement_id ORDER BY sev DESC, dt ASC, measurement_id ASC LIMIT 1)[OFFSET(0)] mid, MAX(IF(hadm_id IS NULL,0,1)) inpt FROM grp2 GROUP BY 1,2)
SELECT *, ROW_NUMBER() OVER (PARTITION BY subject_id ORDER BY dt, epi) rn, COUNT(*) OVER (PARTITION BY subject_id) n_ep, LEAD(sev) OVER (PARTITION BY subject_id ORDER BY dt, epi) nxt, LEAD(dt) OVER (PARTITION BY subject_id ORDER BY dt, epi) nxt_dt FROM ep;
CREATE TEMP TABLE all1 AS
WITH pt AS (SELECT subject_id, anchor_age, anchor_year, anchor_year_group, gender, dod FROM `physionet-data.mimiciv_3_1_hosp.patients`)
SELECT s.subject_id, s.mid, s.dt d1, s.dt_first d1_first, s.dt_last d1_end, s.n_stud idx_n_stud, s.sev sev1, s.n_ep, p.anchor_age+(EXTRACT(YEAR FROM s.dt_last)-p.anchor_year) age0, IF(p.gender="M",1,0) male, p.anchor_year_group era, p.dod
       FROM seq s JOIN pt p USING(subject_id) WHERE s.rn=1;
CREATE TEMP TABLE adult AS SELECT * FROM all1 WHERE n_ep>=2 AND age0>=18;
CREATE TEMP TABLE pros AS SELECT DISTINCT s.subject_id FROM `physionet-data.mimiciv_echo.structured_measurement` s JOIN adult a USING(subject_id)
   WHERE LOWER(s.test_type)="tte" AND s.measurement_datetime<=a.d1_end AND ((s.measurement LIKE "mvr_%" AND s.result IS NOT NULL AND TRIM(s.result)!="") OR (s.measurement="mv_leaflets" AND REGEXP_CONTAINS(LOWER(s.result), r"prosth|mechanical|bioprosth|annuloplasty|ring")));
CREATE TEMP TABLE rheum_any AS SELECT s.subject_id, MIN(s.measurement_datetime) first_rheum FROM `physionet-data.mimiciv_echo.structured_measurement` s JOIN adult a USING(subject_id)
   WHERE LOWER(s.test_type)="tte" AND s.measurement IN ("rheumatic_mv","mv_prolapse") AND LOWER(TRIM(s.result))="rheumatic deformity" GROUP BY 1;
CREATE TEMP TABLE coh AS
WITH rheum_pre AS (SELECT r.subject_id FROM rheum_any r JOIN adult a USING(subject_id) WHERE r.first_rheum<=a.d1_end),
rheum_post AS (SELECT r.subject_id FROM rheum_any r JOIN adult a USING(subject_id) WHERE r.first_rheum>a.d1_end)
SELECT a.*, IF(a.subject_id IN (SELECT subject_id FROM rheum_post),1,0) rheum_post FROM adult a WHERE a.subject_id NOT IN (SELECT subject_id FROM pros) AND a.subject_id NOT IN (SELECT subject_id FROM rheum_pre);
CREATE TEMP TABLE counts AS SELECT
  (SELECT COUNT(*) FROM all1) n_any, (SELECT COUNT(*) FROM all1 WHERE age0>=18) n_adult_any, (SELECT COUNT(*) FROM all1 WHERE age0>=18 AND n_ep=1) n_single_episode, (SELECT COUNT(*) FROM all1 WHERE n_ep>=2) n_s0, (SELECT COUNT(*) FROM adult) n_adult,
  (SELECT COUNT(*) FROM adult WHERE subject_id IN (SELECT subject_id FROM pros)) n_pros,
  (SELECT COUNT(*) FROM adult a JOIN rheum_any r USING(subject_id) WHERE a.subject_id NOT IN (SELECT subject_id FROM pros) AND r.first_rheum<=a.d1_end) n_rheum,
  (SELECT COUNT(*) FROM coh WHERE rheum_post=1) n_rheum_post;
WITH later AS (SELECT c.subject_id, s.dt, s.dt_first, s.dt_last, s.sev, s.nxt, s.nxt_dt, s.rn FROM coh c JOIN seq s USING(subject_id) WHERE s.rn>1),
ev AS (SELECT subject_id, MIN(IF(sev>=2,dt,NULL)) t_first, MIN(IF(sev>=2 AND nxt>=2,dt,NULL)) t_conf, MAX(dt) t_last, MAX(dt_last) t_last_study,
       ARRAY_AGG(IF(sev>=2, IF(nxt IS NULL, "unconfirmable", IF(nxt>=2,"confirmed","refuted")), NULL) IGNORE NULLS ORDER BY dt, rn LIMIT 1)[SAFE_OFFSET(0)] first_status,
       ARRAY_AGG(IF(sev>=2, nxt_dt, NULL) IGNORE NULLS ORDER BY dt, rn LIMIT 1)[SAFE_OFFSET(0)] t_first_next,      -- date of the episode after the first moderate or severe episode (classification echo)
       ARRAY_AGG(IF(sev>=2 AND nxt>=2, nxt_dt, NULL) IGNORE NULLS ORDER BY dt, rn LIMIT 1)[SAFE_OFFSET(0)] t_conf_next,   -- the echo that confirms the first confirmed pair (differs from t_first_next when a first read is refuted and a later pair confirms)
       ARRAY_AGG(IF(sev>=2, dt_last, NULL) IGNORE NULLS ORDER BY dt, rn LIMIT 1)[SAFE_OFFSET(0)] t_first_end,     -- last study of the first moderate or severe episode
       ARRAY_AGG(IF(sev>=2, dt_first, NULL) IGNORE NULLS ORDER BY dt, rn LIMIT 1)[SAFE_OFFSET(0)] t_first_start,   -- first study of the first moderate or severe episode
       MIN(IF(sev>=1,dt,NULL)) t_first_mild_or_worse FROM later GROUP BY 1),
-- index-echo measured exposures
echo AS (SELECT s.measurement_id,
   MAX(IF(s.measurement="av_leaflets", CASE WHEN REGEXP_CONTAINS(LOWER(s.result),r"severe thick") THEN 3 WHEN REGEXP_CONTAINS(LOWER(s.result),r"mod thick") THEN 2 WHEN REGEXP_CONTAINS(LOWER(s.result),r"mild thick") THEN 1 WHEN REGEXP_CONTAINS(LOWER(s.result),r"^nl") THEN 0 END, NULL)) av,
   MAX(IF(s.measurement="lvef", SAFE_CAST(s.result AS FLOAT64), NULL)) lvef,
   MAX(IF(s.measurement="mv_peak_e", SAFE_CAST(s.result AS FLOAT64), NULL)) E,
   MAX(IF(s.measurement="sept_e_prime", SAFE_CAST(s.result AS FLOAT64), NULL)) sept_ep,
   MAX(IF(s.measurement="lat_e_prime", SAFE_CAST(s.result AS FLOAT64), NULL)) lat_ep,
   MAX(IF(s.measurement="la_dimen", SAFE_CAST(s.result AS FLOAT64), NULL)) la,
   MAX(IF(s.measurement="septal_thickness", SAFE_CAST(s.result AS FLOAT64), NULL)) ivs,
   COALESCE(MAX(IF(s.measurement="height_cm", SAFE_CAST(s.result AS FLOAT64), NULL)), MAX(IF(s.measurement="height_inches", SAFE_CAST(s.result AS FLOAT64)*2.54, NULL))) ht,
   COALESCE(MAX(IF(s.measurement="weight_kg", SAFE_CAST(s.result AS FLOAT64), NULL)), MAX(IF(s.measurement="weight_pounds", SAFE_CAST(s.result AS FLOAT64)*0.4536, NULL))) wt,
   MAX(IF(s.measurement="mitral_regurg" AND REGEXP_CONTAINS(LOWER(s.result),r"^mod|sever"),1,0)) mr0,
   MAX(IF((s.measurement="mitral_stenosis" AND REGEXP_CONTAINS(LOWER(s.result),r"^mod|sever")) OR (s.measurement="mv_mean_grad" AND SAFE_CAST(s.result AS FLOAT64)>=5),1,0)) ms0,
   MAX(IF(s.measurement="mv_mean_grad", SAFE_CAST(s.result AS FLOAT64), NULL)) grad0
   FROM `physionet-data.mimiciv_echo.structured_measurement` s JOIN coh c ON c.mid=s.measurement_id GROUP BY 1),
-- baseline valve status for the stenosis and regurgitation outcomes is judged on EVERY study of the index episode (not only the defining study)
idxvalve AS (SELECT c.subject_id,
   MAX(IF(s.measurement="mitral_regurg" AND REGEXP_CONTAINS(LOWER(s.result),r"^mod|sever"),1,0)) mr0,
   MAX(IF((s.measurement="mitral_stenosis" AND REGEXP_CONTAINS(LOWER(s.result),r"^mod|sever")) OR (s.measurement="mv_mean_grad" AND SAFE_CAST(s.result AS FLOAT64)>=5),1,0)) ms0
   FROM `physionet-data.mimiciv_echo.structured_measurement` s JOIN coh c ON s.subject_id=c.subject_id AND s.measurement_datetime BETWEEN c.d1_first AND c.d1_end WHERE LOWER(s.test_type)="tte" GROUP BY 1),
-- last hospital discharge: MIMIC-IV records deaths (hospital and state records) up to one year after the last hospital discharge, so vital status is ascertained to that date (or to the last echocardiogram if later)
lastadm AS (SELECT c.subject_id, MAX(a.dischtime) last_disch FROM coh c JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a USING(subject_id) GROUP BY 1),
-- study-level mitral function after index (secondary outcomes)
stud AS (SELECT s.subject_id, s.measurement_id, s.measurement_datetime sdt,
   MAX(IF(s.measurement="mitral_regurg" AND REGEXP_CONTAINS(LOWER(s.result),r"^mod|sever"),1,0)) mr,
   MAX(IF((s.measurement="mitral_stenosis" AND REGEXP_CONTAINS(LOWER(s.result),r"^mod|sever")) OR (s.measurement="mv_mean_grad" AND SAFE_CAST(s.result AS FLOAT64)>=5),1,0)) ms
   FROM `physionet-data.mimiciv_echo.structured_measurement` s JOIN coh c USING(subject_id) WHERE LOWER(s.test_type)="tte" AND s.measurement_datetime>c.d1_end GROUP BY 1,2,3),
dys AS (SELECT subject_id, MIN(IF(ms=1,sdt,NULL)) t_ms, MIN(IF(mr=1,sdt,NULL)) t_mr FROM stud GROUP BY 1),
-- labs in 365d before index
idxadm AS (SELECT c.subject_id, MIN(a.hadm_id) idx_hadm FROM coh c JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a ON a.subject_id=c.subject_id AND c.d1 BETWEEN a.admittime AND a.dischtime GROUP BY 1),
lv0 AS (SELECT c.subject_id, d.label, l.valuenum v, l.charttime ct, IF(l.hadm_id IS NOT NULL AND l.hadm_id=i.idx_hadm,1,0) idxadm
   FROM coh c LEFT JOIN idxadm i USING(subject_id) JOIN `physionet-data.mimiciv_3_1_hosp.labevents` l ON l.subject_id=c.subject_id AND DATETIME(l.charttime) BETWEEN DATETIME_SUB(c.d1_end, INTERVAL 365 DAY) AND c.d1_end
   JOIN `physionet-data.mimiciv_3_1_hosp.d_labitems` d USING(itemid) WHERE d.label IN ("Creatinine","Phosphate","Calcium, Total","Alkaline Phosphatase") AND l.valuenum IS NOT NULL AND l.valuenum>0),
lv AS (SELECT *, PERCENTILE_DISC(v,0.5) OVER (PARTITION BY subject_id,label) med, FIRST_VALUE(v) OVER (PARTITION BY subject_id,label ORDER BY ct DESC, v DESC) lastv,
   FIRST_VALUE(v) OVER (PARTITION BY subject_id,label,idxadm ORDER BY ct DESC, v DESC) lastv_pre FROM lv0),
lab AS (SELECT subject_id,
   COUNTIF(label="Creatinine") n_cr, MIN(IF(label="Creatinine",v,NULL)) cr_min, MAX(IF(label="Creatinine",med,NULL)) cr_median, MAX(IF(label="Creatinine",lastv,NULL)) cr_last,
   COUNTIF(label="Phosphate") n_phos, MAX(IF(label="Phosphate",med,NULL)) phos_median, MAX(IF(label="Phosphate",lastv,NULL)) phos_last,
   MAX(IF(label="Phosphate" AND idxadm=0,lastv_pre,NULL)) phos_preidx_last,
   MAX(IF(label="Calcium, Total",med,NULL)) ca_median, MAX(IF(label="Calcium, Total",lastv,NULL)) ca_last,
   MAX(IF(label="Alkaline Phosphatase",med,NULL)) alp_median, MAX(IF(label="Alkaline Phosphatase",lastv,NULL)) alp_last
   FROM lv GROUP BY 1),
-- diagnoses at or before index, ICD-version restricted; code families verified against d_icd_diagnoses
dxdt AS (SELECT d.subject_id, d.icd_code, d.icd_version, a.admittime, a.dischtime FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd` d JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a USING(hadm_id) JOIN coh c2 ON c2.subject_id=d.subject_id),
cm AS (SELECT c.subject_id,
   MAX(IF((x.icd_version=10 AND REGEXP_CONTAINS(x.icd_code,r"^(Z992|N186)")) OR (x.icd_version=9 AND REGEXP_CONTAINS(x.icd_code,r"^(V4511|V4512|5856)")),1,0)) esrd_dx,
   MAX(IF((x.icd_version=10 AND REGEXP_CONTAINS(x.icd_code,r"^(E08|E09|E1[0-4])")) OR (x.icd_version=9 AND REGEXP_CONTAINS(x.icd_code,r"^(249|250)")),1,0)) dm,
   MAX(IF((x.icd_version=10 AND REGEXP_CONTAINS(x.icd_code,r"^I1[0-6]")) OR (x.icd_version=9 AND REGEXP_CONTAINS(x.icd_code,r"^40[1-5]")),1,0)) htn,
   MAX(IF((x.icd_version=10 AND REGEXP_CONTAINS(x.icd_code,r"^(I2[0-5]|Z951|Z955)")) OR (x.icd_version=9 AND REGEXP_CONTAINS(x.icd_code,r"^(41[0-4]|V4581|V4582)")),1,0)) cad,
   MAX(IF((x.icd_version=10 AND REGEXP_CONTAINS(x.icd_code,r"^I48")) OR (x.icd_version=9 AND REGEXP_CONTAINS(x.icd_code,r"^4273")),1,0)) af_icd
   FROM coh c LEFT JOIN dxdt x ON x.subject_id=c.subject_id AND x.admittime<=c.d1_end GROUP BY 1),
-- MIMIC-IV diagnosis codes are assigned per admission at discharge and carry no date within the stay, so a code from an admission
-- that was still open at time zero may describe a condition documented later in that stay: the main definition (cm) takes codes
-- from admissions that BEGAN by the end of the index episode; the sensitivity definition (cmp) takes only admissions COMPLETED by then
cmp AS (SELECT c.subject_id,
   MAX(IF((x.icd_version=10 AND REGEXP_CONTAINS(x.icd_code,r"^(Z992|N186)")) OR (x.icd_version=9 AND REGEXP_CONTAINS(x.icd_code,r"^(V4511|V4512|5856)")),1,0)) esrd_dx_prior,
   MAX(IF((x.icd_version=10 AND REGEXP_CONTAINS(x.icd_code,r"^(E08|E09|E1[0-4])")) OR (x.icd_version=9 AND REGEXP_CONTAINS(x.icd_code,r"^(249|250)")),1,0)) dm_prior,
   MAX(IF((x.icd_version=10 AND REGEXP_CONTAINS(x.icd_code,r"^I1[0-6]")) OR (x.icd_version=9 AND REGEXP_CONTAINS(x.icd_code,r"^40[1-5]")),1,0)) htn_prior,
   MAX(IF((x.icd_version=10 AND REGEXP_CONTAINS(x.icd_code,r"^(I2[0-5]|Z951|Z955)")) OR (x.icd_version=9 AND REGEXP_CONTAINS(x.icd_code,r"^(41[0-4]|V4581|V4582)")),1,0)) cad_prior,
   MAX(IF((x.icd_version=10 AND REGEXP_CONTAINS(x.icd_code,r"^I48")) OR (x.icd_version=9 AND REGEXP_CONTAINS(x.icd_code,r"^4273")),1,0)) af_icd_prior
   FROM coh c LEFT JOIN dxdt x ON x.subject_id=c.subject_id AND x.dischtime<=c.d1_end GROUP BY 1),
-- dialysis from ICU derived table
ihd AS (SELECT DISTINCT c.subject_id FROM coh c JOIN `physionet-data.mimiciv_3_1_icu.icustays` ie USING(subject_id) JOIN `physionet-data.mimiciv_3_1_derived.rrt` r USING(stay_id) WHERE r.dialysis_active=1 AND r.dialysis_type IN ("IHD","Peritoneal") AND DATETIME(r.charttime)<=c.d1_end),
crrt AS (SELECT DISTINCT c.subject_id FROM coh c JOIN `physionet-data.mimiciv_3_1_icu.icustays` ie USING(subject_id) JOIN `physionet-data.mimiciv_3_1_derived.rrt` r USING(stay_id) WHERE r.dialysis_type IN ("CRRT","CVVHDF","CVVHD","CVVH","SCUF") AND DATETIME(r.charttime)<=c.d1_end),
-- AF on ECG machine report between 30 days before index and the end of the index episode
ecgaf AS (SELECT c.subject_id, MAX(1) ecg30, MAX(IF(REGEXP_CONTAINS(LOWER(CONCAT(IFNULL(m.report_0,''),' ',IFNULL(m.report_1,''),' ',IFNULL(m.report_2,''),' ',IFNULL(m.report_3,''),' ',IFNULL(m.report_4,''),' ',IFNULL(m.report_5,''))), r"atrial fib|a-?fib|atrial flutter"),1,0)) af_ecg
   FROM coh c JOIN `physionet-data.mimiciv_ecg.machine_measurements` m ON m.subject_id=c.subject_id AND DATETIME(m.ecg_time) BETWEEN DATETIME_SUB(c.d1, INTERVAL 30 DAY) AND c.d1_end GROUP BY 1),   -- ECG rhythm from 30 days before index to the end of the index episode
prior AS (SELECT c.subject_id, LEAST(IFNULL(a.fa, c.d1), IFNULL(e.fe, c.d1), IFNULL(r.fr, c.d1)) fa FROM coh c
   LEFT JOIN (SELECT subject_id, MIN(admittime) fa FROM `physionet-data.mimiciv_3_1_hosp.admissions` GROUP BY 1) a USING(subject_id)
   LEFT JOIN (SELECT subject_id, MIN(intime) fe FROM `physionet-data.mimiciv_ed.edstays` GROUP BY 1) e USING(subject_id)
   LEFT JOIN (SELECT subject_id, MIN(dt) fr FROM raw GROUP BY 1) r USING(subject_id)),   -- first contact of any kind: admission, ED stay or echocardiogram
setting AS (SELECT c.subject_id, CASE WHEN icu.subject_id IS NOT NULL THEN 3 WHEN ad.subject_id IS NOT NULL THEN 2 WHEN ed.subject_id IS NOT NULL THEN 1 ELSE 0 END st FROM coh c
   LEFT JOIN (SELECT DISTINCT c.subject_id FROM coh c JOIN `physionet-data.mimiciv_3_1_icu.icustays` i ON i.subject_id=c.subject_id AND c.d1 BETWEEN i.intime AND i.outtime) icu USING(subject_id)
   LEFT JOIN (SELECT DISTINCT c.subject_id FROM coh c JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a ON a.subject_id=c.subject_id AND c.d1 BETWEEN a.admittime AND a.dischtime) ad USING(subject_id)
   LEFT JOIN (SELECT DISTINCT c.subject_id FROM coh c JOIN `physionet-data.mimiciv_ed.edstays` e ON e.subject_id=c.subject_id AND c.d1 BETWEEN e.intime AND e.outtime) ed USING(subject_id)),
-- mitral intervention: ICD-9-CM / ICD-10-PCS procedure codes verified against d_icd_procedures
prc AS (SELECT p.subject_id, p.hadm_id, p.chartdate, p.icd_code, p.icd_version FROM `physionet-data.mimiciv_3_1_hosp.procedures_icd` p JOIN coh USING(subject_id)),
ms_dx_adm AS (SELECT DISTINCT d.hadm_id FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd` d WHERE (d.icd_version=10 AND REGEXP_CONTAINS(d.icd_code,r"^(I050|I052|I342)")) OR (d.icd_version=9 AND REGEXP_CONTAINS(d.icd_code,r"^(3940|3942)"))),
mv AS (SELECT subject_id, hadm_id, MIN(chartdate) mvdate,
   MAX(CASE WHEN icd_version=10 AND REGEXP_CONTAINS(icd_code,r"^02RG3") THEN "tmvr" WHEN icd_version=10 AND REGEXP_CONTAINS(icd_code,r"^02RG") THEN "smvr" WHEN icd_version=9 AND icd_code IN ("3523","3524") THEN "smvr"
            WHEN (icd_version=10 AND REGEXP_CONTAINS(icd_code,r"^02UG3")) OR (icd_version=9 AND icd_code="3597") THEN "teer"
            WHEN (icd_version=10 AND REGEXP_CONTAINS(icd_code,r"^027G")) OR (icd_version=9 AND icd_code="3596") THEN "balloon"
            ELSE "repair" END) modality
   FROM prc WHERE (icd_version=10 AND REGEXP_CONTAINS(icd_code, r"^02[7BCNQRUVW]G")) OR (icd_version=9 AND (icd_code IN ("3502","3512","3523","3524","3532","3533","3597") OR (icd_code="3596" AND hadm_id IN (SELECT hadm_id FROM ms_dx_adm)))) GROUP BY 1,2),
conc AS (SELECT m.subject_id, m.hadm_id, MAX(IF((p.icd_version=10 AND REGEXP_CONTAINS(p.icd_code,r"^021[0-3]")) OR (p.icd_version=9 AND p.icd_code LIKE "361%"),1,0)) cabg,
   MAX(IF((p.icd_version=10 AND REGEXP_CONTAINS(p.icd_code,r"^02[QRU]F")) OR (p.icd_version=9 AND p.icd_code IN ("3511","3521","3522","3505","3506")),1,0)) avs FROM mv m JOIN prc p USING(subject_id, hadm_id) GROUP BY 1,2),
firstmv AS (SELECT m.subject_id, MIN(m.mvdate) mvdate, ARRAY_AGG(m.modality ORDER BY m.mvdate, m.hadm_id LIMIT 1)[OFFSET(0)] modality, ARRAY_AGG(IF(c.cabg=1 OR c.avs=1,0,1) ORDER BY m.mvdate, m.hadm_id LIMIT 1)[OFFSET(0)] isolated FROM mv m JOIN conc c USING(subject_id, hadm_id) GROUP BY 1),
-- moderate or greater mitral dysfunction: first study, from the first study of the qualifying (first moderate or severe) episode onward, with moderate or greater MR or MS, or mean gradient >= 5 mmHg;
-- dysfunction on any study of that episode counts as already present at the qualifying episode
modany AS (SELECT c.subject_id, IF(c.sev1>=2, c.d1, e.t_first) t_mod, IF(c.sev1>=2, c.d1_end, e.t_first_end) t_mod_end, IF(c.sev1>=2, c.d1_first, e.t_first_start) t_mod_start,
                  IF(c.sev1>=2, c.d1, IF(e.t_conf IS NOT NULL, e.t_conf_next, NULL)) t_conf_start   -- confirmed definition: index for moderate or severe at index, otherwise the echo that confirms the pair
           FROM coh c LEFT JOIN ev e USING(subject_id)),
stud2 AS (SELECT s.subject_id, s.measurement_datetime sdt,
   MAX(IF(s.measurement="mitral_regurg" AND REGEXP_CONTAINS(LOWER(s.result),r"^mod|sever"),1,0)) mr,
   MAX(IF((s.measurement="mitral_stenosis" AND REGEXP_CONTAINS(LOWER(s.result),r"^mod|sever")) OR (s.measurement="mv_mean_grad" AND SAFE_CAST(s.result AS FLOAT64)>=5),1,0)) ms
   FROM `physionet-data.mimiciv_echo.structured_measurement` s JOIN modany m USING(subject_id) WHERE LOWER(s.test_type)="tte" AND m.t_mod IS NOT NULL AND s.measurement_datetime>=m.t_mod_start GROUP BY 1,2),   -- from the first study of the qualifying episode
pheno AS (SELECT s.subject_id, MIN(IF(s.mr=1 OR s.ms=1, s.sdt, NULL)) t_pheno, ARRAY_AGG(IF(s.mr=1 OR s.ms=1, IF(s.ms=1,"sten","mr"), NULL) IGNORE NULLS ORDER BY s.sdt, s.ms DESC LIMIT 1)[SAFE_OFFSET(0)] pheno_type,
          MAX(IF((s.mr=1 OR s.ms=1) AND s.sdt<=m.t_mod_end,1,0)) pheno_same_episode,   -- dysfunction already documented within the qualifying episode
          MIN(IF((s.mr=1 OR s.ms=1) AND s.sdt>m.t_mod_end, s.sdt, NULL)) t_pheno_after,   -- first dysfunction on a study after that episode
          MIN(IF((s.mr=1 OR s.ms=1) AND m.t_conf_start IS NOT NULL AND s.sdt>=m.t_conf_start, s.sdt, NULL)) t_pheno_conf,   -- confirmed definition: first dysfunction at or after the confirming episode
          ARRAY_AGG(IF((s.mr=1 OR s.ms=1) AND m.t_conf_start IS NOT NULL AND s.sdt>=m.t_conf_start, IF(s.ms=1,"sten","mr"), NULL) IGNORE NULLS ORDER BY s.sdt, s.ms DESC LIMIT 1)[SAFE_OFFSET(0)] pheno_type_conf
   FROM stud2 s JOIN modany m USING(subject_id) GROUP BY 1)
SELECT TO_HEX(MD5(CAST(c.subject_id AS STRING))) pid, c.sev1, c.n_ep, c.age0, c.male, c.era,
  DATE_DIFF(DATE(e.t_first),DATE(c.d1_end),DAY) t_first, DATE_DIFF(DATE(e.t_conf),DATE(c.d1_end),DAY) t_conf, DATE_DIFF(DATE(e.t_last),DATE(c.d1_end),DAY) t_last, DATE_DIFF(DATE(e.t_last_study),DATE(c.d1_end),DAY) t_last_study, e.first_status,
  DATE_DIFF(DATE(e.t_first_next),DATE(c.d1_end),DAY) t_first_next, DATE_DIFF(DATE(e.t_conf_next),DATE(c.d1_end),DAY) t_conf_next, DATE_DIFF(DATE(e.t_first_end),DATE(c.d1_end),DAY) t_first_end,
  c.idx_n_stud, DATE_DIFF(DATE(c.d1_end),DATE(c.d1_first),DAY) idx_span_days, IF(c.d1=c.d1_first,1,0) idx_def_is_first, DATE_DIFF(DATE(c.d1_end),DATE(c.d1),DAY) idx_def_before_end, c.rheum_post,
  DATE_DIFF(DATE(e.t_first_mild_or_worse),DATE(c.d1_end),DAY) t_first_mild,
  IF(c.dod IS NOT NULL, DATE_DIFF(DATE(c.dod),DATE(c.d1_end),DAY), NULL) t_death,
  IF(DATE_DIFF(DATE(c.d1_end),DATE(p.fa),DAY)>=365,1,0) yr1, st.st setting, IF(ia.idx_hadm IS NULL,0,1) idx_inpt,
  ec.av, ec.lvef, SAFE_DIVIDE(ec.E, ec.sept_ep) E_sept, SAFE_DIVIDE(ec.E, ec.lat_ep) E_lat, ec.la, ec.ivs, ec.ht, ec.wt, IFNULL(ivx.mr0, ec.mr0) mr0, IFNULL(ivx.ms0, ec.ms0) ms0, ec.grad0, DATE_DIFF(DATE(la.last_disch), DATE(c.d1_end), DAY) t_last_disch,
  DATE_DIFF(DATE(dy.t_ms),DATE(c.d1_end),DAY) t_ms, DATE_DIFF(DATE(dy.t_mr),DATE(c.d1_end),DAY) t_mr,
  l.n_cr, l.cr_min, l.cr_median, l.cr_last, l.n_phos, l.phos_median, l.phos_last, l.phos_preidx_last, l.ca_median, l.ca_last, l.alp_median, l.alp_last,
  IF(cm.esrd_dx=1 OR c.subject_id IN (SELECT subject_id FROM ihd),1,0) esrd,
  IF(c.subject_id IN (SELECT subject_id FROM crrt) AND cm.esrd_dx=0 AND c.subject_id NOT IN (SELECT subject_id FROM ihd),1,0) aki_rrt,
  cm.dm, cm.htn, cm.cad, cm.af_icd, IFNULL(ea.ecg30,0) ecg30, IFNULL(ea.af_ecg,0) af_ecg, IF(cm.af_icd=1 OR IFNULL(ea.af_ecg,0)=1,1,0) af,
  IFNULL(cmp.dm_prior,0) dm_prior, IFNULL(cmp.htn_prior,0) htn_prior, IFNULL(cmp.cad_prior,0) cad_prior, IF(IFNULL(cmp.af_icd_prior,0)=1 OR IFNULL(ea.af_ecg,0)=1,1,0) af_prior,
  IF(IFNULL(cmp.esrd_dx_prior,0)=1 OR c.subject_id IN (SELECT subject_id FROM ihd),1,0) esrd_prior,   /* sensitivity: codes only from admissions completed by time zero (dialysis records and ECGs carry exact times and are unchanged) */
  GREATEST(DATE_DIFF(DATE(ma.t_mod),DATE(c.d1_end),DAY),0) t_modany, GREATEST(DATE_DIFF(DATE(ph.t_pheno),DATE(c.d1_end),DAY),0) t_pheno, ph.pheno_type, ph.pheno_same_episode, DATE_DIFF(DATE(ph.t_pheno_after),DATE(c.d1_end),DAY) t_pheno_after, GREATEST(DATE_DIFF(DATE(ph.t_pheno_conf),DATE(c.d1_end),DAY),0) t_pheno_conf, ph.pheno_type_conf, GREATEST(DATE_DIFF(DATE(ma.t_conf_start),DATE(c.d1_end),DAY),0) t_conf_start,
  DATE_DIFF(fm.mvdate,DATE(c.d1_end),DAY) t_interv, fm.modality interv_modality, fm.isolated interv_isolated,
  k.n_any, k.n_adult_any, k.n_single_episode, k.n_s0, k.n_adult, k.n_pros, k.n_rheum, k.n_rheum_post
FROM coh c CROSS JOIN counts k LEFT JOIN ev e USING(subject_id) LEFT JOIN echo ec ON ec.measurement_id=c.mid LEFT JOIN dys dy USING(subject_id) LEFT JOIN lab l USING(subject_id) LEFT JOIN cm USING(subject_id) LEFT JOIN cmp USING(subject_id) LEFT JOIN idxvalve ivx USING(subject_id) LEFT JOIN lastadm la USING(subject_id)
 LEFT JOIN ecgaf ea USING(subject_id) LEFT JOIN prior p USING(subject_id) LEFT JOIN setting st USING(subject_id) LEFT JOIN idxadm ia USING(subject_id) LEFT JOIN modany ma USING(subject_id) LEFT JOIN pheno ph USING(subject_id) LEFT JOIN firstmv fm USING(subject_id)
