-- Analysis frame: one row per patient (adults with >=2 transthoracic echocardiography episodes carrying the annular
-- calcification field). Sources: MIMIC-IV v3.1 hospital and ICU modules, MIMIC-IV-ECHO v1.0 structured measurements,
-- MIMIC-IV-ED, MIMIC-IV-ECG (all `physionet-data` on BigQuery; credentialed access required).
-- Episodes: studies within one admission, or within one calendar week outside an admission, form one episode carrying
-- the highest grade. Index = first episode. Exclusions: mitral prosthesis or ring at or before index; rheumatic mitral
-- valve on any report. Times are day offsets from the index episode; patient identifier is a hash of subject_id.
WITH raw AS (
  SELECT subject_id, measurement_id, measurement_datetime dt,
    CASE TRIM(LOWER(result)) WHEN "severe" THEN 3 WHEN "mod mac" THEN 2 WHEN "mild" THEN 1 ELSE 0 END sev
  FROM `physionet-data.mimiciv_echo.structured_measurement` WHERE measurement="mac_severity" AND LOWER(test_type)="tte"),
adm AS (SELECT r.*, a.hadm_id FROM raw r LEFT JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a ON a.subject_id=r.subject_id AND r.dt BETWEEN a.admittime AND a.dischtime),
grp AS (SELECT *, COALESCE(CAST(hadm_id AS STRING), CONCAT("out_", CAST(DATE_TRUNC(DATE(dt), WEEK) AS STRING))) epi FROM adm),
ep AS (SELECT subject_id, epi, MIN(dt) dt, MAX(sev) sev, ARRAY_AGG(measurement_id ORDER BY sev DESC, dt ASC, measurement_id ASC LIMIT 1)[OFFSET(0)] mid FROM grp GROUP BY 1,2),
seq AS (SELECT *, ROW_NUMBER() OVER (PARTITION BY subject_id ORDER BY dt, epi) rn, COUNT(*) OVER (PARTITION BY subject_id) n_ep, LEAD(sev) OVER (PARTITION BY subject_id ORDER BY dt, epi) nxt FROM ep),
pt AS (SELECT subject_id, anchor_age, anchor_year, anchor_year_group, gender, dod FROM `physionet-data.mimiciv_3_1_hosp.patients`),
s0 AS (SELECT s.subject_id, s.mid, s.dt d1, s.sev sev1, s.n_ep, p.anchor_age+(EXTRACT(YEAR FROM s.dt)-p.anchor_year) age0, IF(p.gender="M",1,0) male, p.anchor_year_group era, p.dod
       FROM seq s JOIN pt p USING(subject_id) WHERE s.rn=1 AND s.n_ep>=2),
adult AS (SELECT * FROM s0 WHERE age0>=18),
-- exclusions: prosthesis/ring at or before index; rheumatic morphology on any TTE
pros AS (SELECT DISTINCT s.subject_id FROM `physionet-data.mimiciv_echo.structured_measurement` s JOIN adult a USING(subject_id)
   WHERE LOWER(s.test_type)="tte" AND s.measurement_datetime<=a.d1 AND ((s.measurement LIKE "mvr_%" AND s.result IS NOT NULL AND TRIM(s.result)!="") OR (s.measurement="mv_leaflets" AND REGEXP_CONTAINS(LOWER(s.result), r"prosth|mechanical|bioprosth|annuloplasty|ring")))),
rheum_echo AS (SELECT DISTINCT s.subject_id FROM `physionet-data.mimiciv_echo.structured_measurement` s JOIN adult a USING(subject_id)
   WHERE LOWER(s.test_type)="tte" AND s.measurement IN ("rheumatic_mv","mv_prolapse") AND LOWER(TRIM(s.result))="rheumatic deformity"),
coh AS (SELECT * FROM adult WHERE subject_id NOT IN (SELECT subject_id FROM pros) AND subject_id NOT IN (SELECT subject_id FROM rheum_echo)),
-- outcomes from later episodes
later AS (SELECT c.subject_id, s.dt, s.sev, s.nxt, s.rn FROM coh c JOIN seq s USING(subject_id) WHERE s.rn>1),
ev AS (SELECT subject_id, MIN(IF(sev>=2,dt,NULL)) t_first, MIN(IF(sev>=2 AND nxt>=2,dt,NULL)) t_conf, MAX(dt) t_last,
       ARRAY_AGG(IF(sev>=2, IF(nxt IS NULL, "unconfirmable", IF(nxt>=2,"confirmed","refuted")), NULL) IGNORE NULLS ORDER BY dt, rn LIMIT 1)[SAFE_OFFSET(0)] first_status,
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
-- study-level mitral function after index (secondary outcomes)
stud AS (SELECT s.subject_id, s.measurement_id, s.measurement_datetime sdt,
   MAX(IF(s.measurement="mitral_regurg" AND REGEXP_CONTAINS(LOWER(s.result),r"^mod|sever"),1,0)) mr,
   MAX(IF((s.measurement="mitral_stenosis" AND REGEXP_CONTAINS(LOWER(s.result),r"^mod|sever")) OR (s.measurement="mv_mean_grad" AND SAFE_CAST(s.result AS FLOAT64)>=5),1,0)) ms
   FROM `physionet-data.mimiciv_echo.structured_measurement` s JOIN coh c USING(subject_id) WHERE LOWER(s.test_type)="tte" AND s.measurement_datetime>c.d1 GROUP BY 1,2,3),
dys AS (SELECT subject_id, MIN(IF(ms=1,sdt,NULL)) t_ms, MIN(IF(mr=1,sdt,NULL)) t_mr FROM stud GROUP BY 1),
-- labs in 365d before index
idxadm AS (SELECT c.subject_id, MIN(a.hadm_id) idx_hadm FROM coh c JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a ON a.subject_id=c.subject_id AND c.d1 BETWEEN a.admittime AND a.dischtime GROUP BY 1),
lv0 AS (SELECT c.subject_id, d.label, l.valuenum v, l.charttime ct, IF(l.hadm_id IS NOT NULL AND l.hadm_id=i.idx_hadm,1,0) idxadm
   FROM coh c LEFT JOIN idxadm i USING(subject_id) JOIN `physionet-data.mimiciv_3_1_hosp.labevents` l ON l.subject_id=c.subject_id AND DATETIME(l.charttime) BETWEEN DATETIME_SUB(c.d1, INTERVAL 365 DAY) AND c.d1
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
dxdt AS (SELECT d.subject_id, d.icd_code, d.icd_version, a.admittime FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd` d JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a USING(hadm_id) JOIN coh c2 ON c2.subject_id=d.subject_id),
cm AS (SELECT c.subject_id,
   MAX(IF((x.icd_version=10 AND REGEXP_CONTAINS(x.icd_code,r"^(Z992|N186)")) OR (x.icd_version=9 AND REGEXP_CONTAINS(x.icd_code,r"^(V4511|V4512|5856)")),1,0)) esrd_dx,
   MAX(IF((x.icd_version=10 AND REGEXP_CONTAINS(x.icd_code,r"^(E08|E09|E1[0-4])")) OR (x.icd_version=9 AND REGEXP_CONTAINS(x.icd_code,r"^(249|250)")),1,0)) dm,
   MAX(IF((x.icd_version=10 AND REGEXP_CONTAINS(x.icd_code,r"^I1[0-6]")) OR (x.icd_version=9 AND REGEXP_CONTAINS(x.icd_code,r"^40[1-5]")),1,0)) htn,
   MAX(IF((x.icd_version=10 AND REGEXP_CONTAINS(x.icd_code,r"^(I2[0-5]|Z951|Z955)")) OR (x.icd_version=9 AND REGEXP_CONTAINS(x.icd_code,r"^(41[0-4]|V4581|V4582)")),1,0)) cad,
   MAX(IF((x.icd_version=10 AND REGEXP_CONTAINS(x.icd_code,r"^I48")) OR (x.icd_version=9 AND REGEXP_CONTAINS(x.icd_code,r"^4273")),1,0)) af_icd
   FROM coh c LEFT JOIN dxdt x ON x.subject_id=c.subject_id AND x.admittime<=c.d1 GROUP BY 1),
-- dialysis from ICU derived table
ihd AS (SELECT DISTINCT c.subject_id FROM coh c JOIN `physionet-data.mimiciv_3_1_icu.icustays` ie USING(subject_id) JOIN `physionet-data.mimiciv_3_1_derived.rrt` r USING(stay_id) WHERE r.dialysis_active=1 AND r.dialysis_type IN ("IHD","Peritoneal") AND DATETIME(r.charttime)<=c.d1),
crrt AS (SELECT DISTINCT c.subject_id FROM coh c JOIN `physionet-data.mimiciv_3_1_icu.icustays` ie USING(subject_id) JOIN `physionet-data.mimiciv_3_1_derived.rrt` r USING(stay_id) WHERE r.dialysis_type IN ("CRRT","CVVHDF","CVVHD","CVVH","SCUF") AND DATETIME(r.charttime)<=c.d1),
-- AF on ECG machine report within 30 days of index
ecgaf AS (SELECT c.subject_id, MAX(1) ecg30, MAX(IF(REGEXP_CONTAINS(LOWER(CONCAT(IFNULL(m.report_0,''),' ',IFNULL(m.report_1,''),' ',IFNULL(m.report_2,''),' ',IFNULL(m.report_3,''),' ',IFNULL(m.report_4,''),' ',IFNULL(m.report_5,''))), r"atrial fib|a-?fib|atrial flutter"),1,0)) af_ecg
   FROM coh c JOIN `physionet-data.mimiciv_ecg.machine_measurements` m ON m.subject_id=c.subject_id AND ABS(DATETIME_DIFF(DATETIME(m.ecg_time), c.d1, DAY))<=30 GROUP BY 1),
prior AS (SELECT c.subject_id, MIN(a.admittime) fa FROM coh c JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a USING(subject_id) GROUP BY 1),
setting AS (SELECT c.subject_id,
   CASE WHEN EXISTS(SELECT 1 FROM `physionet-data.mimiciv_3_1_icu.icustays` i WHERE i.subject_id=c.subject_id AND c.d1 BETWEEN i.intime AND i.outtime) THEN 3
        WHEN EXISTS(SELECT 1 FROM `physionet-data.mimiciv_3_1_hosp.admissions` a WHERE a.subject_id=c.subject_id AND c.d1 BETWEEN a.admittime AND a.dischtime) THEN 2
        WHEN EXISTS(SELECT 1 FROM `physionet-data.mimiciv_ed.edstays` e WHERE e.subject_id=c.subject_id AND c.d1 BETWEEN e.intime AND e.outtime) THEN 1 ELSE 0 END st FROM coh c),
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
-- moderate or greater mitral dysfunction: first study at or after the first moderate or greater MAC read with moderate or greater MR or MS, or mean gradient >= 5 mmHg
modany AS (SELECT c.subject_id, IF(c.sev1>=2, c.d1, e.t_first) t_mod FROM coh c LEFT JOIN ev e USING(subject_id)),
stud2 AS (SELECT s.subject_id, s.measurement_datetime sdt,
   MAX(IF(s.measurement="mitral_regurg" AND REGEXP_CONTAINS(LOWER(s.result),r"^mod|sever"),1,0)) mr,
   MAX(IF((s.measurement="mitral_stenosis" AND REGEXP_CONTAINS(LOWER(s.result),r"^mod|sever")) OR (s.measurement="mv_mean_grad" AND SAFE_CAST(s.result AS FLOAT64)>=5),1,0)) ms
   FROM `physionet-data.mimiciv_echo.structured_measurement` s JOIN modany m USING(subject_id) WHERE LOWER(s.test_type)="tte" AND m.t_mod IS NOT NULL AND s.measurement_datetime>=m.t_mod GROUP BY 1,2),
pheno AS (SELECT subject_id, MIN(IF(mr=1 OR ms=1, sdt, NULL)) t_pheno, ARRAY_AGG(IF(mr=1 OR ms=1, IF(ms=1,"sten","mr"), NULL) IGNORE NULLS ORDER BY sdt, ms DESC LIMIT 1)[SAFE_OFFSET(0)] pheno_type FROM stud2 GROUP BY 1)
SELECT TO_HEX(MD5(CAST(c.subject_id AS STRING))) pid, c.sev1, c.n_ep, c.age0, c.male, c.era,
  DATE_DIFF(DATE(e.t_first),DATE(c.d1),DAY) t_first, DATE_DIFF(DATE(e.t_conf),DATE(c.d1),DAY) t_conf, DATE_DIFF(DATE(e.t_last),DATE(c.d1),DAY) t_last, e.first_status,
  DATE_DIFF(DATE(e.t_first_mild_or_worse),DATE(c.d1),DAY) t_first_mild,
  IF(c.dod IS NOT NULL, DATE_DIFF(DATE(c.dod),DATE(c.d1),DAY), NULL) t_death,
  IF(DATE_DIFF(DATE(c.d1),DATE(p.fa),DAY)>=365,1,0) yr1, st.st setting, IF(ia.idx_hadm IS NULL,0,1) idx_inpt,
  ec.av, ec.lvef, SAFE_DIVIDE(ec.E, ec.sept_ep) E_sept, SAFE_DIVIDE(ec.E, ec.lat_ep) E_lat, ec.la, ec.ivs, ec.ht, ec.wt, ec.mr0, ec.ms0, ec.grad0,
  DATE_DIFF(DATE(dy.t_ms),DATE(c.d1),DAY) t_ms, DATE_DIFF(DATE(dy.t_mr),DATE(c.d1),DAY) t_mr,
  l.n_cr, l.cr_min, l.cr_median, l.cr_last, l.n_phos, l.phos_median, l.phos_last, l.phos_preidx_last, l.ca_median, l.ca_last, l.alp_median, l.alp_last,
  IF(cm.esrd_dx=1 OR c.subject_id IN (SELECT subject_id FROM ihd),1,0) esrd,
  IF(c.subject_id IN (SELECT subject_id FROM crrt) AND cm.esrd_dx=0 AND c.subject_id NOT IN (SELECT subject_id FROM ihd),1,0) aki_rrt,
  cm.dm, cm.htn, cm.cad, cm.af_icd, IFNULL(ea.ecg30,0) ecg30, IFNULL(ea.af_ecg,0) af_ecg, IF(cm.af_icd=1 OR IFNULL(ea.af_ecg,0)=1,1,0) af,
  DATE_DIFF(DATE(ma.t_mod),DATE(c.d1),DAY) t_modany, DATE_DIFF(DATE(ph.t_pheno),DATE(c.d1),DAY) t_pheno, ph.pheno_type,
  DATE_DIFF(fm.mvdate,DATE(c.d1),DAY) t_interv, fm.modality interv_modality, fm.isolated interv_isolated,
  (SELECT COUNT(*) FROM s0) n_s0, (SELECT COUNT(*) FROM adult) n_adult, (SELECT COUNT(*) FROM adult WHERE subject_id IN (SELECT subject_id FROM pros)) n_pros, (SELECT COUNT(*) FROM adult WHERE subject_id NOT IN (SELECT subject_id FROM pros) AND subject_id IN (SELECT subject_id FROM rheum_echo)) n_rheum
FROM coh c LEFT JOIN ev e USING(subject_id) LEFT JOIN echo ec ON ec.measurement_id=c.mid LEFT JOIN dys dy USING(subject_id) LEFT JOIN lab l USING(subject_id) LEFT JOIN cm USING(subject_id)
 LEFT JOIN ecgaf ea USING(subject_id) LEFT JOIN prior p USING(subject_id) LEFT JOIN setting st USING(subject_id) LEFT JOIN idxadm ia USING(subject_id) LEFT JOIN modany ma USING(subject_id) LEFT JOIN pheno ph USING(subject_id) LEFT JOIN firstmv fm USING(subject_id)
