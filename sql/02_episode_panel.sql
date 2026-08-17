-- Episode panel: one row per echocardiography episode for the analysis cohort (patient hash, episode order, day offset of the
-- episode's defining study from time zero, which is the last study of the index episode (the index episode itself is at 0),
-- day offsets of the episode's first and last studies, number of studies, inpatient flag, reported grade), used for the
-- person-interval models, the visit-process model and the hidden Markov model. Episode construction is identical to
-- 01_analysis_frame.sql (same temporary-table chain), so the two extracts always agree.
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
SELECT TO_HEX(MD5(CAST(c.subject_id AS STRING))) pid, s.rn, IF(s.rn=1, 0, DATE_DIFF(DATE(s.dt),DATE(c.d1_end),DAY)) t, DATE_DIFF(DATE(s.dt_first),DATE(c.d1_end),DAY) t_start, DATE_DIFF(DATE(s.dt_last),DATE(c.d1_end),DAY) t_end, s.n_stud, s.inpt, s.sev, s.nxt, c.sev1, c.n_ep, c.rheum_post
FROM coh c JOIN seq s USING(subject_id) ORDER BY 1,2
