-- Episode panel: one row per echocardiography episode for the analysis cohort (patient hash, episode order, day
-- offset from index, reported grade), used for the person-interval models and the hidden Markov model.
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
coh AS (SELECT * FROM adult WHERE subject_id NOT IN (SELECT subject_id FROM pros) AND subject_id NOT IN (SELECT subject_id FROM rheum_echo))
SELECT TO_HEX(MD5(CAST(c.subject_id AS STRING))) pid, s.rn, DATE_DIFF(DATE(s.dt),DATE(c.d1),DAY) t, s.sev, s.nxt, c.sev1, c.n_ep
FROM coh c JOIN seq s USING(subject_id) ORDER BY 1,2
