-- Step 01: find BioProjects containing an Alzheimer's-related signal.
--
-- EXPENSIVE: scans the nested `attributes` column across all of SRA.
-- Run ONCE. Writes a tiny table (just accessions) and returns no result grid.
-- Check the byte estimate in the editor before running.
--
-- Why anchor on disease rather than on controls: disease terms are rare,
-- control terms are not. A regex for `control` would match millions of
-- unrelated runs. Controls are identified in Step 03, within the already
-- narrowed set of studies.
--
-- See docs/query-walkthrough.md for the negation-rewriting logic.

CREATE OR REPLACE TABLE sra_ad.ad_anchor_runs AS
WITH norm AS (
  SELECT
    m.bioproject,
    -- Scope matching to label-ish attribute KEYS, then normalise the text:
    --   1. 'ad libitum' -> blank  (diet descriptions are not Alzheimer's)
    --   2. negated phrases -> ' negctrl '  so that 'non-demented' does not
    --      match the AD pattern via the substring 'demented', and instead
    --      becomes a positive control signal
    REGEXP_REPLACE(
      REGEXP_REPLACE(
        LOWER(IFNULL((SELECT STRING_AGG(CONCAT(a.k, '=', a.v), ' | ')
                      FROM UNNEST(m.attributes) a
                      WHERE REGEXP_CONTAINS(LOWER(a.k),
                        r'disease|diagnos|phenotype|condition|group|genotype|strain|status|cohort|treatment|braak|cerad|subject|clinical|histolog|neuropath|description|title')
                     ), '')),
        r'ad libitum', ' '),
      r'non[-_ ]?(alzheimer[a-z\']*|demented|dementia|tg|transgenic|carrier|ad)\b', ' negctrl '
    ) AS label_blob
  FROM `nih-sra-datastore.sra.metadata` AS m
  WHERE m.consent = 'public'
    AND m.organism IN ('Homo sapiens', 'Mus musculus')
)
SELECT DISTINCT bioproject
FROM norm
WHERE REGEXP_CONTAINS(label_blob,
  -- human diagnosis and neuropathology vocabulary
  r'alzheimer|\bdementia\b|\bbraak\b|\bcerad\b|nia[-_ ]?reagan|neurofibrillary|senile plaque|amyloid|tauopath|\bptau\b'
  -- standard AD / tauopathy mouse models
  r'|app[-_/ ]?ps1|appswe|ps1de9|psen1de9|5x[-_ ]?fad|e4fad|3x[-_ ]?tg|tg2576|tgcrnd8|\bj20\b|r?tg4510|p301[sl]|\bps19\b|\bhtau\b|thy[-_ ]?tau22|tastpm|nl[-_ ]?g[-_ ]?f|samp8'
);

-- CAVEAT: `dementia` matches frontotemporal dementia, and the tau/amyloid
-- terms match Parkinson's and ALS work. This is intentional at the anchor
-- stage (better to over-collect and filter later) but it means the output of
-- Step 03 MUST be checked against study titles. See scripts/enrich_bioprojects.py
