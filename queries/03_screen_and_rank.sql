-- Step 03: label runs, classify tissue, aggregate to studies, filter and rank.
--
-- CHEAP: reads only sra_ad.ad_runs. Iterate on this file freely.
--
-- Output is one row per study with a cohort line, a tissue verdict, tiered
-- small-RNA evidence counts, and direct NCBI links.
--
-- IMPORTANT: this cannot distinguish Alzheimer's from other dementias, or a
-- disease study from a methods paper. Export the results and run
-- scripts/enrich_bioprojects.py before trusting anything here.

WITH pat AS (
  SELECT
    r'alzheimer|\bdementia\b|\bbraak\b|\bcerad\b|nia[-_ ]?reagan|neurofibrillary|senile plaque|amyloid|tauopath|\bptau\b|app[-_/ ]?ps1|appswe|ps1de9|psen1de9|5x[-_ ]?fad|e4fad|3x[-_ ]?tg|tg2576|tgcrnd8|\bj20\b|r?tg4510|p301[sl]|\bps19\b|\bhtau\b|thy[-_ ]?tau22|tastpm|nl[-_ ]?g[-_ ]?f|samp8' AS ad_re,
    r'negctrl|\bcontrol\b|\bctrl\b|\bhealthy\b|\bnormal\b|unaffected|non[-_ ]?demented|\bnd\b|wild[-_ ]?type|\bwt\b|\bntg\b|\bnon[-_ ]?tg\b|littermate|cognitively (?:normal|intact)|\bcn\b' AS ctrl_re,
    r'mirna|micro[-_ ]?rna|small[-_ ]?rna|\bsrna\b|\bsmrna\b|size fractionation|small ncrna' AS srna_re,
    r'single[-_ ]?cell|single[-_ ]?nucle|10x genomics|chromium|smart[-_ ]?seq|drop[-_ ]?seq|cel[-_ ]?seq|scrna|snrna[-_ ]?seq|10xv[23]|plate[-_ ]?based' AS sc_re
),

base AS (
  SELECT
    r.*,
    -- tissue text scoped to tissue-ish KEYS. Unscoped matching lets
    -- 'fetal bovine serum' in a media description masquerade as a serum
    -- biomarker sample.
    IFNULL((SELECT STRING_AGG(kv, ' | ' ORDER BY kv)
            FROM UNNEST(SPLIT(r.attr_blob, ' | ')) AS kv
            WHERE REGEXP_CONTAINS(kv,
              r'^[^=]*(tissue|source_?name|organ|body_?site|biomaterial|cell_?type|cell_?line|isolate|sample_?type|anatom|region|brain)[^=]*=')
           ), '') AS tissue_blob
  FROM sra_ad.ad_runs r
),

labeled AS (
  SELECT
    b.*,
    REGEXP_CONTAINS(b.label_blob, p.ad_re)   AS ad_hit,
    REGEXP_CONTAINS(b.label_blob, p.ctrl_re) AS ctrl_hit,
    REGEXP_EXTRACT(b.label_blob, CONCAT('(', p.ad_re, ')'))   AS ad_term,
    REGEXP_EXTRACT(b.label_blob, CONCAT('(', p.ctrl_re, ')')) AS ctrl_term,
    REGEXP_CONTAINS(CONCAT(b.attr_blob,' ',b.lib_lc,' ',b.libsel_lc), p.sc_re) AS single_cell,

    -- tiered small-RNA evidence, strongest first.
    -- 'readlen_inferred' is a guess with a high false-positive rate; prefer
    -- the two declared tiers.
    CASE
      WHEN LOWER(b.assay_type) IN ('mirna-seq','ncrna-seq') THEN 'assay_declared'
      WHEN REGEXP_CONTAINS(CONCAT(b.lib_lc,' ',b.libsel_lc), p.srna_re) THEN 'library_declared'
      WHEN REGEXP_CONTAINS(b.attr_blob, p.srna_re) THEN 'attribute_mentioned'
      WHEN LOWER(b.assay_type) IN ('rna-seq','other')
           AND b.avgspotlen BETWEEN 15 AND 51 THEN 'readlen_inferred'
      ELSE NULL
    END AS srna_evidence,

    -- tissue class. Order matters: prefrontal is tested first, so a
    -- multi-region study classifies as prefrontal (check n_hippocampal for
    -- the split).
    CASE
      WHEN REGEXP_CONTAINS(b.tissue_blob,
        r'prefrontal|pre[-_ ]?frontal|\bpfc\b|\bdlpfc\b|brodmann|\bba\s?[0-9]|frontal') THEN 'brain_prefrontal'
      WHEN REGEXP_CONTAINS(b.tissue_blob,
        r'hippocamp|entorhinal|dentate|subiculum|\bca1\b|\bca3\b') THEN 'brain_hippocampal'
      WHEN REGEXP_CONTAINS(b.tissue_blob,
        r'temporal|parietal|occipital|cingulate|insula|cortex|cortical|amygdala|cerebell|striat|substantia|thalam|hypothalam|\bbrain\b|cerebr') THEN 'brain_other'
      WHEN REGEXP_CONTAINS(b.tissue_blob, r'\bcsf\b|cerebrospinal') THEN 'csf'
      WHEN REGEXP_CONTAINS(b.tissue_blob,
        r'plasma|\bserum\b|whole blood|\bblood\b|exosom|extracellular vesicle|\bevs?\b|urine|saliva|platelet') THEN 'biofluid'
      WHEN REGEXP_CONTAINS(b.tissue_blob,
        r'ipsc|hipsc|\bsh[-_ ]?sy5y\b|\bhek\b|\bn2a\b|neuro[-_ ]?2a|cell line|primary (?:culture|neuron)|organoid|\bnpc\b') THEN 'cell_model'
      WHEN b.tissue_blob = '' THEN 'unannotated'
      ELSE 'other'
    END AS tissue_class

  FROM base b CROSS JOIN pat p
),

agg AS (
  SELECT
    bioproject, sra_study, organism,
    COUNT(*)                  AS n_runs,
    COUNT(DISTINCT biosample) AS n_biosamples,

    -- disjoint counts: a run matching both patterns lands in `ambiguous`
    -- rather than being counted twice. Nonzero ambiguous means the study's
    -- labels need reading by hand.
    COUNTIF(ad_hit AND NOT ctrl_hit)     AS case_only,
    COUNTIF(ctrl_hit AND NOT ad_hit)     AS control_only,
    COUNTIF(ad_hit AND ctrl_hit)         AS ambiguous,
    COUNTIF(NOT ad_hit AND NOT ctrl_hit) AS unlabeled,

    COUNTIF(srna_evidence = 'assay_declared')      AS n_assay_declared,
    COUNTIF(srna_evidence = 'library_declared')    AS n_library_declared,
    COUNTIF(srna_evidence = 'attribute_mentioned') AS n_attr_mentioned,
    COUNTIF(srna_evidence = 'readlen_inferred')    AS n_readlen_inferred,
    COUNTIF(single_cell) AS n_single_cell,

    COUNTIF(tissue_class = 'brain_prefrontal')  AS n_prefrontal,
    COUNTIF(tissue_class = 'brain_hippocampal') AS n_hippocampal,
    COUNTIF(tissue_class = 'brain_other')       AS n_brain_other,
    COUNTIF(tissue_class LIKE 'brain%')         AS n_brain_any,
    COUNTIF(tissue_class = 'csf')               AS n_csf,
    COUNTIF(tissue_class = 'biofluid')          AS n_biofluid,
    COUNTIF(tissue_class = 'cell_model')        AS n_cell_model,
    COUNTIF(tissue_class = 'unannotated')       AS n_tissue_unannotated,
    COUNTIF(tissue_class = 'other')             AS n_tissue_other,

    ARRAY_AGG(DISTINCT ad_term      IGNORE NULLS) AS ad_terms,
    ARRAY_AGG(DISTINCT ctrl_term    IGNORE NULLS) AS ctrl_terms,
    ARRAY_AGG(DISTINCT tissue_class IGNORE NULLS) AS tissue_classes,
    ARRAY_AGG(DISTINCT assay_type   IGNORE NULLS) AS assay_types,

    APPROX_QUANTILES(avgspotlen, 2)[OFFSET(1)] AS median_readlen,
    ROUND(SUM(mbytes) / COUNT(*) / 1024, 3)    AS gb_per_run,
    ROUND(SUM(mbytes) / 1024, 1)               AS total_gb,
    MAX(releasedate)                           AS latest_release
  FROM labeled
  GROUP BY bioproject, sra_study, organism
),

keep AS (
  SELECT * FROM agg
  WHERE case_only >= 3                                  -- real case arm
    AND control_only >= 3                               -- real control arm
    AND (n_assay_declared + n_library_declared) > 0      -- not just a read-length guess
    AND median_readlen BETWEEN 15 AND 120                -- excludes long-read / lncRNA
    AND gb_per_run > 0.015                               -- excludes single-cell plates
    AND n_single_cell = 0
),

tis_detail AS (
  SELECT bioproject, sra_study,
    STRING_AGG(DISTINCT kv, ', ' ORDER BY kv LIMIT 6) AS tissue_attrs
  FROM labeled l, UNNEST(SPLIT(l.tissue_blob, ' | ')) AS kv
  WHERE l.bioproject IN (SELECT bioproject FROM keep)
    AND kv != '' AND LENGTH(kv) < 70
  GROUP BY 1, 2
)

SELECT
  k.bioproject, k.sra_study, k.organism,

  -- NB: sample count is biosamples, case/control counts are runs. These
  -- differ when a study has replicates per sample. Trust n_biosamples for
  -- cohort size.
  CONCAT(CAST(k.n_biosamples AS STRING), ' samples: ',
         CAST(k.case_only AS STRING), ' AD-like / ',
         CAST(k.control_only AS STRING), ' control',
         IF(k.unlabeled > 0, CONCAT(' / ', CAST(k.unlabeled AS STRING), ' other'), '')
  ) AS cohort_line,

  CASE
    WHEN k.n_prefrontal > 0 THEN 'PREFRONTAL'
    WHEN k.n_brain_any  > 0 THEN 'brain (other region)'
    WHEN k.n_tissue_unannotated = k.n_runs THEN 'unannotated - CHECK STUDY'
    WHEN k.n_csf > 0 OR k.n_biofluid > 0 THEN 'biofluid/CSF only'
    ELSE 'other'
  END AS tissue_verdict,

  k.n_prefrontal, k.n_hippocampal, k.n_brain_other, k.n_brain_any,
  k.n_csf, k.n_biofluid, k.n_cell_model, k.n_tissue_unannotated, k.n_tissue_other,
  t.tissue_attrs,

  LEAST(k.case_only, k.control_only) AS smaller_arm,   -- what actually bounds power
  k.case_only, k.control_only, k.ambiguous, k.unlabeled,
  k.n_runs, k.n_biosamples,
  k.n_assay_declared, k.n_library_declared,
  k.ad_terms, k.ctrl_terms, k.assay_types,
  k.median_readlen, k.gb_per_run, k.total_gb, k.latest_release,

  CONCAT('https://www.ncbi.nlm.nih.gov/bioproject/', k.bioproject) AS bioproject_url,
  CONCAT('https://www.ncbi.nlm.nih.gov/Traces/study/?acc=', k.sra_study) AS run_selector_url

FROM keep k
LEFT JOIN tis_detail t USING (bioproject, sra_study)
ORDER BY k.n_prefrontal DESC, k.n_brain_any DESC, smaller_arm DESC;

-- If results come back too thin, relax the filters in this order:
--   1. drop gb_per_run
--   2. raise median_readlen upper bound to 150
--   3. allow n_attr_mentioned > 0
--   4. lower case_only/control_only to 2
