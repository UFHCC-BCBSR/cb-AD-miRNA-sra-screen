-- Step 05: screen the SAME study set for total / long-read RNA-seq with
-- random priming -- the assay the noncoding target actually needs.
--
-- CHEAP: reads only sra_ad.ad_runs. No re-scan of nih-sra-datastore.
--
-- This works without regenerating anything because Step 01 anchored on
-- disease vocabulary alone and Step 02 pulled EVERY run in each matched
-- BioProject with no assay filter. The total-RNA, polyA and long-read runs
-- have been sitting in sra_ad.ad_runs the whole time; Step 03 simply
-- filtered them out at the last stage.
--
-- TWO LIMITS INHERITED FROM 01/02, both important here:
--
--   1. consent = 'public' was applied upstream. The major AD brain bulk
--      RNA-seq resources (ROSMAP, MSBB, Mayo -- the AMP-AD studies) are
--      dbGaP controlled-access and are therefore ABSENT BY CONSTRUCTION.
--      A thin result set is partly an artifact of that, not evidence that
--      no such data exists.
--   2. Recall is still bounded by sample-level attribute text. Studies that
--      mention Alzheimer's only in the abstract were never anchored.
--
-- Unlike Step 03 this does NOT hard-filter to >=3 vs >=3. Expected yield is
-- low, and a hard filter returning zero rows tells you nothing about why.
-- Every study with at least one candidate run is returned, with a `verdict`
-- column saying whether it clears the cohort bar.

WITH pat AS (
  SELECT
    -- unchanged from Step 03
    r'alzheimer|\bdementia\b|\bbraak\b|\bcerad\b|nia[-_ ]?reagan|neurofibrillary|senile plaque|amyloid|tauopath|\bptau\b|app[-_/ ]?ps1|appswe|ps1de9|psen1de9|5x[-_ ]?fad|e4fad|3x[-_ ]?tg|tg2576|tgcrnd8|\bj20\b|r?tg4510|p301[sl]|\bps19\b|\bhtau\b|thy[-_ ]?tau22|tastpm|nl[-_ ]?g[-_ ]?f|samp8' AS ad_re,
    r'negctrl|\bcontrol\b|\bctrl\b|\bhealthy\b|\bnormal\b|unaffected|non[-_ ]?demented|\bnd\b|wild[-_ ]?type|\bwt\b|\bntg\b|\bnon[-_ ]?tg\b|littermate|cognitively (?:normal|intact)|\bcn\b' AS ctrl_re,
    r'mirna|micro[-_ ]?rna|small[-_ ]?rna|\bsrna\b|\bsmrna\b|size fractionation|small ncrna' AS srna_re,
    r'single[-_ ]?cell|single[-_ ]?nucle|10x genomics|chromium|smart[-_ ]?seq|drop[-_ ]?seq|cel[-_ ]?seq|scrna|snrna[-_ ]?seq|10xv[23]|plate[-_ ]?based' AS sc_re,

    -- STRICT random-priming / rRNA-depletion evidence.
    --
    -- Note what is deliberately NOT in here: bare 'total rna'. Almost every
    -- submitter writes "total RNA was extracted" in the attributes, including
    -- for polyA libraries -- it describes the input, not the library. Matching
    -- it would flag essentially every RNA study in SRA. Only kit names,
    -- depletion language, and explicit priming language count.
    r'random[-_ ]?(?:hexamer|primer|priming)|ribo[-_ ]?zero|ribozero|ribo[-_ ]?erase|riboerase|ribo[-_ ]?minus|rrna[-_ ]?deplet|rrna[-_ ]?remov|depletion of (?:cytoplasmic )?rrna|ribosomal rna deplet|inverse rrna|whole[-_ ]?transcriptome|total[-_ ]?rna[-_ ]?seq|smarter[-_ ]?stranded|truseq stranded total' AS totalrna_re,

    -- polyA / oligo-dT: disqualifying for the noncoding target, since polyA
    -- selection strips most lncRNA and all non-polyadenylated species.
    r'poly[-_ ]?a[-_ ]?(?:select|enrich|captur|plus)|polyadenylat|oligo[-_ ]?dt|mrna captur|truseq (?:stranded )?mrna|nebnext.{0,20}poly' AS polya_re,

    -- long-read platforms and protocols
    r'sequel|revio|pacbio|\bsmrt\b|minion|gridion|promethion|nanopore|iso[-_ ]?seq|isoseq|direct[-_ ]?rna[-_ ]?seq|full[-_ ]?length (?:cdna|transcript)|long[-_ ]?read' AS long_re,

    -- Background-strain names, used as a CONDITIONAL control signal only.
    -- Mouse submitters routinely record the control arm as its strain
    -- (genotype=B6129SF2/J) rather than as 'wild type', which ctrl_re misses.
    -- This must never fire on a run that already matched an AD model term:
    -- 5XFAD is maintained ON a C57BL/6 background, so an unconditional match
    -- would push real cases into `ambiguous`. See the ctrl_hit logic below.
    r'c57bl|c57black|\bb6129|129s[0-9]?[/ ]|balb[-_/ ]?c|\bfvb\b|\bicr\b|swiss webster' AS strain_ctrl_re
),

base AS (
  SELECT
    r.*,
    -- tissue text scoped to tissue-ish KEYS, identical to Step 03
    IFNULL((SELECT STRING_AGG(kv, ' | ' ORDER BY kv)
            FROM UNNEST(SPLIT(r.attr_blob, ' | ')) AS kv
            WHERE REGEXP_CONTAINS(kv,
              r'^[^=]*(tissue|source_?name|organ|body_?site|biomaterial|cell_?type|cell_?line|isolate|sample_?type|anatom|region|brain)[^=]*=')
           ), '') AS tissue_blob
  FROM sra_ad.ad_runs r
  WHERE r.organism IN ('Homo sapiens', 'Mus musculus')
),

-- Priming and read length are decided first, because rna_class depends on them.
primed AS (
  SELECT
    b.*,
    p.ad_re, p.ctrl_re, p.srna_re, p.sc_re, p.long_re, p.strain_ctrl_re,

    -- Effective label text. When a run's disease/genotype attributes carry
    -- NEITHER a case nor a control signal, fall back to library_name.
    -- Observed in real output: PRJNA894145 names its files WT3.R1.fq.gz /
    -- AD3.R1.fq.gz, PRJNA1482722 uses WT-107, PRJNA720779 uses
    -- PREC_S11_control / PREC_S12_AD. In each case the study design is in the
    -- filename and nowhere else, and the control arm counted as zero.
    --
    -- Gated on 'attributes said nothing' so that a study with real attribute
    -- labels is never overridden by an incidental filename token.
    IF(REGEXP_CONTAINS(b.label_blob, p.ad_re) OR REGEXP_CONTAINS(b.label_blob, p.ctrl_re),
       b.label_blob,
       CONCAT(b.label_blob, ' ', b.lib_lc)) AS label_eff,

    -- Tiered, strongest first. Declared libraryselection beats free text.
    -- 'Inverse rRNA' is kept as its own tier: it is a depletion method, and
    -- in practice always implies random priming, but it is worth seeing
    -- which studies got there by depletion vs by declaring RANDOM outright.
    CASE
      WHEN b.libsel_lc IN ('random', 'random pcr', 'cdna_randompriming')
        THEN 'random_declared'
      WHEN b.libsel_lc LIKE 'inverse rrna%'
        THEN 'ribodep_declared'
      WHEN b.libsel_lc IN ('polya', 'oligo-dt', 'cdna_oligo_dt')
        THEN 'polya_declared'
      -- free-text fallback: library_name is a deliberate submitter label,
      -- attr_blob is noisier but the regex above is strict enough to use it
      WHEN REGEXP_CONTAINS(CONCAT(b.lib_lc, ' ', b.attr_blob), p.totalrna_re)
        THEN 'random_or_ribodep_text'
      WHEN REGEXP_CONTAINS(CONCAT(b.lib_lc, ' ', b.attr_blob), p.polya_re)
        THEN 'polya_text'
      WHEN b.libsel_lc IN ('cdna', 'unspecified', 'other', '')
        THEN 'unspecified'
      ELSE 'other_selection'
    END AS priming_class,

    -- Long read. Platform is authoritative; avgspotlen is a weak tiebreak
    -- only, since 2x250 Illumina also lands around 500.
    CASE
      WHEN UPPER(IFNULL(b.platform, '')) IN ('PACBIO_SMRT', 'OXFORD_NANOPORE')
        THEN 'long_declared'
      WHEN REGEXP_CONTAINS(CONCAT(LOWER(IFNULL(b.instrument, '')), ' ', b.lib_lc, ' ', b.attr_blob), p.long_re)
        THEN 'long_text'
      WHEN b.avgspotlen >= 600 THEN 'long_inferred'
      ELSE 'short_read'
    END AS read_class

  FROM base b CROSS JOIN pat p
),

labeled AS (
  SELECT
    q.*,
    REGEXP_CONTAINS(q.label_eff, q.ad_re) AS ad_hit,

    -- Two ways to be a control: the explicit vocabulary, or -- only when no
    -- AD model term is present on the same run -- a bare background-strain
    -- name. The NOT(ad_re) guard is what keeps 5XFAD-on-C57BL/6 counted as a
    -- case rather than dumped into `ambiguous`.
    (REGEXP_CONTAINS(q.label_eff, q.ctrl_re)
     OR (NOT REGEXP_CONTAINS(q.label_eff, q.ad_re)
         AND REGEXP_CONTAINS(q.label_eff, q.strain_ctrl_re))) AS ctrl_hit,

    REGEXP_EXTRACT(q.label_eff, CONCAT('(', q.ad_re, ')'))   AS ad_term,
    COALESCE(REGEXP_EXTRACT(q.label_eff, CONCAT('(', q.ctrl_re, ')')),
             IF(NOT REGEXP_CONTAINS(q.label_eff, q.ad_re),
                REGEXP_EXTRACT(q.label_eff, CONCAT('(', q.strain_ctrl_re, ')')),
                NULL)) AS ctrl_term,

    -- librarysource is the reliable single-cell tell. The keyword regex alone
    -- let PRJEB54589 (10x microglia) and PRJNA1117576 (snRNA-seq of FFPE)
    -- through as 'PASSES' in the first run of this query.
    (REGEXP_CONTAINS(CONCAT(q.attr_blob, ' ', q.lib_lc, ' ', q.libsel_lc), q.sc_re)
     OR LOWER(IFNULL(q.librarysource, '')) = 'transcriptomic single cell') AS single_cell,

    CASE
      -- small-RNA claims a run first, so a size-fractionated library is never
      -- mistaken for a total-RNA one
      WHEN LOWER(q.assay_type) IN ('mirna-seq', 'ncrna-seq')
           OR REGEXP_CONTAINS(CONCAT(q.lib_lc, ' ', q.libsel_lc), q.srna_re)
        THEN 'small_rna'
      WHEN LOWER(IFNULL(q.librarysource, '')) NOT LIKE 'transcriptomic%'
           AND LOWER(IFNULL(q.assay_type, '')) NOT IN ('rna-seq', 'other')
        THEN 'not_rna'
      WHEN q.priming_class IN ('random_declared', 'ribodep_declared', 'random_or_ribodep_text')
        THEN 'total_rna_random'
      WHEN q.priming_class IN ('polya_declared', 'polya_text')
        THEN 'polya_mrna'
      ELSE 'rna_unspecified'
    END AS rna_class,

    -- tissue class, identical to Step 03. Prefrontal tested first.
    CASE
      WHEN REGEXP_CONTAINS(q.tissue_blob,
        r'prefrontal|pre[-_ ]?frontal|\bpfc\b|\bdlpfc\b|brodmann|\bba\s?[0-9]|frontal') THEN 'brain_prefrontal'
      WHEN REGEXP_CONTAINS(q.tissue_blob,
        r'hippocamp|entorhinal|dentate|subiculum|\bca1\b|\bca3\b') THEN 'brain_hippocampal'
      WHEN REGEXP_CONTAINS(q.tissue_blob,
        r'temporal|parietal|occipital|cingulate|insula|cortex|cortical|amygdala|cerebell|striat|substantia|thalam|hypothalam|\bbrain\b|cerebr') THEN 'brain_other'
      WHEN REGEXP_CONTAINS(q.tissue_blob, r'\bcsf\b|cerebrospinal') THEN 'csf'
      WHEN REGEXP_CONTAINS(q.tissue_blob,
        r'plasma|\bserum\b|whole blood|\bblood\b|exosom|extracellular vesicle|\bevs?\b|urine|saliva|platelet') THEN 'biofluid'
      WHEN REGEXP_CONTAINS(q.tissue_blob,
        r'ipsc|hipsc|\bsh[-_ ]?sy5y\b|\bhek\b|\bn2a\b|neuro[-_ ]?2a|cell line|primary (?:culture|neuron)|organoid|\bnpc\b') THEN 'cell_model'
      WHEN q.tissue_blob = '' THEN 'unannotated'
      ELSE 'other'
    END AS tissue_class

  FROM primed q
),

-- THE TARGET DEFINITION. Two ways in:
--   a) a short-read library that is randomly primed or rRNA-depleted, or
--   b) any long-read RNA library that is not small-RNA
-- Long-read cDNA is random-primed or full-length by construction, so it
-- qualifies even when libraryselection is left unspecified -- which it
-- usually is on PacBio and Nanopore submissions.
flagged AS (
  SELECT
    l.*,
    (
      l.rna_class = 'total_rna_random'
      OR (l.read_class IN ('long_declared', 'long_text')
          AND l.rna_class NOT IN ('small_rna', 'not_rna'))
    ) AND NOT l.single_cell AS is_target
  FROM labeled l
),

agg AS (
  SELECT
    bioproject, sra_study, organism,

    COUNT(*)                              AS n_runs_total,
    COUNTIF(is_target)                    AS n_target_runs,
    COUNT(DISTINCT IF(is_target, biosample, NULL)) AS n_target_biosamples,

    -- Case/control counted WITHIN the target arm only. Counting across the
    -- whole study would let a project whose AD labels live entirely on a
    -- small-RNA arm pass on the strength of runs we are not going to use.
    COUNTIF(is_target AND ad_hit AND NOT ctrl_hit) AS case_only,
    COUNTIF(is_target AND ctrl_hit AND NOT ad_hit) AS control_only,
    COUNTIF(is_target AND ad_hit AND ctrl_hit)     AS ambiguous,
    COUNTIF(is_target AND NOT ad_hit AND NOT ctrl_hit) AS unlabeled,

    -- full library inventory, so you can see what else the study contains
    COUNTIF(rna_class = 'total_rna_random')  AS n_total_rna_random,
    COUNTIF(rna_class = 'polya_mrna')        AS n_polya_mrna,
    COUNTIF(rna_class = 'small_rna')         AS n_small_rna,
    COUNTIF(rna_class = 'rna_unspecified')   AS n_rna_unspecified,
    COUNTIF(rna_class = 'not_rna')           AS n_not_rna,
    COUNTIF(read_class IN ('long_declared', 'long_text')) AS n_long_read,
    COUNTIF(read_class = 'long_inferred')    AS n_long_inferred,
    COUNTIF(single_cell)                     AS n_single_cell,

    COUNTIF(is_target AND priming_class = 'random_declared')        AS n_random_declared,
    COUNTIF(is_target AND priming_class = 'ribodep_declared')       AS n_ribodep_declared,
    COUNTIF(is_target AND priming_class = 'random_or_ribodep_text') AS n_random_text,
    COUNTIF(is_target AND priming_class = 'unspecified')            AS n_priming_unspecified,

    COUNTIF(is_target AND tissue_class = 'brain_prefrontal')  AS n_prefrontal,
    COUNTIF(is_target AND tissue_class = 'brain_hippocampal') AS n_hippocampal,
    COUNTIF(is_target AND tissue_class = 'brain_other')       AS n_brain_other,
    COUNTIF(is_target AND tissue_class LIKE 'brain%')         AS n_brain_any,
    COUNTIF(is_target AND tissue_class = 'csf')               AS n_csf,
    COUNTIF(is_target AND tissue_class = 'biofluid')          AS n_biofluid,
    COUNTIF(is_target AND tissue_class = 'cell_model')        AS n_cell_model,
    COUNTIF(is_target AND tissue_class = 'unannotated')       AS n_tissue_unannotated,

    -- Instrumentation. n_single_cell was computed but never selected in the
    -- first version of this query, which is exactly why PRJEB54589 (a 10x
    -- microglia study) could rank as 'PASSES - real contrast' unnoticed.
    -- librarysource is surfaced for the same reason: it is the field the
    -- single-cell exclusion relies on, and on ENA-originated submissions it
    -- is not always set to 'TRANSCRIPTOMIC SINGLE CELL'.
    ARRAY_AGG(DISTINCT IF(is_target, librarysource, NULL) IGNORE NULLS) AS target_librarysources,
    ARRAY_AGG(DISTINCT IF(is_target, libraryselection, NULL) IGNORE NULLS) AS target_selections,
    ARRAY_AGG(DISTINCT IF(is_target, assay_type,       NULL) IGNORE NULLS) AS target_assay_types,
    ARRAY_AGG(DISTINCT IF(is_target, platform,         NULL) IGNORE NULLS) AS target_platforms,
    ARRAY_AGG(DISTINCT IF(is_target, instrument,       NULL) IGNORE NULLS) AS target_instruments,
    ARRAY_AGG(DISTINCT IF(is_target, library_name,     NULL) IGNORE NULLS LIMIT 8) AS target_library_names,
    ARRAY_AGG(DISTINCT ad_term   IGNORE NULLS) AS ad_terms,
    ARRAY_AGG(DISTINCT ctrl_term IGNORE NULLS) AS ctrl_terms,

    APPROX_QUANTILES(IF(is_target, avgspotlen, NULL), 2)[OFFSET(1)] AS median_target_readlen,
    ROUND(SUM(IF(is_target, mbytes, 0)) / 1024, 1) AS target_gb,
    MAX(releasedate) AS latest_release
  FROM flagged
  GROUP BY bioproject, sra_study, organism
),

keep AS (
  SELECT * FROM agg WHERE n_target_runs > 0
),

tis_detail AS (
  SELECT bioproject, sra_study,
    STRING_AGG(DISTINCT kv, ', ' ORDER BY kv LIMIT 6) AS tissue_attrs
  FROM flagged f, UNNEST(SPLIT(f.tissue_blob, ' | ')) AS kv
  WHERE f.is_target
    AND f.bioproject IN (SELECT bioproject FROM keep)
    AND kv != '' AND LENGTH(kv) < 70
  GROUP BY 1, 2
)

SELECT
  -- Read this column first. Nothing below it means much until you know
  -- whether the study actually has two arms.
  CASE
    WHEN k.case_only >= 3 AND k.control_only >= 3 THEN 'PASSES - real contrast'
    WHEN k.case_only >= 2 AND k.control_only >= 2 THEN 'thin - 2v2, underpowered'
    WHEN k.unlabeled = k.n_target_runs           THEN 'target arm unlabeled - CHECK BY HAND'
    ELSE 'single-arm or too small'
  END AS verdict,

  k.bioproject, k.sra_study, k.organism,

  CONCAT(CAST(k.n_target_biosamples AS STRING), ' target samples: ',
         CAST(k.case_only AS STRING), ' AD-like / ',
         CAST(k.control_only AS STRING), ' control',
         IF(k.unlabeled > 0, CONCAT(' / ', CAST(k.unlabeled AS STRING), ' other'), '')
  ) AS cohort_line,

  -- how we decided this was total / long-read, in descending trustworthiness
  CASE
    WHEN k.n_random_declared > 0  THEN 'libraryselection = RANDOM'
    WHEN k.n_ribodep_declared > 0 THEN 'libraryselection = Inverse rRNA'
    WHEN k.n_long_read > 0        THEN 'long-read platform'
    WHEN k.n_random_text > 0      THEN 'free text only - VERIFY'
    ELSE 'weak'
  END AS evidence,

  CASE
    WHEN k.n_prefrontal > 0 THEN 'PREFRONTAL'
    WHEN k.n_brain_any  > 0 THEN 'brain (other region)'
    WHEN k.n_tissue_unannotated = k.n_target_runs THEN 'unannotated - CHECK STUDY'
    WHEN k.n_csf > 0 OR k.n_biofluid > 0 THEN 'biofluid/CSF only'
    ELSE 'other'
  END AS tissue_verdict,

  LEAST(k.case_only, k.control_only) AS smaller_arm,
  k.n_target_runs, k.n_target_biosamples, k.n_runs_total,
  k.n_total_rna_random, k.n_long_read, k.n_long_inferred,
  k.n_polya_mrna, k.n_small_rna, k.n_rna_unspecified, k.n_not_rna,
  k.n_random_declared, k.n_ribodep_declared, k.n_random_text, k.n_priming_unspecified,
  k.n_prefrontal, k.n_hippocampal, k.n_brain_other, k.n_brain_any,
  k.n_csf, k.n_biofluid, k.n_cell_model, k.n_tissue_unannotated,
  t.tissue_attrs,
  k.case_only, k.control_only, k.ambiguous, k.unlabeled,
  k.n_single_cell,
  k.target_librarysources, k.target_selections, k.target_assay_types, k.target_platforms,
  k.target_instruments, k.target_library_names,
  k.ad_terms, k.ctrl_terms,
  k.median_target_readlen,

  -- total RNA runs are 10-100x the size of small-RNA runs. This is the
  -- download you are signing up for, per study.
  k.target_gb,
  k.latest_release,

  CONCAT('https://www.ncbi.nlm.nih.gov/bioproject/', k.bioproject) AS bioproject_url,
  CONCAT('https://www.ncbi.nlm.nih.gov/Traces/study/?acc=', k.sra_study) AS run_selector_url

FROM keep k
LEFT JOIN tis_detail t USING (bioproject, sra_study)
ORDER BY
  CASE WHEN k.case_only >= 3 AND k.control_only >= 3 THEN 0 ELSE 1 END,
  k.n_prefrontal DESC, k.n_brain_any DESC, smaller_arm DESC, k.n_target_runs DESC;

-- READ BEFORE TRUSTING:
--
-- 1. `evidence = 'free text only - VERIFY'` is the tier to distrust. Confirm
--    against the study's library protocol before committing.
-- 2. `n_priming_unspecified > 0` on a long-read study is expected, not a
--    problem -- PacBio and Nanopore submitters routinely leave
--    libraryselection blank.
-- 3. A study with both n_small_rna > 0 and n_target_runs > 0 is a matched
--    design (same subjects, both assays). Those are the most valuable hits
--    here and are worth checking by hand even if the verdict is 'thin'.
-- 4. Titles still are not in BigQuery. Export and run
--    scripts/enrich_bioprojects.py, exactly as for Step 03.
--
-- If this returns too little, relax in this order:
--   1. allow read_class = 'long_inferred' into is_target
--   2. treat rna_class = 'rna_unspecified' as a candidate when
--      librarysource = 'TRANSCRIPTOMIC' and assay_type = 'RNA-Seq'
--      (this is where most under-annotated total-RNA studies hide)
--   3. drop the organism filter to include other model species
--   4. accept that the answer is upstream: the AMP-AD bulk RNA-seq cohorts
--      are dbGaP-controlled and no public-consent query will ever see them
