-- Step 02: pull EVERY run belonging to the studies found in Step 01.
--
-- EXPENSIVE: second and last full scan. Run ONCE.
--
-- Deliberately does NOT filter on small-RNA evidence or case/control status.
-- We want the complete run inventory per study so Step 03 can see the whole
-- design, including control and unlabeled arms the anchor scan never touched.
--
-- This also materialises the parsed text blobs so Step 03 never re-parses
-- the nested `attributes` column. Everything downstream reads this table
-- and is therefore effectively free.

CREATE OR REPLACE TABLE sra_ad.ad_runs AS
SELECT
  m.bioproject, m.sra_study, m.acc, m.biosample, m.sample_acc, m.organism,
  m.assay_type, m.librarysource, m.libraryselection, m.librarylayout,
  m.library_name, m.platform, m.instrument, m.avgspotlen, m.mbases, m.mbytes,
  m.releasedate,

  LOWER(IFNULL(m.library_name, ''))     AS lib_lc,
  LOWER(IFNULL(m.libraryselection, '')) AS libsel_lc,

  -- full attribute text, for tissue and assay keyword matching
  LOWER(IFNULL((SELECT STRING_AGG(CONCAT(a.k, '=', a.v), ' | ' ORDER BY a.k)
                FROM UNNEST(m.attributes) a), '')) AS attr_blob,

  -- label text: scoped to disease/genotype-ish keys and normalised exactly as
  -- in Step 01, so the same patterns behave identically here
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
  AND m.bioproject IN (SELECT bioproject FROM sra_ad.ad_anchor_runs);
