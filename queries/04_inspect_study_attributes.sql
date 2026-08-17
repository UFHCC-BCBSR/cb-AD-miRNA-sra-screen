-- Step 04: read the submitter's literal attribute strings for specific studies.
--
-- CHEAP: reads only sra_ad.ad_runs.
--
-- Aggregates lie. Before committing to a study, look at the raw key=value
-- text. This is how you find out whether disease stage, tissue, or brain
-- region is recoverable from metadata at all, or whether it exists only in
-- the abstract (usually the latter, in which case it becomes a manual
-- review column rather than a filter).

SELECT
  bioproject, sra_study, acc, biosample,
  assay_type, libraryselection, library_name, avgspotlen,
  ROUND(mbytes / 1024, 3) AS gb,
  label_blob,
  attr_blob
FROM sra_ad.ad_runs
WHERE bioproject IN ('PRJNA1333830', 'PRJNA641912', 'PRJNA201039')
ORDER BY bioproject, acc;
