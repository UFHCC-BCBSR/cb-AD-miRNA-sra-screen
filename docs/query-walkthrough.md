# How the queries work

The screen is split across four files rather than written as one query. The
reason is cost, not style.

## Why three stages

Referencing the nested `attributes` column scans it across all of SRA. If the
whole screen were one query, every regex tweak would trigger another full
scan and a handful of iterations would exhaust the monthly free tier.

So:

| Query | Scans | Run how often |
|---|---|---|
| `01_ad_anchor_scan.sql` | All of SRA (expensive) | Once |
| `02_extract_study_runs.sql` | All of SRA (expensive) | Once |
| `03_screen_and_rank.sql` | Your own small table (free) | As often as you like |
| `04_inspect_study_attributes.sql` | Your own small table (free) | Ad hoc |
| `05_total_longread_screen.sql` | Your own small table (free) | As often as you like |

Steps 01 and 02 write tables into `sra_ad` and return no result grid. After
they've run, all iteration happens against a few hundred MB instead of tens of
TB.

## Step 01 — anchor on disease, not on controls

The query finds runs with an Alzheimer's-related signal and returns only the
distinct BioProject accessions.

Anchoring on disease is what bounds the output. A regex for `control` or
`healthy` would match millions of unrelated runs across every field of
biology. Disease terms are rare; control terms are not. So the disease side
does the filtering, and controls are identified later *within* the already
narrowed set of studies.

Two text-normalisation tricks happen before matching:

```sql
REGEXP_REPLACE(..., r'ad libitum', ' ')
```

Without this, any bare `\bad\b` pattern matches diet descriptions. (We
ultimately dropped bare `ad` anyway — it is too promiscuous even scoped to
label keys.)

```sql
REGEXP_REPLACE(..., r'non[-_ ]?(alzheimer[a-z\']*|demented|dementia|tg|...)\b', ' negctrl ')
```

This is the important one. `disease_state=non-demented` contains the substring
`demented`, so a naive AD regex marks controls as cases. Rewriting negated
phrases to a single `negctrl` token before matching removes the false positive
*and* converts it into a positive control signal, since `negctrl` is in the
control pattern.

Matching is scoped to values under label-ish attribute **keys**
(`disease`, `diagnos`, `genotype`, `phenotype`, `condition`, `group`,
`treatment`, `braak`, `cerad`, ...) rather than the whole attribute blob. This
is the single biggest precision gain, and it's what makes short tokens like
`\bwt\b` and `\bcn\b` safe enough to use at all.

## Step 02 — pull every run from the matched studies

Deliberately does *not* filter on small-RNA evidence or on case/control
status. We want the complete run inventory for each candidate study, so that
Step 03 can see the whole design — including control arms and unlabeled arms
that the anchor scan never touched.

This is also where `label_blob`, `attr_blob`, `lib_lc`, and `libsel_lc` are
materialised, so Step 03 never has to re-parse the nested column.

## Step 03 — label, classify, rank

### Disjoint case/control counts

```sql
COUNTIF(ad_hit AND NOT ctrl_hit) AS case_only,
COUNTIF(ctrl_hit AND NOT ad_hit) AS control_only,
COUNTIF(ad_hit AND ctrl_hit)     AS ambiguous
```

An earlier version used `COUNTIF(ad_hit)` and `COUNTIF(ctrl_hit)`, which
overlap: a run matching both patterns was counted twice. You could detect the
problem by summing the two columns and comparing against `n_runs` — several
studies summed to more than their run count.

The usual cause is a study-level annotation stamped onto every sample
(`disease=Alzheimer's disease` present on controls too, or `strain=APP/PS1` on
the WT littermates) alongside a per-sample `disease_state=control`. Disjoint
counting makes the `HAVING`/`WHERE` clause enforce a genuine split, and
`ambiguous` becomes a diagnostic column: nonzero means the study's labels need
manual reading.

### Tiered small-RNA evidence

```
assay_declared      assay_type is miRNA-Seq or ncRNA-Seq       -- trust
library_declared    library_name/libraryselection says small RNA -- trust
attribute_mentioned attributes mention miRNA                    -- weak
readlen_inferred    RNA-Seq/OTHER with median read ≤51 bp       -- guess
```

Many small-RNA libraries are deposited as plain `RNA-Seq`, so read length is
tempting as a proxy: mature miRNAs are ~22 nt, so libraries get sequenced
~50 bp single-end. In practice `readlen_inferred` has a high false-positive
rate — it catches plate-based single-cell data, degraded libraries, and
short-read bulk RNA-seq. Treat it as a lead, not a result.

`libraryselection = 'size fractionation'` is a strong and underused signal.

### Tissue classification, reported not filtered

Tissue is bucketed into `brain_prefrontal`, `brain_hippocampal`,
`brain_other`, `csf`, `biofluid`, `cell_model`, `unannotated`, `other`.

Two design notes:

- Matching is scoped to tissue-ish keys, because an unscoped `serum` pattern
  matches **fetal bovine serum** in culture-media descriptions, making cell
  culture look like a serum biomarker study.
- `prefrontal` is tested before `hippocampal`, so a multi-region study
  classifies as prefrontal. Check `n_hippocampal` to see the split. Note that
  this makes `tissue_verdict` read "brain (other region)" for purely
  hippocampal studies, which is confusing but harmless — the per-region counts
  are correct.

Tissue is **not** in the filter clause. Submitters often record tissue only in
the study description, so a hard tissue filter silently drops real studies.

### Filters worth understanding

```sql
median_readlen BETWEEN 15 AND 120
```
The upper bound matters even for `assay_declared` studies. Several 300 bp
libraries labeled `ncRNA-Seq` appeared in output — honest labels, but long-read
or lncRNA work, not miRNA.

```sql
gb_per_run > 0.015
```
Catches plate-based single-cell data by file size when keyword matching misses
it (a 1152-run study at ~13 MB/run was one cell per "biosample"). Set this too
high and you lose real plasma miRNA libraries, which are legitimately ~20 MB
because the input transcriptome is tiny. 0.05 was too aggressive; 0.015 works.

## Step 04 — read the actual attributes

Aggregates lie. Before committing to a study, look at the submitter's literal
key=value strings. This is cheap and it's how you find out whether disease
stage, tissue, or region is recoverable from metadata at all, or whether it
lives only in the abstract.

## What SQL cannot do

`sra.metadata` has no title, abstract, or description. So the BigQuery screen
cannot tell:

- Alzheimer's from frontotemporal dementia (both match `dementia`)
- a disease study from a methods/platform evaluation using disease samples
- AD dementia from MCI or preclinical AD
- which brain region, when tissue is described only at study level

`scripts/enrich_bioprojects.py` fetches titles and abstracts from Entrez for
whatever the screen returns. **Run it and read the output.** In one iteration,
6 of 19 surviving studies were the wrong target, including the top-ranked one.

Note that the script's `title_mentions_alzheimer` flag keys on a keyword list
that includes `dementia`, `amyloid`, and `neurodegener`, so it returns TRUE for
FTD and Parkinson's studies. It is a triage aid, not a verdict.

## Step 05 — the same tables, a different assay

Steps 01 and 02 contain no assay, `libraryselection` or `library_name` term
anywhere — 01 anchored on *disease* alone, and 02 pulled every run in every
matched BioProject. The only library-construction filter in the workflow is the
`keep` block of Step 03.

That makes the Step 05 pool the whole of `ad_anchor_runs`, not the 19 studies
that survived Step 03. A bulk RNA-seq AD study was anchored and materialised
into `sra_ad.ad_runs` like any other, then dropped for failing
`(n_assay_declared + n_library_declared) > 0`. Step 05 is therefore not a hunt
for a hidden total-RNA arm inside the small-RNA hits — it recovers what Step 03
discarded, and the discard criterion was precisely *absence of small-RNA
evidence*.

Worth noting too that Step 03 filters `agg`, which is `GROUP BY bioproject`.
Those thresholds are study-level. `median_readlen BETWEEN 15 AND 120` is the
median over *all* runs in a study, so a real matched design — small-RNA plus
total-RNA on the same subjects — could have been discarded whole because its
combined median landed outside the window.

Target definition, two ways in:

- a short-read library declared `RANDOM` / `cDNA_randomPriming` /
  `Inverse rRNA`, or with explicit depletion-kit language in the free text; or
- any long-read RNA library that is not small-RNA. Long-read cDNA is
  random-primed or full-length by construction, and PacBio/Nanopore submitters
  routinely leave `libraryselection` blank, so platform alone qualifies.

Two things it does differently from Step 03:

**Case/control are counted within the target arm only.** A study whose AD
labels all sit on its small-RNA runs would otherwise pass on the strength of
runs we are not going to use.

**Nothing is hard-filtered on cohort size.** Expected yield is low, and a
`WHERE case >= 3` returning zero rows tells you nothing about why. Every study
with at least one candidate run comes back, carrying a `verdict` column.

One deliberate omission in the regex: bare `total rna`. Nearly every submitter
writes "total RNA was extracted" in the attributes, polyA libraries included —
it describes the input, not the library. Only kit names, depletion language and
explicit priming language count as evidence.

### The ceiling this cannot break

`consent = 'public'` was applied in Steps 01 and 02. For small-RNA that cost
little. For bulk/total AD brain RNA-seq it is decisive: ROSMAP, MSBB, Mayo —
the AMP-AD studies, which are where the well-powered AD brain transcriptomes
actually live — are dbGaP controlled-access. No public-consent query will ever
return them. A thin Step 05 result is substantially an artifact of that, and
the fix is a dbGaP data access request, not a better regex.

## Complementary approach

Because the metadata screen is precision-biased, the honest complement is a
literature-side search: find AD + miRNA papers in PubMed/GEO, collect their
`PRJNA` accessions, and filter this table with

```sql
WHERE bioproject IN UNNEST(['PRJNA...', 'PRJNA...'])
```

That lets Entrez do the free-text disease matching it is good at, and BigQuery
do the run-level assay and cohort-structure analysis it is good at. Not done
in this repo yet.
