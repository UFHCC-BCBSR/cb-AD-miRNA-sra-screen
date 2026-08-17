# ad-noncoding-expression

Finding public small-RNA / miRNA sequencing datasets that contain both
Alzheimer's disease cases and matched controls, by querying the complete NCBI
SRA metadata catalog in Google BigQuery.

The point is a shortlist of datasets we can download and reprocess, to test
whether specific miRNAs of interest differ in abundance between AD and control.

## Why BigQuery

NCBI mirrors all SRA submitter-supplied metadata into BigQuery as
`nih-sra-datastore.sra.metadata` — one row per sequencing run. That allows
filtering tens of millions of runs on assay type, read length, organism, and
sample attributes in a single query, which the SRA web interface can't do.

The catch, which shapes everything here:

> **`sra.metadata` has no study title, abstract, or description column.**

Every text match is against *sample-level* attributes that submitters typed in
by hand. High precision, poor recall. It finds studies where someone annotated
disease state per sample, and misses studies where Alzheimer's is mentioned
only in the abstract.

## Reproducing it

1. Set up a BigQuery sandbox project — [`docs/bigquery-setup.md`](docs/bigquery-setup.md).
   Free, no credit card, and the free tier covers this whole workflow.
2. Run `queries/00` through `queries/03` in order. 00–02 run once and write
   tables into your own project; 03 is the one you re-run while tuning.
3. Export the Step 03 results to CSV.
4. Add study titles and abstracts, which BigQuery cannot provide:

   ```bash
   python scripts/enrich_bioprojects.py results/03_tissue_classified.csv \
       results/04_enriched_with_titles.csv --email you@example.edu
   ```

   Standard library only — no install step.

5. **Read the titles.** See "Known failure modes" for why this is not optional.

## Layout

```
queries/   BigQuery SQL, numbered in run order
docs/      Setup guide, and an annotated walkthrough of the query logic
scripts/   enrich_bioprojects.py — pulls titles/abstracts from NCBI Entrez
results/   Exported CSVs, one per iteration of the screen
```

## What the screen filters on

`queries/03_screen_and_rank.sql` keeps a study when all of these hold:

| Filter | Rationale |
|---|---|
| Human or mouse, public consent | Scope; dbGaP needs separate authorization |
| ≥3 case-labeled and ≥3 control-labeled samples, non-overlapping | A real contrast, not a single-arm cohort |
| Small-RNA evidence from assay type or library fields | Excludes bulk mRNA-seq that merely mentions AD |
| Median read length 15–120 bp | Mature miRNAs are ~22 nt |
| Not single-cell | Plate-based scRNA-seq is also short-read and high-N |

Case labels cover human diagnosis and neuropathology terms (Braak, CERAD,
amyloid, tauopathy) and the standard mouse models (5XFAD, E4FAD, APP/PS1,
PS19, P301S, 3xTg, Tg2576, others).

Tissue is classified into `brain_prefrontal`, `brain_hippocampal`,
`brain_other`, `csf`, `biofluid`, `cell_model`, or `unannotated`, and is
**reported rather than filtered** — tissue is often recorded only at study
level, so a hard tissue filter silently drops good studies.

## Known failure modes

Every one of these was observed in real output. Read before trusting a result set.

1. **Wrong disease passes.** `dementia` matches frontotemporal dementia; a
   postmortem FTD frontal-lobe study ranked first in one iteration until its
   title was read. Parkinson's studies pass via tau and amyloid vocabulary.
2. **Methods papers pass.** A sequencing-platform evaluation using dementia
   blood samples is indistinguishable from a dementia biomarker study at the
   metadata level.
3. **Attribute keys leak.** A `brainid` field — a brain-bank identifier —
   matched a tissue pattern and made a biofluid study look like brain tissue.
4. **Disease stage is invisible.** "MCI due to AD" and "preclinical AD" both
   contain the AD string. One plasma study read as 28 cases vs 20 controls but
   contained no dementia patients. Biofluid cohorts skew early-stage;
   postmortem brain is late-stage by construction.
5. **Read length is a weak proxy.** Studies caught only by the
   `readlen_inferred` tier are often not small-RNA at all. Prefer
   `assay_declared` and `library_declared`.
6. **Runs ≠ samples.** Check `n_biosamples` against `n_runs`; replicates
   inflate apparent cohort size.
7. **Companion submissions.** One experiment is sometimes deposited as two
   BioProjects (e.g. brain and plasma-EV arms). Not independent replication.

## Results

| File | Studies | What changed |
|---|---|---|
| `results/01_permissive_ad_smallrna.csv` | 47 | Tiered small-RNA evidence |
| `results/02_strict_screen.csv` | 17 | Declared assay only, read-length and size filters |
| `results/03_tissue_classified.csv` | 19 | Tissue split into brain regions vs biofluid |
| `results/04_enriched_with_titles.csv` | 19 | Titles and abstracts from Entrez |

Reading the titles in step 04 removed 6 of the 19 as wrong-target (FTD,
Parkinson's, a platform evaluation, spatial transcriptomics, an ALS/FTD/PD
panel, and a cell line).

Surviving candidates worth pursuing:

**Brain tissue** — PRJNA1333830 (human temporal cortex, 10 AD vs 9 control,
explicitly a miRNA/tau study), PRJNA641912 (human prefrontal cortex, 6 vs 4),
PRJNA700664 (APP/PS1 mouse hippocampus, 6 vs 6, matched miRNA and mRNA),
PRJNA1470739/PRJNA1470741 (PS19 mouse, 8 vs 8, cortex and plasma EV arms).

**Blood / plasma** — PRJNA201039 (48 AD vs 22 control whole blood; the paper
reports 82 miRNAs up and 58 down, so we have a comparison point),
PRJNA1221222 (plasma exosomes, 8 vs 8), PRJNA931267 (7 vs 4), PRJEB52506
(48 samples, but cases are MCI and preclinical), PRJEB49130, PRJNA1262979
(hiPSC-neuron EVs).

Brain tissue is the weak spot — the largest genuine AD brain cohort here is
10 vs 9. PRJNA201039's 48 vs 22 is the only dataset with real statistical
power, and blood miRNA correlates poorly with brain.

## Links

- [SRA in the cloud](https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud/)
- [SRA BigQuery setup](https://www.ncbi.nlm.nih.gov/sra/docs/sra-bigquery/)
- [SRA metadata column reference](https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-based-metadata-table/)
- [BigQuery sandbox](https://cloud.google.com/bigquery/docs/sandbox)
- [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [SRA Toolkit](https://github.com/ncbi/sra-tools)
