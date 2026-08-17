# Setting up BigQuery for SRA metadata search

## Cost model, first

BigQuery on-demand pricing bills for **bytes scanned**, not rows returned.
The free tier covers 1 TiB of scanned data per month, which is enough for this
whole workflow if you follow one rule:

> Never `SELECT *` from `nih-sra-datastore.sra.metadata`.

It is a columnar store, so naming columns is what keeps the scan small.
`LIMIT` does **not** reduce bytes scanned — it truncates output after the scan
has already happened. The `attributes` column is the expensive one; touching
it scans that nested field across all of SRA.

The byte estimate appears in the top-right of the query editor before you run
anything. Check it. That is your cost guardrail.

## Sandbox vs billed project

The [BigQuery sandbox](https://cloud.google.com/bigquery/docs/sandbox) needs no
credit card and no billing account, and grants the same free-tier limits
(1 TiB scanned/month, 10 GB storage). It is sufficient for everything in this
repo, because we only *read* public tables and write small intermediate ones.

You need real billing only for:

- Downloading raw sequence data from GCS (requester-pays buckets)
- dbGaP controlled-access data

Sandbox caveat: **tables you create expire after 60 days.** The intermediate
tables in `queries/01` and `queries/02` are cheap to regenerate, so this
matters little, but don't treat them as durable storage.

## Steps

1. Go to [console.cloud.google.com](https://console.cloud.google.com) and
   create a project.

   **If you see a required "Parent resource" field you can't fill:** your
   Google account is managed by Workspace or Cloud Identity (typically an
   institutional or employer address), and you likely lack
   `resourcemanager.projectCreator` on the org. Personal Gmail accounts never
   show this field. Either sign in with a personal account, or ask your
   institution's cloud admin for Project Creator or a folder you can create in.
   Institutional accounts are worth pursuing if you'll eventually download
   data, since many universities have negotiated billing or NIH STRIDES credits.

2. Open BigQuery from the console navigation menu. A "Sandbox" badge confirms
   you're on the free path. No separate API enablement step is needed —
   opening BigQuery turns the API on.

3. Pin the SRA project: in the Explorer panel, click **+ Add data** →
   **Star a project by name** → enter `nih-sra-datastore`. Two datasets appear:

   - `sra` — submitter metadata (`sra.metadata`, one row per run)
   - `sra_tax_analysis_tool` — k-mer taxonomy results, joins on `acc`

4. Smoke test:

   ```sql
   SELECT acc, organism, assay_type, mbases
   FROM `nih-sra-datastore.sra.metadata`
   WHERE organism = 'Homo sapiens' AND consent = 'public'
   LIMIT 10;
   ```

5. Create a dataset to hold intermediate tables. It **must** be in the US
   multi-region to join against `nih-sra-datastore`:

   ```sql
   CREATE SCHEMA IF NOT EXISTS sra_ad OPTIONS(location = 'US');
   ```

## Column names to know

The BigQuery table renames some fields relative to the SRA Run Selector web
UI. The one that catches everyone:

| Run Selector label | BigQuery column |
|---|---|
| LibraryStrategy | `assay_type` |
| LibraryLayout | `librarylayout` |
| LibrarySelection | `libraryselection` |
| LibrarySource | `librarysource` |

There is no `librarystrategy` column. Full reference:
[SRA cloud-based metadata table](https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-based-metadata-table/).

That page was last updated in 2020, so treat the live schema as authoritative:
expand `nih-sra-datastore` → `sra` → `metadata` in the Explorer and open the
**Schema** tab.

## Working effectively in the console

- The red squiggle under a column name is a free dry-run against the real
  schema. If it clears, your column names are valid and you haven't spent
  anything yet.
- Use the **+** button to open a new query tab rather than overwriting one.
  Keep one tab per query file.
- **Save → Save query** prompts for a name *after* you click it. The tab is
  called "Untitled query" because it hasn't been saved yet; saving is what
  names it.
- Every query you run is in **Job history** in the left panel permanently,
  even if you overwrote the editor text. Nothing is lost.
- `sra_ad.ad_runs` resolves against your current project. The fully qualified
  form is `` `your-project-id.sra_ad.ad_runs` ``. Third-party tables like
  `nih-sra-datastore.sra.metadata` always need all three parts.

## Command line

```bash
# install the SDK: https://cloud.google.com/sdk/docs/install
gcloud auth login
gcloud config set project YOUR_PROJECT_ID

# dry run to check cost before committing
bq query --dry_run --nouse_legacy_sql < queries/02_extract_study_runs.sql

# produce an accession list for the SRA Toolkit
bq --format=csv query --nouse_legacy_sql --max_rows=10000 \
  'SELECT acc FROM `nih-sra-datastore.sra.metadata`
   WHERE bioproject = "PRJNA1333830"' | tail -n +2 > accessions.txt
```
