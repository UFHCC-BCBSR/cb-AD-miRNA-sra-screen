#!/usr/bin/env python3
"""
Enrich a BigQuery-exported study table with BioProject titles and descriptions
from NCBI Entrez, producing a GREIN-style review sheet.

Usage:
    python enrich_bioprojects.py results.csv enriched.csv
    python enrich_bioprojects.py results.csv enriched.csv --email you@ufl.edu --api-key KEY

Input CSV must have a 'bioproject' column (e.g. the export from your
BigQuery query). All existing columns are preserved; new columns are appended.

Notes on rate limits: NCBI allows 3 requests/sec without an API key,
10/sec with one. Get a free key at
https://www.ncbi.nlm.nih.gov/account/settings/ and pass --api-key.
"""

import argparse
import csv
import json
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CHUNK = 100  # accessions per esearch/esummary round trip


def fetch_json(endpoint, params, retries=3):
    """GET an E-utilities endpoint and parse the JSON response."""
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "sra-ad-screen/1.0"}
            )
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except Exception as err:  # noqa: BLE001 - network errors vary
            last_err = err
            time.sleep(2 * (attempt + 1))
    print(f"  ! failed after {retries} tries: {last_err}", file=sys.stderr)
    return None


def base_params(args):
    p = {"retmode": "json", "tool": "sra-ad-screen"}
    if args.email:
        p["email"] = args.email
    if args.api_key:
        p["api_key"] = args.api_key
    return p


def accessions_to_uids(accessions, args, pause):
    """Map BioProject accessions (PRJNA…/PRJEB…) to Entrez UIDs."""
    uids = []
    for i in range(0, len(accessions), CHUNK):
        chunk = accessions[i : i + CHUNK]
        term = " OR ".join(f"{acc}[Project Accession]" for acc in chunk)
        params = base_params(args)
        params.update({"db": "bioproject", "term": term, "retmax": str(len(chunk) * 2)})
        data = fetch_json("esearch.fcgi", params)
        if data:
            found = data.get("esearchresult", {}).get("idlist", [])
            uids.extend(found)
            print(f"  esearch {i + 1}-{i + len(chunk)}: {len(found)} UIDs")
        time.sleep(pause)
    return uids


def uids_to_records(uids, args, pause):
    """Pull title/description per UID, keyed by project accession."""
    out = {}
    for i in range(0, len(uids), CHUNK):
        chunk = uids[i : i + CHUNK]
        params = base_params(args)
        params.update({"db": "bioproject", "id": ",".join(chunk)})
        data = fetch_json("esummary.fcgi", params)
        if data:
            result = data.get("result", {})
            for uid in result.get("uids", []):
                rec = result.get(uid, {})
                acc = rec.get("project_acc", "")
                if not acc:
                    continue
                out[acc] = {
                    "study_title": (rec.get("project_title") or "").strip(),
                    "study_description": " ".join(
                        (rec.get("project_description") or "").split()
                    )[:1200],
                    "submitter": (rec.get("submitter_organization") or "").strip(),
                    "registration_date": (rec.get("registration_date") or "").strip(),
                    "project_data_type": (rec.get("project_data_type") or "").strip(),
                }
            print(f"  esummary {i + 1}-{i + len(chunk)}: {len(out)} records so far")
        time.sleep(pause)
    return out


NEW_COLS = [
    "study_title",
    "study_description",
    "submitter",
    "registration_date",
    "project_data_type",
    "title_mentions_alzheimer",
    "title_mentions_mirna",
]

AD_WORDS = ("alzheimer", "dementia", "amyloid", "tauopath", "neurodegener")
MIR_WORDS = ("mirna", "microrna", "micro-rna", "small rna", "smallrna", "ncrna")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile", help="CSV with a 'bioproject' column")
    ap.add_argument("outfile", help="output CSV path")
    ap.add_argument("--email", default="", help="contact email (NCBI courtesy)")
    ap.add_argument("--api-key", default="", help="NCBI API key (raises rate limit)")
    args = ap.parse_args()

    pause = 0.11 if args.api_key else 0.34

    with open(args.infile, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit("input CSV is empty")
    if "bioproject" not in rows[0]:
        sys.exit("input CSV has no 'bioproject' column")

    accessions = sorted({r["bioproject"].strip() for r in rows if r["bioproject"].strip()})
    print(f"{len(rows)} rows, {len(accessions)} distinct BioProjects")

    print("Resolving accessions to UIDs...")
    uids = accessions_to_uids(accessions, args, pause)
    print(f"Fetching summaries for {len(uids)} UIDs...")
    records = uids_to_records(uids, args, pause)

    missing = [a for a in accessions if a not in records]
    if missing:
        print(f"No summary returned for: {', '.join(missing)}", file=sys.stderr)

    fieldnames = list(rows[0].keys()) + [c for c in NEW_COLS if c not in rows[0]]
    with open(args.outfile, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            rec = records.get(row["bioproject"].strip(), {})
            row.update({c: rec.get(c, "") for c in NEW_COLS[:5]})
            blob = f"{rec.get('study_title', '')} {rec.get('study_description', '')}".lower()
            row["title_mentions_alzheimer"] = str(any(w in blob for w in AD_WORDS))
            row["title_mentions_mirna"] = str(any(w in blob for w in MIR_WORDS))
            writer.writerow(row)

    print(f"\nWrote {args.outfile}")
    print("Check the two title_mentions_* columns first: a study whose own")
    print("title never mentions Alzheimer's or small RNA is the likeliest")
    print("false positive from attribute-only matching.")


if __name__ == "__main__":
    main()
