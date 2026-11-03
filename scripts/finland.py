"""
Finland pipeline: BioColour taxa -> GBIF records in Finland -> biogeographical provinces.

Answers Q1 (which provinces has a taxon been found in?) and Q2 (which taxa have been
found in a province?) for the 21 eliömaakunnat. Replaces the worldwide WGSRPD outputs
of join_wcvp.py / fungi_gbif.py for this project. Plan: plans/next_steps_finland_region.md.

Run the stages in order (from anywhere; paths default to the project's data/ folder):

  0. python scripts/get_provinces.py          -> data/boundaries/provinces.gpkg (once)

  1. python scripts/finland.py match
       Matches every taxon in data/taxa.csv against the GBIF Backbone (no account needed):
       vascular plants in Plantae, fungi + lichens in Fungi, the alga with no kingdom.
       Genus-only taxa ("Salix spp.") match the genus; infraspecific and cultivar taxa match
       their species. Insects are excluded. PROXIES adds extra stand-in species rows.
       Also asks GBIF how many Finnish records each key has (n_fi_records, for review only).
       Writes data/out/finland_name_matches.csv. Check rows flagged REVIEW. To fix one, edit
       gbif_key (key to download), match_on (speciesKey | taxonKey | genus) and match_value
       (how downloaded records are assigned), then set check=manual; manual rows are kept
       when 'match' is rerun. LUMPED and SYNONYM_NOTES below handle known synonym issues.

  2. python scripts/finland.py download
       One GBIF occurrence download (SIMPLE_CSV) for all ok/manual keys, Finland + Åland
       (COUNTRY in FI, AX; GBIF files Åland records under AX), with coordinates, no
       geospatial issues, status PRESENT. Needs GBIF_USER, GBIF_PWD, GBIF_EMAIL in .env at
       the project root (or the environment). Polls until ready, then saves
       data/out/finland_gbif_occurrences.zip and data/out/finland_gbif_download.json (DOI to
       cite). If interrupted, rerun: it resumes the same download.

  3. python scripts/finland.py regions
       Places records in provinces (EPSG:3067; points just outside any polygon are snapped
       to the nearest province within --snap-m) and aggregates per taxon x province.
       No record threshold is applied to the tables; --min-records is only used for the
       "ge" columns and the report. Writes to data/out/:
         finland_taxon_regions.csv    one row per taxon x province with >= 1 record
         finland_taxa.csv             one row per BioColour taxon (all of them)
         finland_province_summary.csv one row per province (all 21)
         finland_report.txt           numbers for the write-up
       If data/out/species_regions.csv (join_wcvp.py) exists, plants get wcvp_fi_status
       (native / introduced in Finland per WCVP).
       If the FinBIF checklist outputs exist (python scripts/finbif_checklist.py, run it
       first), its per-province status is outer-joined in: finland_taxon_regions.csv gets
       checklist_status / checklist_category and source (gbif | checklist | both | none),
       and finland_taxa.csv gets Finnish name, origin in Finland, Red List, invasive,
       protection, habitats, etc.

Options for stage 3: --min-records 3  --max-uncertainty-m 10000  --snap-m 5000  --use-proxies
"""
import argparse
import csv
import json
import os
import sys
import time
import zipfile
from collections import Counter
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"
API = "https://api.gbif.org/v1"

# BioColour group -> GBIF kingdom to match in (None: don't constrain, don't check).
KINGDOM = {"vascular plant": "Plantae", "fungus": "Fungi", "lichen": "Fungi", "alga": None}
EXCLUDED_GROUPS = {"insect": "insect, not a plant or fungus"}
# taxon_id -> stand-in species, used only with `regions --use-proxies`.
PROXIES = {14: "Aronia melanocarpa"}   # Aronia mitschurinii is a cultivated hybrid
# Names the GBIF Backbone lumps into another species. Matching on speciesKey would count
# the other species too, so these count only records identified under this name (taxonKey).
LUMPED = {110: "GBIF lumps Tagetes patula into T. erecta; only records identified as "
               "T. patula are counted"}
# Notes for synonyms worth explaining in the write-up.
SYNONYM_NOTES = {91: "Finnish bracken is P. pinetorum (P. aquilinum subsp. pinetorum); "
                     "P. aquilinum records in Finland are used"}
USABLE = ["ok", "manual"]
FUNGAL = ["fungus", "lichen"]

MATCH_FILE = "finland_name_matches.csv"
ZIP_FILE = "finland_gbif_occurrences.zip"
META_FILE = "finland_gbif_download.json"


def display_name(t):
    if t["rank"] in ("cultivar", "infraspecific") or not isinstance(t["binomial"], str):
        return t["sci_name_raw"]
    return t["binomial"]


# ---------------------------------------------------------------- stage 1
def gbif_match(session, name, kingdom, rank):
    params = {"name": name, "rank": rank}   # without rank, bare genus names match HIGHERRANK
    if kingdom:
        params["kingdom"] = kingdom
    r = session.get(f"{API}/species/match", params=params, timeout=30)
    r.raise_for_status()
    return r.json()


def fi_count(session, key):
    """Finnish + Åland records with clean coordinates, or None if GBIF doesn't answer."""
    params = [("taxonKey", key), ("country", "FI"), ("country", "AX"),
              ("hasCoordinate", "true"), ("hasGeospatialIssue", "false"),
              ("occurrenceStatus", "PRESENT"), ("limit", 0)]
    for attempt in range(2):
        try:
            r = session.get(f"{API}/occurrence/search", params=params, timeout=20)
            if r.ok:
                return r.json().get("count")
        except Exception:
            pass
        time.sleep(2)
    return None


def match_row(t, m, query, level, expect_rank, kingdom):
    """One row of finland_name_matches.csv.

    gbif_key is the accepted key used in the download request (TAXON_KEY includes
    synonyms and descendants). match_on / match_value say how records in the download
    are assigned to this taxon: speciesKey (default), taxonKey (LUMPED names) or genus.
    """
    rank = m.get("rank")
    if expect_rank == "GENUS":
        key = m.get("genusKey") or m.get("acceptedUsageKey") or m.get("usageKey")
        accepted = m.get("genus")
        match_on, match_value = "genus", accepted
    else:
        key = m.get("speciesKey") or m.get("acceptedUsageKey") or m.get("usageKey")
        accepted = m.get("species")
        match_on, match_value = "speciesKey", key
    ok = (m.get("matchType") == "EXACT" and rank == expect_rank
          and m.get("status") in ("ACCEPTED", "SYNONYM")
          and (kingdom is None or m.get("kingdom") == kingdom))
    notes = []
    if m.get("status") == "SYNONYM":
        notes.append(f"GBIF synonym of {accepted}")
    if level == "species" and t["taxon_id"] in LUMPED:
        match_on, match_value = "taxonKey", m.get("usageKey")
        notes.append(LUMPED[t["taxon_id"]])
    if level == "species" and t["taxon_id"] in SYNONYM_NOTES:
        notes.append(SYNONYM_NOTES[t["taxon_id"]])
    return {
        "taxon_id": t["taxon_id"], "name": display_name(t), "group": t["group"],
        "query_name": query, "match_level": level, "gbif_key": key,
        "match_on": match_on, "match_value": match_value,
        "gbif_name": m.get("scientificName"), "gbif_accepted": accepted,
        "gbif_usage_key": m.get("usageKey"), "gbif_rank": rank,
        "gbif_status": m.get("status"), "gbif_match_type": m.get("matchType"),
        "gbif_confidence": m.get("confidence"), "gbif_kingdom": m.get("kingdom"),
        "gbif_family": m.get("family"),
        "check": "ok" if ok else "REVIEW", "note": "; ".join(notes),
    }


def stage_match(taxa_path, out, counts):
    import requests

    taxa = pd.read_csv(taxa_path)
    path = out / MATCH_FILE
    manual = pd.DataFrame()
    if path.exists():
        old = pd.read_csv(path)
        manual = old[old["check"] == "manual"]
        if len(manual):
            print(f"Keeping {len(manual)} manual row(s) from {path.name}.")
    manual_ids = set(zip(manual.get("taxon_id", []), manual.get("match_level", [])))

    s = requests.Session()
    rows = []
    for t in taxa.to_dict("records"):
        if t["group"] in EXCLUDED_GROUPS:
            rows.append({"taxon_id": t["taxon_id"], "name": display_name(t),
                         "group": t["group"], "match_level": "excluded",
                         "check": "excluded", "note": EXCLUDED_GROUPS[t["group"]]})
            continue
        kingdom = KINGDOM.get(t["group"])
        if t["rank"] == "genus":
            jobs = [(t["genus"], "genus", "GENUS")]
        else:
            jobs = [(t["binomial"], "species", "SPECIES")]
        if t["taxon_id"] in PROXIES:
            jobs.append((PROXIES[t["taxon_id"]], "proxy", "SPECIES"))
        for query, level, expect in jobs:
            if (t["taxon_id"], level) in manual_ids:
                continue
            m = gbif_match(s, query, kingdom, expect)
            row = match_row(t, m, query, level, expect, kingdom)
            if t["rank"] in ("infraspecific", "cultivar") and level == "species":
                row["note"] = (row["note"] + "; " if row["note"] else "") + \
                    f"{t['rank']} '{t['sci_name_raw']}' matched at species level"
            if level == "proxy":
                row["note"] = (row["note"] + "; " if row["note"] else "") + \
                    "stand-in species, used only with --use-proxies"
            rows.append(row)
            time.sleep(0.2)

    df = pd.concat([pd.DataFrame(rows), manual], ignore_index=True)
    df["gbif_key"] = pd.to_numeric(df["gbif_key"], errors="coerce").astype("Int64")
    if counts:
        uniq = df.loc[df["gbif_key"].notna(), "gbif_key"].unique()
        print(f"Counting Finnish records for {len(uniq)} keys ...", flush=True)
        n, fails = {}, 0
        for k in uniq:
            n[k] = fi_count(s, int(k))
            fails = fails + 1 if n[k] is None else 0
            if fails >= 3:                     # occurrence search is down; don't wait it out
                print("  GBIF occurrence search is not responding; skipping counts.")
                break
        df["n_fi_records"] = df["gbif_key"].map(n).astype("Int64")
        missing = len(uniq) - sum(v is not None for v in n.values())
        if missing:
            print(f"  GBIF occurrence search did not answer for {missing} keys "
                  "(n_fi_records left blank; rerun 'match' later to fill them).")
    df = df.sort_values(["taxon_id", "match_level"]).reset_index(drop=True)
    df.to_csv(path, index=False)

    show = ["taxon_id", "name", "match_level", "gbif_status", "gbif_accepted",
            "match_on", "match_value", "check"]
    if "n_fi_records" in df:
        show.insert(-1, "n_fi_records")
    print(df[show].to_string(index=False))
    c = df["check"].value_counts().to_dict()
    print(f"\n{path}\nchecks: {c}")
    real = df[df["check"].isin(USABLE) & (df["match_level"] != "proxy")]
    shared = real.groupby(["match_on", "match_value"])["taxon_id"].unique()
    for (on, v), ids in shared.items():
        if len(ids) > 1:
            print(f"  note: taxa {[int(i) for i in ids]} share {on}={v}; "
                  "records are counted for each.")
    if c.get("REVIEW"):
        print("Fix REVIEW rows (edit gbif_key, set check=manual) before 'download'.")


# ---------------------------------------------------------------- stage 2
def find_recent_download(user, pwd, predicate, minutes=30):
    import requests

    r = requests.get(f"{API}/occurrence/download/user/{user}", params={"limit": 10},
                     auth=(user, pwd), timeout=120)
    r.raise_for_status()
    cutoff = pd.Timestamp.now(tz="UTC") - pd.Timedelta(minutes=minutes)
    want = {str(v) for p in predicate["predicates"] if p["key"] == "TAXON_KEY" for v in p["values"]}
    for d in r.json().get("results", []):
        got = json.dumps(d.get("request", {}).get("predicate"))
        # GBIF stores a normalised copy of the predicate, so compare loosely.
        if (pd.Timestamp(d["created"]) >= cutoff and '"COUNTRY"' in got
                and all(f'"{k}"' in got for k in want)):
            return d["key"]
    return None


def stage_download(out):
    import requests
    from dotenv import load_dotenv

    env_file = DATA.parent / ".env"
    load_dotenv(env_file)  # variables already set in the shell take precedence
    user, pwd, email = (os.environ.get(k) for k in ("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL"))
    if not (user and pwd and email):
        sys.exit(f"Set GBIF_USER, GBIF_PWD and GBIF_EMAIL in {env_file} "
                 "or the environment (free account at gbif.org).")

    meta_path = out / META_FILE
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}

    if "key" not in meta:
        m = pd.read_csv(out / MATCH_FILE)
        if (m["check"] == "REVIEW").any():
            sys.exit(f"Resolve REVIEW rows in {MATCH_FILE} first.")
        m = m[m["check"].isin(USABLE)]
        keys = sorted(m["gbif_key"].dropna().astype(int).astype(str).unique().tolist())
        if not keys:
            sys.exit(f"No usable keys in {MATCH_FILE}.")
        body = {
            "creator": user,
            "notificationAddresses": [email],
            "sendNotification": True,
            "format": "SIMPLE_CSV",
            "predicate": {"type": "and", "predicates": [
                {"type": "in", "key": "TAXON_KEY", "values": keys},
                {"type": "in", "key": "COUNTRY", "values": ["FI", "AX"]},
                {"type": "equals", "key": "HAS_COORDINATE", "value": "true"},
                {"type": "equals", "key": "HAS_GEOSPATIAL_ISSUE", "value": "false"},
                {"type": "equals", "key": "OCCURRENCE_STATUS", "value": "PRESENT"},
            ]},
        }
        key = None
        for attempt in range(6):
            try:
                r = requests.post(f"{API}/occurrence/download/request", json=body,
                                  auth=(user, pwd), timeout=300)
                if r.status_code < 500:
                    r.raise_for_status()
                    key = r.text.strip()
                    break
                problem = f"HTTP {r.status_code}"
            except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
                problem = type(e).__name__
            # GBIF may have accepted it anyway; adopt it rather than request a duplicate.
            time.sleep(20)
            try:
                key = find_recent_download(user, pwd, body["predicate"])
            except Exception:
                key = None
            if key:
                break
            wait = 60 * (attempt + 1)
            print(f"{time.strftime('%H:%M:%S')}  request failed ({problem}); "
                  f"retrying in {wait} s", flush=True)
            time.sleep(wait)
        if not key:
            sys.exit("GBIF kept failing the download request (likely an outage; see "
                     "https://status.gbif.org). Rerun 'download' later.")
        meta = {"key": key, "n_taxon_keys": len(keys), "taxon_keys": keys}
        meta_path.write_text(json.dumps(meta, indent=2))
        print(f"Requested GBIF download {meta['key']} for {len(keys)} taxon keys.")

    while True:
        try:
            s = requests.get(f"{API}/occurrence/download/{meta['key']}", timeout=60).json()
        except Exception as e:                       # transient API hiccup: keep polling
            print(f"{time.strftime('%H:%M:%S')}  status check failed ({e}); retrying")
            time.sleep(60)
            continue
        status = s.get("status")
        print(f"{time.strftime('%H:%M:%S')}  {status}", flush=True)
        if status == "SUCCEEDED":
            break
        if status in ("KILLED", "CANCELLED", "FAILED"):
            meta_path.unlink()
            sys.exit(f"Download {status}. Rerun 'download' to request a new one.")
        time.sleep(60)

    meta.update(doi=s.get("doi"), total_records=s.get("totalRecords"),
                created=s.get("created"), link=s.get("downloadLink"))
    meta_path.write_text(json.dumps(meta, indent=2))
    zip_path = out / ZIP_FILE
    tmp = zip_path.with_suffix(".part")
    with requests.get(meta["link"], stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for block in r.iter_content(1 << 20):
                f.write(block)
    tmp.replace(zip_path)
    print(f"Saved {zip_path} ({meta['total_records']} records). "
          f"Cite: GBIF.org ({meta['created'][:10]}) GBIF Occurrence Download "
          f"https://doi.org/{meta['doi']}")


# ---------------------------------------------------------------- stage 3
PT_KEYS = ["taxonKey", "speciesKey", "genus", "lat", "lon", "year"]


def read_points(zip_path, max_unc):
    """Stream the SIMPLE_CSV; collapse to (taxon keys, genus, ~10 m point, year) with counts."""
    cols = ["datasetKey", "basisOfRecord", "countryCode", "taxonKey", "speciesKey", "genus",
            "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters", "year"]
    stats, basis, datasets, countries = Counter(), Counter(), Counter(), Counter()
    parts = []
    with zipfile.ZipFile(zip_path) as z:
        name = next(n for n in z.namelist() if n.endswith(".csv"))
        with z.open(name) as fh:
            for ch in pd.read_csv(fh, sep="\t", usecols=cols, dtype=str,
                                  quoting=csv.QUOTE_NONE, on_bad_lines="skip",
                                  chunksize=250_000):
                stats["rows"] += len(ch)
                keep = ~ch["basisOfRecord"].isin(["FOSSIL_SPECIMEN", "LIVING_SPECIMEN"])
                stats["dropped_basis"] += int((~keep).sum())
                ch = ch[keep]
                unc = pd.to_numeric(ch["coordinateUncertaintyInMeters"], errors="coerce")
                keep = unc.isna() | (unc <= max_unc)
                stats["dropped_uncertainty"] += int((~keep).sum())
                ch = ch[keep]
                lat = pd.to_numeric(ch["decimalLatitude"], errors="coerce").round(4)
                lon = pd.to_numeric(ch["decimalLongitude"], errors="coerce").round(4)
                keep = lat.notna() & lon.notna() & ((lat != 0) | (lon != 0))
                stats["dropped_coords"] += int((~keep).sum())
                ch, lat, lon = ch[keep], lat[keep], lon[keep]
                basis.update(ch["basisOfRecord"].fillna("").value_counts().to_dict())
                datasets.update(ch["datasetKey"].fillna("").value_counts().to_dict())
                countries.update(ch["countryCode"].fillna("").value_counts().to_dict())
                df = pd.DataFrame({
                    "taxonKey": pd.to_numeric(ch["taxonKey"], errors="coerce")
                                  .fillna(-1).astype("int64"),
                    "speciesKey": pd.to_numeric(ch["speciesKey"], errors="coerce")
                                    .fillna(-1).astype("int64"),
                    "genus": ch["genus"].fillna(""),
                    "lat": lat, "lon": lon,
                    "year": pd.to_numeric(ch["year"], errors="coerce").fillna(-1).astype(int),
                })
                parts.append(df.groupby(PT_KEYS).size().rename("n").reset_index())
    pts = pd.concat(parts).groupby(PT_KEYS, as_index=False)["n"].sum()
    stats["kept"] = int(pts["n"].sum())
    return pts, stats, basis, datasets, countries


def place_locations(pts, provinces, snap_m):
    """Assign each unique location a province code (within, else nearest <= snap_m)."""
    import geopandas as gpd

    locs = pts[["lat", "lon"]].drop_duplicates().reset_index(drop=True)
    locs["loc_id"] = locs.index
    g = gpd.GeoDataFrame(locs, geometry=gpd.points_from_xy(locs.lon, locs.lat),
                         crs=4326).to_crs(provinces.crs)
    prov = provinces[["code", "geometry"]]
    hit = gpd.sjoin(g, prov, how="left", predicate="within")
    hit = hit[~hit.index.duplicated()]
    code = hit["code"].copy()
    miss = code.isna()
    n_snapped = 0
    if miss.any():
        near = gpd.sjoin_nearest(g[miss], prov, how="left", max_distance=snap_m)
        near = near[~near.index.duplicated()]
        code.loc[near.index] = near["code"]
        n_snapped = int(near["code"].notna().sum())
    locs["region_code"] = code
    return locs.drop(columns="geometry", errors="ignore"), n_snapped


def wcvp_fi_status(out, taxa):
    """native / introduced / not_in_wcvp_fin / no_wcvp_range for plants; blank otherwise."""
    status = pd.Series(pd.NA, index=taxa["taxon_id"], dtype="object")
    path = out / "species_regions.csv"
    if not path.exists():
        return status, False
    sr = pd.read_csv(path)
    fin = sr[sr["area_code_l3"] == "FIN"].drop_duplicates("taxon_id").set_index("taxon_id")
    has_range = set(sr["taxon_id"])
    for t in taxa.itertuples():
        if t.group != "vascular plant":
            continue
        if t.taxon_id in fin.index:
            status[t.taxon_id] = "introduced" if fin.loc[t.taxon_id, "introduced"] == 1 else "native"
        elif t.taxon_id in has_range:
            status[t.taxon_id] = "not_in_wcvp_fin"
        else:
            status[t.taxon_id] = "no_wcvp_range"
    return status, True


CHECKLIST_PRESENT = ["occurs_expert", "occurs_records"]
CHECKLIST_TAXON_COLS = ["finbif_id", "finbif_name", "name_fi", "fi_origin",
                        "cultivation_escape", "red_list", "red_list_year", "invasive",
                        "sensitive", "protection", "regionally_threatened_zones",
                        "habitats", "finbif_obs_fi"]


def load_checklist(out):
    """finbif_taxon_provinces.csv and finbif_taxa.csv, or (None, None) if not built."""
    pv, tx = out / "finbif_taxon_provinces.csv", out / "finbif_taxa.csv"
    if not (pv.exists() and tx.exists()):
        return None, None
    # keep_default_na=False: the Red List code "NA" (not applicable) must stay a string.
    t = pd.read_csv(tx, keep_default_na=False, na_values=[""])
    for c in CHECKLIST_TAXON_COLS:
        if c not in t:
            t[c] = pd.NA
    return pd.read_csv(pv), t


def checklist_report(tr, tx):
    m = tx["finbif_name"].notna()
    agree = pd.crosstab(tr["group"], tr["source"])
    lines = [
        "FinBIF checklist (laji.fi)",
        f"  taxa matched: {int(m.sum())} of {len(tx)}; not in checklist: "
        + ", ".join(tx.loc[~m, "name"]),
        "  taxon x province rows by source (both = GBIF records and checklist presence):",
        *["    " + l for l in agree.to_string().splitlines()],
        "  plant origin in Finland: " + ", ".join(
            f"{k} {v}" for k, v in tx["fi_origin"].value_counts().items()),
        "  Red List: " + ", ".join(f"{k} {v}" for k, v in tx["red_list"].value_counts().items()),
    ]
    flag = lambda c: tx[tx[c].fillna(False).astype(str).str.lower().eq("true")]["name"].tolist()
    lines += [
        "  threatened / near threatened: " + ", ".join(
            f"{r.name} ({r.red_list})" for r in tx[tx["red_list"].isin(
                ["NT", "VU", "EN", "CR", "RE"])].itertuples()),
        "  invasive: " + ", ".join(flag("invasive")),
        "  legally protected: " + (", ".join(
            f"{r.name} ({r.protection})" for r in tx[tx["protection"].fillna("") != ""]
            .itertuples()) or "none"),
        "  sensitive (locations coarsened): " + (", ".join(flag("sensitive")) or "none"),
    ]
    odd = tr[~tr["checklist_category"].isna() & ~tr["in_checklist"]]
    if len(odd):
        lines += ["  checklist rows that are not presence (extinct / uncertain / old):",
                  *[f"    {r.name} {r.region_code}: {r.checklist_status} "
                    f"(GBIF records: {r.n_records})" for r in odd.itertuples()]]
    only = tr[tr["source"] == "checklist"].groupby("name")["region_code"].apply(
        lambda s: ", ".join(sorted(s)))
    if len(only):
        lines += ["  provinces in the checklist with no GBIF records (after filters):",
                  *[f"    {n}: {v}" for n, v in only.items()]]
    return lines


def stage_regions(taxa_path, samples_path, provinces_path, out,
                  min_records, max_unc, snap_m, use_proxies):
    import geopandas as gpd

    taxa = pd.read_csv(taxa_path)
    taxa["name"] = taxa.apply(display_name, axis=1)
    samples = pd.read_csv(samples_path)
    provinces = gpd.read_file(provinces_path)
    matches = pd.read_csv(out / MATCH_FILE)
    if (matches["check"] == "REVIEW").any():
        sys.exit(f"Resolve REVIEW rows in {MATCH_FILE} first.")
    meta_path = out / META_FILE
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}

    pts, stats, basis, datasets, countries = read_points(out / ZIP_FILE, max_unc)
    locs, n_snapped = place_locations(pts, provinces, snap_m)
    pts = pts.merge(locs, on=["lat", "lon"])
    n_unplaced = int(pts.loc[pts["region_code"].isna(), "n"].sum())

    # Records -> taxa, by each match row's match_on column (speciesKey, taxonKey or genus).
    use = matches[matches["check"].isin(USABLE)].copy()
    if not use_proxies:
        use_rec = use[use["match_level"] != "proxy"]
    else:
        use_rec = use
    parts = []
    for on, grp in use_rec.groupby("match_on"):
        v = grp[["taxon_id", "match_value"]].copy()
        if on != "genus":
            v["match_value"] = pd.to_numeric(v["match_value"]).astype("int64")
        parts.append(pts.merge(v, left_on=on, right_on="match_value"))
    rec = pd.concat(parts, ignore_index=True)
    rec = rec.drop_duplicates(["taxon_id", "taxonKey", "speciesKey", "genus", "loc_id", "year"])
    rec = rec[rec["region_code"].notna()]
    rec["yr"] = rec["year"].where(rec["year"] > 0)
    rec["n2000"] = rec["n"].where(rec["year"] >= 2000, 0)

    tr = (rec.groupby(["taxon_id", "region_code"])
             .agg(n_records=("n", "sum"), n_locations=("loc_id", "nunique"),
                  first_year=("yr", "min"), last_year=("yr", "max"),
                  n_records_since_2000=("n2000", "sum"))
             .reset_index())

    # FinBIF checklist (scripts/finbif_checklist.py), if present: outer-join its
    # per-province status so a row exists when either source places the taxon there.
    cl_prov, cl_taxa = load_checklist(out)
    have_cl = cl_prov is not None
    if have_cl:
        tr = tr.merge(cl_prov, on=["taxon_id", "region_code"], how="outer")
    else:
        tr["checklist_status"] = pd.NA
        tr["checklist_category"] = pd.NA
    for c in ["n_records", "n_locations", "n_records_since_2000"]:
        tr[c] = tr[c].fillna(0).astype(int)
    for c in ["first_year", "last_year"]:
        tr[c] = tr[c].astype("Int64")
    tr["in_gbif"] = tr["n_records"] > 0
    tr["in_checklist"] = tr["checklist_category"].isin(CHECKLIST_PRESENT)
    tr["source"] = ["both" if g and c else "gbif" if g else "checklist" if c else "none"
                    for g, c in zip(tr["in_gbif"], tr["in_checklist"])]

    level = (use.sort_values("match_level", key=lambda s: s.map({"species": 0, "genus": 0, "proxy": 1}))
                .drop_duplicates("taxon_id").set_index("taxon_id")["match_level"])
    pnames = provinces.set_index("code")[["name_fi", "name_lat"]]
    tr = (tr.merge(taxa[["taxon_id", "name", "group"]], on="taxon_id")
            .merge(pnames, left_on="region_code", right_index=True)
            .rename(columns={"name_fi": "region_name_fi", "name_lat": "region_name_lat"}))
    tr["match_level"] = tr["taxon_id"].map(level)
    tr = tr[["taxon_id", "name", "group", "match_level", "region_code", "region_name_fi",
             "region_name_lat", "source", "n_records", "n_locations", "first_year",
             "last_year", "n_records_since_2000", "checklist_status", "checklist_category",
             "in_gbif", "in_checklist"]].sort_values(["taxon_id", "region_code"])
    tr.to_csv(out / "finland_taxon_regions.csv", index=False)
    tg = tr[tr["in_gbif"]]                           # GBIF-only view for the record stats

    # Per taxon (all 121).
    fi_status, have_wcvp = wcvp_fi_status(out, taxa)
    ge_col = f"n_provinces_ge{min_records}"
    per = (tg.groupby("taxon_id")
             .agg(n_provinces=("region_code", "nunique"),
                  **{ge_col: ("n_records", lambda s: int((s >= min_records).sum()))},
                  n_records=("n_records", "sum"),
                  n_records_since_2000=("n_records_since_2000", "sum"),
                  first_year=("first_year", "min"), last_year=("last_year", "max")))
    per["n_locations"] = rec.groupby("taxon_id")["loc_id"].nunique()
    tx = taxa[["taxon_id", "name", "common_names", "group", "rank", "n_samples"]].merge(
        per, left_on="taxon_id", right_index=True, how="left")
    for c in ["n_provinces", ge_col, "n_records", "n_records_since_2000", "n_locations"]:
        tx[c] = tx[c].fillna(0).astype(int)
    tx["first_year"] = tx["first_year"].astype("Int64")
    tx["last_year"] = tx["last_year"].astype("Int64")
    tx["wcvp_fi_status"] = tx["taxon_id"].map(fi_status)
    tx["in_finland_gbif"] = tx["n_records"] > 0
    tx["n_provinces_checklist"] = tx["taxon_id"].map(
        tr[tr["in_checklist"]].groupby("taxon_id").size()).fillna(0).astype(int)
    tx["n_provinces_any"] = tx["taxon_id"].map(
        tr[tr["source"] != "none"].groupby("taxon_id").size()).fillna(0).astype(int)
    tx["match_level"] = tx["taxon_id"].map(level)
    notes = matches.groupby("taxon_id")["note"].apply(
        lambda s: "; ".join(sorted({x for x in s.dropna() if x})))
    tx["note"] = tx["taxon_id"].map(notes).fillna("")
    excluded = matches.loc[matches["check"] == "excluded", "taxon_id"]
    tx.loc[tx["taxon_id"].isin(excluded), "match_level"] = "excluded"
    cols = ["taxon_id", "name", "common_names", "group", "rank", "n_samples", "match_level",
            "wcvp_fi_status", "in_finland_gbif", "n_provinces", ge_col, "n_records",
            "n_records_since_2000", "n_locations", "first_year", "last_year",
            "n_provinces_checklist", "n_provinces_any"]
    if have_cl:
        tx = tx.merge(cl_taxa[["taxon_id", *CHECKLIST_TAXON_COLS]], on="taxon_id", how="left")
        cols += CHECKLIST_TAXON_COLS
    tx = tx[cols + ["note"]]
    tx.to_csv(out / "finland_taxa.csv", index=False)

    # Per province (all 21). GBIF counts, plus taxa present per the checklist / either.
    a = tg.merge(taxa[["taxon_id", "n_samples"]], on="taxon_id")
    ps = (a.groupby("region_code")
           .agg(n_taxa=("taxon_id", "nunique"),
                **{f"n_taxa_ge{min_records}": ("n_records", lambda s: int((s >= min_records).sum()))},
                n_plants=("group", lambda s: int((s == "vascular plant").sum())),
                n_fungi_lichens=("group", lambda s: int(s.isin(FUNGAL).sum())),
                n_other=("group", lambda s: int((~s.isin(FUNGAL + ["vascular plant"])).sum())),
                n_samples=("n_samples", "sum"),
                n_records=("n_records", "sum")))
    ps["n_taxa_checklist"] = tr[tr["in_checklist"]].groupby("region_code").size()
    ps["n_taxa_any"] = tr[tr["source"] != "none"].groupby("region_code").size()
    ps = (provinces.drop(columns="geometry")
                   .merge(ps, left_on="code", right_index=True, how="left")
                   .rename(columns={"code": "region_code"}))
    num = [c for c in ps.columns if c.startswith("n_")]
    ps[num] = ps[num].fillna(0).astype(int)
    ps = ps.sort_values("n_taxa", ascending=False)
    ps.to_csv(out / "finland_province_summary.csv", index=False)

    # Report.
    n_samples_linked = int(samples["taxon_id"].isin(tx.loc[tx["in_finland_gbif"], "taxon_id"]).sum())
    uniq_rec = rec.drop_duplicates(["taxonKey", "speciesKey", "genus", "loc_id", "year"])
    yrs = rec.loc[rec["yr"].notna()]
    lines = []
    if meta.get("doi"):
        lines += [f"Source: GBIF.org ({str(meta.get('created', ''))[:10]}) GBIF Occurrence Download "
                  f"https://doi.org/{meta['doi']} ({meta.get('total_records')} records, "
                  f"{meta.get('n_taxon_keys')} taxon keys, Finland + Åland)", ""]
    lines += [
        "Records",
        f"  rows read: {stats['rows']}",
        f"  dropped: fossil/living specimen {stats['dropped_basis']}, "
        f"uncertainty > {max_unc:g} m {stats['dropped_uncertainty']}, "
        f"bad/zero coordinates {stats['dropped_coords']}",
        f"  kept: {stats['kept']} records at {len(locs)} unique locations "
        f"from {len(datasets)} datasets",
        "  by country code: " + ", ".join(f"{k} {v}" for k, v in countries.most_common()),
        "  by basis of record: " + ", ".join(f"{k} {v}" for k, v in basis.most_common()),
        f"  snapped to nearest province (<= {snap_m:g} m): {n_snapped} locations; "
        f"unplaced: {n_unplaced} records",
        f"  records linked to a BioColour taxon and placed: {int(uniq_rec['n'].sum())}"
        + (f" (years {int(yrs['yr'].min())}-{int(yrs['yr'].max())}; "
           f"{int(uniq_rec['n2000'].sum())} since 2000)" if len(yrs) else ""),
        f"  taxon-record links: {int(rec['n'].sum())} (more than records: genus taxa "
        "include their species, and taxa sharing a species each count it)",
        "",
        "Taxa",
        f"  BioColour taxa: {len(taxa)}; with >= 1 Finnish record: {int(tx['in_finland_gbif'].sum())}; "
        f"with >= {min_records} records in >= 1 province: {int((tx[ge_col] > 0).sum())}",
        f"  samples linked to a taxon found in Finland: {n_samples_linked} of {len(samples)}",
        f"  taxon x province links (GBIF): {len(tg)} (>= {min_records} records: "
        f"{int((tg['n_records'] >= min_records).sum())})",
        "  taxa by group (found / total): " + ", ".join(
            f"{g} {int(d['in_finland_gbif'].sum())}/{len(d)}" for g, d in tx.groupby("group")),
    ]
    if have_wcvp:
        ct = pd.crosstab(tx.loc[tx["group"] == "vascular plant", "wcvp_fi_status"],
                         tx.loc[tx["group"] == "vascular plant", "in_finland_gbif"])
        lines += ["  plants, WCVP Finland status x has GBIF Finnish records:",
                  *["    " + l for l in ct.to_string().splitlines()]]
        surprise = tx[(tx["wcvp_fi_status"] == "not_in_wcvp_fin") & tx["in_finland_gbif"]]
        if len(surprise):
            lines += ["  not in Finland per WCVP but recorded in GBIF (garden escapes/casuals?):",
                      *[f"    {r.name}: {r.n_records} records, {r.n_provinces} provinces"
                        for r in surprise.itertuples()]]
    lines += [
        "  no Finnish records:",
        *[f"    {r.name} ({r.group}, {r.n_samples} samples){' - ' + r.note if r.note else ''}"
          for r in tx[~tx["in_finland_gbif"]].itertuples()],
        "  provinces per taxon (fewest first, found taxa only):",
        *[f"    {r.n_provinces:2d} ({r.__getattribute__(ge_col):2d} with >= {min_records})  "
          f"{r.name}: {r.n_records} records"
          for r in tx[tx["in_finland_gbif"]].sort_values(["n_provinces", "n_records"]).itertuples()],
        "",
        "Provinces (taxa found in GBIF / records; checklist; either)",
        *[f"  {r.region_code:4s} {r.name_fi:18s} {r.n_taxa:3d} taxa "
          f"({r.n_plants} plants, {r.n_fungi_lichens} fungi/lichens), {r.n_records} records; "
          f"checklist {r.n_taxa_checklist}; either {r.n_taxa_any}"
          for r in ps.itertuples()],
    ]
    if have_cl:
        lines += ["", *checklist_report(tr, tx)]
    # Proxy comparison, so the Aronia decision can be made from numbers.
    for tid, proxy in PROXIES.items():
        def keys(level):
            r = use[(use["taxon_id"] == tid) & (use["match_level"] == level)]
            return pd.to_numeric(r["match_value"]).astype("int64").tolist()
        own, prx = keys("species"), keys("proxy")

        def summ(keys):
            r = pts[pts["speciesKey"].isin(keys) & pts["region_code"].notna()]
            return int(r["n"].sum()), r["region_code"].nunique()
        n1, p1 = summ(own)
        n2, p2 = summ(prx)
        lines += ["", f"Proxy check: {taxa.set_index('taxon_id').loc[tid, 'name']} own records "
                      f"{n1} in {p1} provinces; stand-in {proxy} {n2} in {p2} provinces "
                      f"({'included' if use_proxies else 'not included; use --use-proxies'})"]

    report = "\n".join(lines)
    (out / "finland_report.txt").write_text(report, encoding="utf-8")
    print(report)


# ---------------------------------------------------------------- cli
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="stage", required=True)
    a = sub.add_parser("match")
    b = sub.add_parser("download")
    c = sub.add_parser("regions")
    for s in (a, b, c):
        s.add_argument("--out", type=Path, default=DATA / "out")
    for s in (a, c):
        s.add_argument("--taxa", type=Path, default=DATA / "taxa.csv")
    a.add_argument("--no-counts", action="store_true",
                   help="skip the per-key Finnish record counts")
    c.add_argument("--samples", type=Path, default=DATA / "samples.csv")
    c.add_argument("--provinces", type=Path, default=DATA / "boundaries" / "provinces.gpkg",
                   help="from scripts/get_provinces.py")
    c.add_argument("--min-records", type=int, default=3)
    c.add_argument("--max-uncertainty-m", type=float, default=10000)
    c.add_argument("--snap-m", type=float, default=5000)
    c.add_argument("--use-proxies", action="store_true",
                   help="count stand-in species (PROXIES) for their taxon")
    args = p.parse_args()

    inputs = [getattr(args, k) for k in ("taxa", "samples", "provinces") if hasattr(args, k)]
    missing = [str(f) for f in inputs if not f.exists()]
    if missing:
        p.error("input not found: " + ", ".join(missing))

    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    if args.stage == "match":
        stage_match(args.taxa, out, counts=not args.no_counts)
    elif args.stage == "download":
        stage_download(out)
    else:
        stage_regions(args.taxa, args.samples, args.provinces, out, args.min_records,
                      args.max_uncertainty_m, args.snap_m, args.use_proxies)


if __name__ == "__main__":
    main()
