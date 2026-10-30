"""
Fungi and lichens: GBIF name matching + occurrence records -> WGSRPD level-3 regions.

Companion to join_wcvp.py. It produces the same taxon x region shape, then merges
both into one table where every link has a status: native, introduced or recorded.
"recorded" means GBIF has at least --min-records clean records there; fungal and
lichen records almost never say whether a species is native.

Three stages, run in order (from anywhere; paths default to the project's data/
folder, override with --taxa, --samples, --wgsrpd, --out):

  1. python scripts/fungi_gbif.py match
       Reads data/taxa.csv. Matches the 32 fungus + lichen binomials against the GBIF Backbone
       (fungal names there come from Species Fungorum). No account needed.
       Writes data/out/fungi_name_matches.csv. Check rows flagged REVIEW before stage 2.

  2. python scripts/fungi_gbif.py download
       Requests one GBIF occurrence download (SIMPLE_CSV) for all matched species:
       has coordinates, no geospatial issues, status PRESENT.
       Needs a free GBIF account. Put the credentials in a .env file at the
       project root (git-ignored), or set them as environment variables:
         GBIF_USER, GBIF_PWD, GBIF_EMAIL
       Polls until GBIF finishes (often 5-30 min), then saves
       data/out/gbif_occurrences.zip and data/out/gbif_download.json (with the DOI
       you must cite). If interrupted, rerun the same command: it resumes the same
       download.

  3. python scripts/fungi_gbif.py regions
       Reads data/taxa.csv, data/samples.csv and data/wgsrpd/, a clone of
       https://github.com/tdwg/wgsrpd:
         git clone https://github.com/tdwg/wgsrpd data/wgsrpd
       Joins points to level-3 polygons, counts records per species x region,
       keeps regions with >= --min-records (default 3).
       Writes data/out/fungi_regions.csv, and if data/out/species_regions.csv from
       join_wcvp.py exists, also data/out/taxon_regions_all.csv,
       data/out/region_summary_all.csv and data/out/fungi_report.txt.

Options for stage 3: --min-records 3  --max-uncertainty-m 20000  --snap-deg 0.1
"""
import argparse
import csv
import json
import os
import sys
import time
import warnings
import zipfile
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"

API = "https://api.gbif.org/v1"
FUNGAL_GROUPS = ["fungus", "lichen"]


# ---------------------------------------------------------------- stage 1
def stage_match(taxa_path, out):
    import requests

    taxa = pd.read_csv(taxa_path)
    fl = taxa[taxa["group"].isin(FUNGAL_GROUPS)].copy()
    rows = []
    for t in fl.itertuples():
        r = requests.get(f"{API}/species/match",
                         params={"name": t.binomial, "kingdom": "Fungi"}, timeout=30)
        r.raise_for_status()
        m = r.json()
        status = m.get("status", "")
        accepted_key = m.get("acceptedUsageKey") if status == "SYNONYM" else m.get("usageKey")
        species_key = m.get("speciesKey") or accepted_key
        ok = (m.get("matchType") == "EXACT" and m.get("rank") == "SPECIES"
              and m.get("kingdom") == "Fungi")
        rows.append({
            "taxon_id": t.taxon_id, "binomial": t.binomial, "group": t.group,
            "gbif_match_type": m.get("matchType"), "gbif_confidence": m.get("confidence"),
            "gbif_rank": m.get("rank"), "gbif_status": status,
            "gbif_usage_key": m.get("usageKey"), "gbif_species_key": species_key,
            "gbif_name": m.get("scientificName"), "gbif_family": m.get("family"),
            "check": "ok" if ok else "REVIEW",
        })
        time.sleep(0.2)  # be polite to the API
    df = pd.DataFrame(rows)
    df.to_csv(out / "fungi_name_matches.csv", index=False)
    print(df[["taxon_id", "binomial", "gbif_match_type", "gbif_status",
              "gbif_name", "check"]].to_string(index=False))
    n_ok = (df["check"] == "ok").sum()
    print(f"\n{n_ok}/{len(df)} exact species matches. "
          "Fix any REVIEW rows in fungi_name_matches.csv (set gbif_species_key, "
          "check=ok) before running 'download'.")


# ---------------------------------------------------------------- stage 2
def stage_download(out):
    import requests
    from dotenv import load_dotenv

    env_file = DATA.parent / ".env"
    load_dotenv(env_file)  # variables already set in the shell take precedence
    user, pwd, email = (os.environ.get(k) for k in ("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL"))
    if not (user and pwd and email):
        sys.exit(f"Set GBIF_USER, GBIF_PWD and GBIF_EMAIL in {env_file} "
                 "or the environment (free account at gbif.org).")

    meta_path = out / "gbif_download.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}

    if "key" not in meta:
        m = pd.read_csv(out / "fungi_name_matches.csv")
        keys = (m.loc[m["check"] == "ok", "gbif_species_key"]
                 .dropna().astype(int).astype(str).unique().tolist())
        if not keys:
            sys.exit("No usable species keys in fungi_name_matches.csv.")
        body = {
            "creator": user,
            "notificationAddresses": [email],
            "sendNotification": True,
            "format": "SIMPLE_CSV",
            "predicate": {"type": "and", "predicates": [
                {"type": "in", "key": "TAXON_KEY", "values": keys},
                {"type": "equals", "key": "HAS_COORDINATE", "value": "true"},
                {"type": "equals", "key": "HAS_GEOSPATIAL_ISSUE", "value": "false"},
                {"type": "equals", "key": "OCCURRENCE_STATUS", "value": "PRESENT"},
            ]},
        }
        r = requests.post(f"{API}/occurrence/download/request", json=body,
                          auth=(user, pwd), timeout=60)
        r.raise_for_status()
        meta = {"key": r.text.strip(), "n_species_keys": len(keys)}
        meta_path.write_text(json.dumps(meta, indent=2))
        print(f"Requested GBIF download {meta['key']} for {len(keys)} species.")

    while True:
        s = requests.get(f"{API}/occurrence/download/{meta['key']}", timeout=60).json()
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
    zip_path = out / "gbif_occurrences.zip"
    with requests.get(meta["link"], stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(zip_path, "wb") as f:
            for block in r.iter_content(1 << 20):
                f.write(block)
    print(f"Saved {zip_path} ({meta['total_records']} records). "
          f"Cite: GBIF.org ({meta['created'][:10]}) GBIF Occurrence Download "
          f"https://doi.org/{meta['doi']}")


# ---------------------------------------------------------------- stage 3
def load_wgsrpd(repo):
    import geopandas as gpd

    repo = Path(repo)
    tbl = repo / "109-488-1-ED" / "2nd Edition"
    rd = dict(sep="*", encoding="latin-1", dtype=str)
    l1 = pd.read_csv(tbl / "tblLevel1.txt", **rd)
    l2 = pd.read_csv(tbl / "tblLevel2.txt", **rd)
    code = lambda s: int(float(str(s).replace(",", ".")))
    cont = {code(r["L1 code"]): r["L1 continent"].upper() for _, r in l1.iterrows()}
    reg = {code(r["L2 code"]): r["L2 region"] for _, r in l2.iterrows()}

    g = gpd.read_file(repo / "level3" / "level3.shp")
    g["geometry"] = g.geometry.make_valid()
    g = g.rename(columns={"LEVEL3_COD": "area_code_l3", "LEVEL3_NAM": "area"})
    g["region"] = g["LEVEL2_COD"].astype(int).map(reg)
    g["continent"] = g["LEVEL1_COD"].astype(int).map(cont)
    return g[["area_code_l3", "area", "region", "continent", "geometry"]]


def read_points(zip_path, max_unc):
    """Stream the SIMPLE_CSV and collapse to unique (species, ~10 m) points with counts."""
    cols = ["speciesKey", "decimalLatitude", "decimalLongitude",
            "coordinateUncertaintyInMeters", "basisOfRecord"]
    parts = []
    with zipfile.ZipFile(zip_path) as z:
        name = next(n for n in z.namelist() if n.endswith(".csv"))
        with z.open(name) as fh:
            for ch in pd.read_csv(fh, sep="\t", usecols=cols, dtype=str,
                                  quoting=csv.QUOTE_NONE, on_bad_lines="skip",
                                  chunksize=250_000):
                ch = ch[ch["speciesKey"].notna()]
                ch = ch[~ch["basisOfRecord"].isin(["FOSSIL_SPECIMEN", "LIVING_SPECIMEN"])]
                unc = pd.to_numeric(ch["coordinateUncertaintyInMeters"], errors="coerce")
                ch = ch[unc.isna() | (unc <= max_unc)]
                lat = pd.to_numeric(ch["decimalLatitude"], errors="coerce").round(4)
                lon = pd.to_numeric(ch["decimalLongitude"], errors="coerce").round(4)
                ch = pd.DataFrame({"speciesKey": ch["speciesKey"].astype(float).astype(int),
                                   "lat": lat, "lon": lon}).dropna()
                ch = ch[(ch.lat != 0) | (ch.lon != 0)]          # null-island records
                parts.append(ch.groupby(["speciesKey", "lat", "lon"]).size()
                               .rename("n").reset_index())
    pts = pd.concat(parts).groupby(["speciesKey", "lat", "lon"], as_index=False)["n"].sum()
    return pts


def stage_regions(taxa_path, samples_path, repo, out, min_records, max_unc, snap_deg):
    import geopandas as gpd

    taxa = pd.read_csv(taxa_path)
    samples = pd.read_csv(samples_path)
    matches = pd.read_csv(out / "fungi_name_matches.csv")
    matches = matches[(matches["check"] == "ok") & matches["gbif_species_key"].notna()]
    matches["gbif_species_key"] = matches["gbif_species_key"].astype(int)

    regions = load_wgsrpd(repo)
    pts = read_points(out / "gbif_occurrences.zip", max_unc)
    gpts = gpd.GeoDataFrame(pts, geometry=gpd.points_from_xy(pts.lon, pts.lat), crs=4326)

    hit = gpd.sjoin(gpts, regions, how="left", predicate="within")
    hit = hit[~hit.index.duplicated()]                     # border points: keep one
    miss = hit["area_code_l3"].isna()
    # Coastal records (shore lichens, sea-cliff fungi) often fall just offshore:
    # snap them to the nearest region within snap_deg degrees.
    if miss.any():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")               # distance in degrees is fine here
            near = gpd.sjoin_nearest(gpts[miss], regions, how="left",
                                     max_distance=snap_deg)
        near = near[~near.index.duplicated()]
        for c in ["area_code_l3", "area", "region", "continent"]:
            hit.loc[near.index, c] = near[c]
    n_unplaced = int(hit.loc[hit["area_code_l3"].isna(), "n"].sum())
    hit = hit.dropna(subset=["area_code_l3"])

    counts = (hit.groupby(["speciesKey", "area_code_l3", "area", "region", "continent"])
                 ["n"].sum().rename("n_records").reset_index())
    kept = counts[counts["n_records"] >= min_records]

    fr = (matches[["taxon_id", "binomial", "group", "gbif_species_key"]]
          .merge(kept, left_on="gbif_species_key", right_on="speciesKey"))
    fr["status"] = "recorded"
    fr["introduced"] = pd.NA
    cols = ["taxon_id", "binomial", "group", "continent", "region",
            "area_code_l3", "area", "status", "introduced", "n_records"]
    fr = fr[cols].sort_values(["taxon_id", "area_code_l3"])
    fr.to_csv(out / "fungi_regions.csv", index=False)

    lines = [
        f"Fungus + lichen taxa matched: {len(matches)}",
        f"GBIF points after filters: {int(pts['n'].sum())} records "
        f"({len(pts)} unique locations); unplaced (open sea etc.): {n_unplaced}",
        f"Species x region links before threshold: {len(counts)}; "
        f">= {min_records} records: {len(fr)}",
        f"Regions with >= 1 fungus or lichen: {fr['area_code_l3'].nunique()} of 369",
        "Regions per taxon (fewest first):",
        *[f"  {b}: {n}" for b, n in fr.groupby("binomial").size().sort_values().items()],
        "Matched but no region passes the threshold:",
        *[f"  {b}" for b in sorted(set(matches["binomial"]) - set(fr["binomial"]))],
    ]

    # Merge with the WCVP plant output, if present.
    plant_path = out / "species_regions.csv"
    if plant_path.exists():
        pl = pd.read_csv(plant_path)
        pl = pl.merge(taxa[["taxon_id", "group"]], on="taxon_id")
        pl["status"] = pl["introduced"].map({0: "native", 1: "introduced"})
        pl["n_records"] = pd.NA
        # Use the same region names for both sources.
        pl = pl.drop(columns=["continent", "region", "area"]).merge(
            regions.drop(columns="geometry"), on="area_code_l3", how="left")
        allr = pd.concat([pl[cols], fr], ignore_index=True)
        allr.to_csv(out / "taxon_regions_all.csv", index=False)

        n_samp = samples.groupby("taxon_id").size().rename("n_samples")
        a = allr.merge(n_samp, on="taxon_id")
        summ = (a.groupby(["area_code_l3", "area", "region", "continent"])
                 .agg(n_taxa=("taxon_id", "nunique"),
                      n_native=("status", lambda s: int((s == "native").sum())),
                      n_introduced=("status", lambda s: int((s == "introduced").sum())),
                      n_recorded_fungi_lichens=("status", lambda s: int((s == "recorded").sum())),
                      n_samples=("n_samples", "sum"))
                 .reset_index().sort_values("n_taxa", ascending=False))
        summ.to_csv(out / "region_summary_all.csv", index=False)
        linked = samples["taxon_id"].isin(allr["taxon_id"]).sum()
        lines += [
            "",
            "Combined with WCVP plants:",
            f"  taxa with a range: {allr['taxon_id'].nunique()} of {len(taxa)}",
            f"  samples linked to a range: {linked} of {len(samples)}",
            "  links by status: " + ", ".join(
                f"{k} {v}" for k, v in allr["status"].value_counts().items()),
            f"  regions with >= 1 dye taxon: {summ['area_code_l3'].nunique()} of 369",
            "  top 10 regions:",
            *[f"    {r.area_code_l3} {r.area}: {r.n_taxa} taxa ({r.n_native} native, "
              f"{r.n_introduced} introduced, {r.n_recorded_fungi_lichens} fungi/lichens)"
              for r in summ.head(10).itertuples()],
        ]
    else:
        lines.append("\n(out/species_regions.csv not found: run join_wcvp.py first "
                     "to get the combined table.)")

    report = "\n".join(lines)
    (out / "fungi_report.txt").write_text(report, encoding="utf-8")
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
    c.add_argument("--samples", type=Path, default=DATA / "samples.csv")
    c.add_argument("--wgsrpd", type=Path, default=DATA / "wgsrpd",
                   help="clone of https://github.com/tdwg/wgsrpd")
    c.add_argument("--min-records", type=int, default=3)
    c.add_argument("--max-uncertainty-m", type=float, default=20000)
    c.add_argument("--snap-deg", type=float, default=0.1)
    args = p.parse_args()

    inputs = [getattr(args, k) for k in ("taxa", "samples", "wgsrpd") if hasattr(args, k)]
    missing = [str(f) for f in inputs if not f.exists()]
    if missing:
        p.error("input not found: " + ", ".join(missing))

    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    if args.stage == "match":
        stage_match(args.taxa, out)
    elif args.stage == "download":
        stage_download(out)
    else:
        stage_regions(args.taxa, args.samples, args.wgsrpd, out,
                      args.min_records, args.max_uncertainty_m, args.snap_deg)


if __name__ == "__main__":
    main()
