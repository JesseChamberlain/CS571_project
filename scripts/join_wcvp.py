"""
Join BioColour taxa to WCVP native/introduced ranges (WGSRPD level 3).

Streams the two big WCVP files in chunks and keeps only the rows for your
species, so memory stays small and it runs in a minute or two.

Usage:
    python join_wcvp.py taxa.csv samples.csv wcvp_names.csv wcvp_distribution.csv out/

Outputs (in out/):
    name_matches.csv     one row per BioColour taxon: how it matched WCVP, accepted name
    species_regions.csv  long table: taxon_id x area_code_l3 x introduced (the core join)
    region_summary.csv   per L3 area: # dye species, # native, # introduced, # samples
    report.txt           counts to quote in Deliverable 1 (Q1 dataset description)
"""
import csv
import re
import sys
from pathlib import Path

import pandas as pd

CHUNK = 250_000
WCVP_READ = dict(sep="|", dtype=str, quoting=csv.QUOTE_NONE,
                 encoding="utf-8", on_bad_lines="skip", chunksize=CHUNK)


def norm(name: str) -> str:
    """Lower-case, drop hybrid signs and extra spaces: 'Aronia × mitschurinii' -> 'aronia mitschurinii'."""
    if not isinstance(name, str):
        return ""
    name = name.replace("×", " ").replace(" x ", " ")
    return re.sub(r"\s+", " ", name).strip().lower()


def main(taxa_path, samples_path, names_path, dist_path, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    taxa = pd.read_csv(taxa_path)
    samples = pd.read_csv(samples_path)

    # WCVP only covers vascular plants, and genus-only names ("Salix spp.") have no range.
    taxa["join_scope"] = "matchable"
    taxa.loc[taxa["group"] != "vascular plant", "join_scope"] = "not vascular (outside WCVP)"
    taxa.loc[(taxa["join_scope"] == "matchable") & (taxa["rank"] == "genus"),
             "join_scope"] = "genus only (no range)"
    # Infraspecific names and cultivars match at species level via the binomial column.
    taxa["key"] = taxa["binomial"].map(norm)
    targets = set(taxa.loc[taxa["join_scope"] == "matchable", "key"])

    # Pass 1: names table -> rows whose name is one of our binomials (species rank only).
    hits = []
    cols = ["plant_name_id", "taxon_rank", "taxon_status", "family",
            "taxon_name", "taxon_authors", "accepted_plant_name_id"]
    for chunk in pd.read_csv(names_path, usecols=cols, **WCVP_READ):
        chunk = chunk[chunk["taxon_rank"] == "Species"]
        chunk["key"] = chunk["taxon_name"].map(norm)
        hits.append(chunk[chunk["key"].isin(targets)])
    names = pd.concat(hits, ignore_index=True)

    # Pick one WCVP record per binomial: Accepted beats Synonym beats anything else.
    rank = {"Accepted": 0, "Synonym": 1, "Orthographic": 2}
    names["pref"] = names["taxon_status"].map(rank).fillna(9)
    names = names[names["accepted_plant_name_id"].notna()]
    best = (names.sort_values(["key", "pref"])
                 .groupby("key", as_index=False).first())
    n_candidates = names.groupby("key").size().rename("n_wcvp_candidates")
    best = best.merge(n_candidates, on="key")
    best["match_type"] = best["taxon_status"].map(
        lambda s: "exact accepted" if s == "Accepted" else f"via {s.lower()}")

    # Look up accepted names for synonyms (second small pass over the names table).
    acc_ids = set(best["accepted_plant_name_id"])
    acc = []
    for chunk in pd.read_csv(names_path, usecols=["plant_name_id", "taxon_name", "family"],
                             **WCVP_READ):
        acc.append(chunk[chunk["plant_name_id"].isin(acc_ids)])
    acc = pd.concat(acc).rename(columns={"plant_name_id": "accepted_plant_name_id",
                                         "taxon_name": "accepted_name",
                                         "family": "accepted_family"})
    best = best.merge(acc, on="accepted_plant_name_id", how="left")

    matches = taxa.merge(
        best[["key", "match_type", "n_wcvp_candidates", "plant_name_id",
              "accepted_plant_name_id", "accepted_name", "accepted_family"]],
        on="key", how="left")
    matches.loc[(matches["join_scope"] == "matchable") & matches["match_type"].isna(),
                "match_type"] = "NO MATCH - fix by hand"
    matches.drop(columns="key").to_csv(out / "name_matches.csv", index=False)

    # Pass 2: distribution table -> rows for our accepted IDs.
    dist = []
    for chunk in pd.read_csv(dist_path, **WCVP_READ):
        dist.append(chunk[chunk["plant_name_id"].isin(acc_ids)])
    dist = pd.concat(dist, ignore_index=True)
    for c in ["introduced", "extinct", "location_doubtful"]:
        dist[c] = pd.to_numeric(dist[c], errors="coerce").fillna(0).astype(int)
    dist = dist[(dist["extinct"] == 0) & (dist["location_doubtful"] == 0)]
    dist = dist.rename(columns={"plant_name_id": "accepted_plant_name_id"})

    sr = (matches[matches["accepted_plant_name_id"].notna()]
          [["taxon_id", "binomial", "accepted_plant_name_id"]]
          .merge(dist[["accepted_plant_name_id", "continent", "region",
                       "area_code_l3", "area", "introduced"]],
                 on="accepted_plant_name_id"))
    sr = sr.drop_duplicates(["taxon_id", "area_code_l3"])
    sr.to_csv(out / "species_regions.csv", index=False)

    # Region summary: what the map will shade.
    n_samp = samples.groupby("taxon_id").size().rename("n_samples")
    srs = sr.merge(n_samp, on="taxon_id")
    summary = (srs.groupby(["area_code_l3", "area", "region", "continent"])
                  .agg(n_species=("taxon_id", "nunique"),
                       n_native=("introduced", lambda x: int((x == 0).sum())),
                       n_introduced=("introduced", lambda x: int((x == 1).sum())),
                       n_samples=("n_samples", "sum"))
                  .reset_index()
                  .sort_values("n_species", ascending=False))
    summary.to_csv(out / "region_summary.csv", index=False)

    # Report for Deliverable 1.
    m = matches
    matchable = m[m["join_scope"] == "matchable"]
    matched = matchable[matchable["accepted_plant_name_id"].notna()]
    lines = [
        f"BioColour taxa: {len(taxa)}  samples: {len(samples)}",
        "Join scope:",
        *[f"  {k}: {v}" for k, v in taxa["join_scope"].value_counts().items()],
        f"Matchable taxa matched to WCVP: {len(matched)}/{len(matchable)} "
        f"({len(matched)/max(len(matchable),1):.0%})",
        "Match types:",
        *[f"  {k}: {v}" for k, v in matchable["match_type"].value_counts().items()],
        f"Samples linked to a WCVP range: "
        f"{samples['taxon_id'].isin(sr['taxon_id']).sum()}/{len(samples)}",
        f"Species x region links: {len(sr)}  "
        f"(native {int((sr.introduced==0).sum())}, introduced {int((sr.introduced==1).sum())})",
        f"L3 areas with >=1 dye plant: {summary['area_code_l3'].nunique()} of 369",
        "Top 10 areas by species count:",
        *[f"  {r.area_code_l3} {r.area}: {r.n_species} spp "
          f"({r.n_native} native / {r.n_introduced} introduced)"
          for r in summary.head(10).itertuples()],
        "Unmatched (fix by hand):",
        *[f"  {r.taxon_id} {r.sci_name_raw}" for r in
          matchable[matchable["accepted_plant_name_id"].isna()].itertuples()],
        "Ambiguous (>1 WCVP record for the name, check):",
        *[f"  {r.taxon_id} {r.binomial} ({int(r.n_wcvp_candidates)} records, used {r.match_type})"
          for r in matched[matched["n_wcvp_candidates"] > 1].itertuples()],
    ]
    report = "\n".join(lines)
    (out / "report.txt").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:])
