"""
Join BioColour taxa to the FinBIF / laji.fi Finnish species checklist.

Inputs: TSV downloads from https://laji.fi/taxon/list (Finnish taxa only, English UI,
"Download the list" -> Text file), one per organism group, in data/finbif/:
    vascular_plants.tsv   informalTaxonGroups=MVL.343
    fungi_lichens.tsv     informalTaxonGroups=MVL.233
Columns used: Identifier, Scientific name, Taxon rank, Synonyms, Misapplied Names,
Finnish name, English name, Status in Finland, Red List category in Finland,
Regulatory Status, Sensitive species, Invasive species, Habitats,
Occurences by biogeographical province (sic), Observation count from Finland.

What the checklist adds per taxon: Finnish name, origin in Finland (native, archaeophyte,
alien established / casual; plants only), Red List category, invasive / sensitive /
legally protected flags, habitats, and a status for each biogeographical province.
For vascular plants every province status is "Occurs on the basis of occurrence data -
not an expert evaluation" (derived from FinBIF records, so not independent of GBIF).
For fungi and lichens many are expert statuses ("Occurs", "Extinct", "Rare", ...).

Matching, in order, within the taxon's own file: scientific name, FinBIF synonyms, then
the GBIF accepted name (from data/out/finland_name_matches.csv, if present) and its
synonyms, then misapplied names (always flagged REVIEW).

Note: red_list uses the official codes, including "NA" (not applicable, e.g. aliens).
Read finbif_taxa.csv with keep_default_na=False so pandas doesn't turn it into NaN. Genus-only taxa ("Salix spp.") take the
union of all checklist rows in that genus. Fix a match by adding it to MANUAL below.

Usage (from anywhere; paths default to the project's data/ folder):
    python scripts/finbif_checklist.py

Outputs (data/out/), merged by `finland.py regions` when present:
    finbif_matches.csv           one row per BioColour taxon: how it matched the checklist
    finbif_taxa.csv              checklist attributes per matched taxon
    finbif_taxon_provinces.csv   taxon x province checklist status
"""
import argparse
import re
import unicodedata
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"

FILES = {"vascular plant": "vascular_plants.tsv",
         "fungus": "fungi_lichens.tsv", "lichen": "fungi_lichens.tsv"}
PROVINCE_COL = "Occurences by biogeographical province"   # spelling as exported

# taxon_id -> FinBIF Identifier, for matches fixed by hand after review.
MANUAL = {}

# Province status -> category; lower rank wins when a genus combines several species.
CATEGORY = {
    "Occurs": "occurs_expert", "Stable population": "occurs_expert",
    "Common": "occurs_expert", "Rare": "occurs_expert", "Very rare": "occurs_expert",
    "Occurs on the basis of occurrence data – not an expert evaluation": "occurs_records",
    "Old records only": "old_records",
    "Uncertain records only": "uncertain", "undefined": "uncertain",
    "Extinct": "extinct",
}
CATEGORY_RANK = {"occurs_expert": 0, "occurs_records": 1, "old_records": 2,
                 "uncertain": 3, "extinct": 4, "other": 5}
PRESENT = ["occurs_expert", "occurs_records"]

RANK_MARKERS = {"subsp.", "var.", "f.", "×", "x", "nothosubsp.", "nothovar.", "agg.",
                "sect.", "ssp."}
PROTECTED_NATIONAL = re.compile(r"protected in the whole country|strict protection", re.I)
PROTECTED_REGIONAL = re.compile(r"protected in the regions? of (.+?) \(", re.I)
REGIONAL_THREAT = re.compile(r"Regionally Threatened \d{4} - (.+)")


def protection(regs):
    """'national', 'regional: <regions>' or ''. Appendix 6 'Threatened species' is a
    listing, not a picking ban, so it is not counted here."""
    if any(PROTECTED_NATIONAL.search(x) for x in regs):
        return "national"
    regions = [m.group(1) for x in regs if (m := PROTECTED_REGIONAL.search(x))]
    return "regional: " + "; ".join(regions) if regions else ""


def norm(s):
    s = unicodedata.normalize("NFKD", str(s))
    s = "".join(c for c in s if not unicodedata.combining(c))
    return re.sub(r"\s+", " ", s.replace("×", " ").lower()).strip()


def canonical(name):
    """Drop authors: 'Lycopodium alpinum L.' -> 'lycopodium alpinum'."""
    toks = str(name).replace("×", " × ").split()
    out = toks[:1]
    for t in toks[1:]:
        if t in RANK_MARKERS or (t[:1].islower() and not t.endswith(".")):
            out.append(t)
        else:
            break
    return norm(" ".join(out))


def split_list(s):
    return [p.strip() for p in str(s).split(";") if p.strip()]


def parse_provinces(s):
    out = {}
    for part in split_list(s):
        if ":" in part:
            code, status = part.split(":", 1)
            out[code.strip()] = status.strip()
    return out


def plant_origin(status):
    """'TNV - alien, new, resident; 3 - spreading ...' -> ('alien, established', '3', False)."""
    codes = re.findall(r"\b([AT])([ON])([VS])([XU]?)\b", status)
    cult = re.findall(r"\b([123]) - spreading in the wild", status)
    if not codes:
        return ("not evaluated" if status.startswith("I -") else ""), ";".join(cult), False
    a, age, est, gone = codes[0]
    if a == "A":
        origin = "native"
    elif age == "O":
        origin = "archaeophyte"
    else:
        origin = "alien, established" if est == "V" else "alien, casual"
    return origin, ";".join(cult), bool(gone)


def load_checklists(folder):
    frames = {}
    for fname in sorted(set(FILES.values())):
        d = pd.read_csv(folder / fname, sep="\t", dtype=str, keep_default_na=False)
        d["canon"] = d["Scientific name"].map(norm)
        d["genus"] = d["canon"].str.split().str[0]
        d["file"] = fname
        frames[fname] = d
    return frames


def build_index(d, column, many):
    idx = {}
    for i, v in zip(d.index, d[column]):
        for name in (split_list(v) if many else [v]):
            idx.setdefault(canonical(name) if many else norm(name), []).append(i)
    return idx


def pick(d, rows):
    """Prefer species rank, then the row with the most Finnish observations."""
    sub = d.loc[rows].copy()
    sub["_sp"] = (sub["Taxon rank"] != "species").astype(int)
    sub["_obs"] = -pd.to_numeric(sub["Observation count from Finland"], errors="coerce").fillna(0)
    return sub.sort_values(["_sp", "_obs"]).index[0], len(set(rows))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--taxa", type=Path, default=DATA / "taxa.csv")
    p.add_argument("--finbif", type=Path, default=DATA / "finbif")
    p.add_argument("--out", type=Path, default=DATA / "out")
    args = p.parse_args()

    taxa = pd.read_csv(args.taxa)
    frames = load_checklists(args.finbif)
    gm = args.out / "finland_name_matches.csv"
    accepted = {}
    if gm.exists():
        m = pd.read_csv(gm)
        m = m[m["match_level"].isin(["species", "genus"])]
        accepted = dict(zip(m["taxon_id"], m["gbif_accepted"]))
    idx = {f: {"sci": build_index(d, "canon", False),
               "syn": build_index(d, "Synonyms", True),
               "mis": build_index(d, "Misapplied Names", True)} for f, d in frames.items()}

    matches, attrs, provs = [], [], []
    for t in taxa.itertuples():
        name = t.binomial if isinstance(t.binomial, str) else t.sci_name_raw
        base = {"taxon_id": t.taxon_id, "name": name, "group": t.group}
        fname = FILES.get(t.group)
        if fname is None:
            matches.append({**base, "finbif_match": "none", "check": "none",
                            "note": f"no checklist downloaded for group '{t.group}'"})
            continue
        d, ix = frames[fname], idx[fname]

        if t.rank == "genus":
            rows = d.index[(d["genus"] == norm(t.genus))].tolist()
            if not rows:
                matches.append({**base, "finbif_match": "none", "check": "none",
                                "note": "genus not in checklist"})
                continue
            matches.append({**base, "finbif_match": "genus", "finbif_name": t.genus,
                            "finbif_rows": len(rows), "check": "ok",
                            "note": f"union of {len(rows)} checklist rows in the genus"})
            best = {}
            for r in rows:
                for code, status in parse_provinces(d.at[r, PROVINCE_COL]).items():
                    cat = CATEGORY.get(status, "other")
                    if code not in best or CATEGORY_RANK[cat] < CATEGORY_RANK[best[code][1]]:
                        best[code] = (status, cat)
            provs += [{"taxon_id": t.taxon_id, "region_code": c, "checklist_status": s,
                       "checklist_category": k} for c, (s, k) in best.items()]
            attrs.append({"taxon_id": t.taxon_id, "finbif_match": "genus",
                          "finbif_name": t.genus,
                          "invasive": bool((d.loc[rows, "Invasive species"] == "Yes").any()),
                          "finbif_obs_fi": int(pd.to_numeric(
                              d.loc[rows, "Observation count from Finland"],
                              errors="coerce").fillna(0).sum())})
            continue

        how, hit, n_hits = None, None, 0
        if t.taxon_id in MANUAL:
            rows = d.index[d["Identifier"] == MANUAL[t.taxon_id]].tolist()
            if rows:
                how, (hit, n_hits) = "manual", pick(d, rows)
        # BioColour's own name first (more specific than GBIF's, which may lump it).
        tries = [("scientific", ix["sci"], norm(t.binomial)),
                 ("synonym", ix["syn"], canonical(t.binomial)),
                 ("gbif_accepted", ix["sci"], norm(accepted.get(t.taxon_id, "")) or None),
                 ("synonym", ix["syn"], canonical(accepted.get(t.taxon_id, "")) or None),
                 ("misapplied", ix["mis"], canonical(t.binomial))]
        for label, index, key in tries:
            if hit is not None:
                break
            if key and key in index:
                how, (hit, n_hits) = label, pick(d, index[key])
        if hit is None:
            matches.append({**base, "finbif_match": "none", "check": "none",
                            "note": "not in the Finnish checklist download"})
            continue
        r = d.loc[hit]
        check = "REVIEW" if how == "misapplied" or n_hits > 1 else "ok"
        note = f"{n_hits} candidate rows; picked the species-rank one" if n_hits > 1 else ""
        matches.append({**base, "finbif_match": how, "finbif_id": r["Identifier"],
                        "finbif_name": r["Scientific name"], "finbif_rank": r["Taxon rank"],
                        "finbif_rows": 1, "check": check, "note": note})

        status = r["Status in Finland"]
        origin, cult, old_only = plant_origin(status) if t.group == "vascular plant" else ("", "", False)
        rl = re.match(r"([A-Z]{2})(?: \((\d{4})\))?", r["Red List category in Finland"])
        regs = [x for x in split_list(r["Regulatory Status"]) if not x.startswith("OUTDATED")]
        attrs.append({
            "taxon_id": t.taxon_id, "finbif_match": how, "finbif_id": r["Identifier"],
            "finbif_name": r["Scientific name"], "name_fi": r["Finnish name"],
            "name_en_finbif": r["English name"], "fi_status_finbif": status,
            "fi_origin": origin, "cultivation_escape": cult, "only_old_records": old_only,
            "red_list": rl.group(1) if rl else "", "red_list_year": rl.group(2) if rl else "",
            "invasive": r["Invasive species"] == "Yes",
            "sensitive": r["Sensitive species"].upper() == "TRUE",
            "protection": protection(regs),
            "regionally_threatened_zones": "; ".join(
                m.group(1) for x in regs if (m := REGIONAL_THREAT.search(x))),
            "regulatory_status": "; ".join(regs), "habitats": r["Habitats"],
            "finbif_obs_fi": pd.to_numeric(r["Observation count from Finland"], errors="coerce"),
        })
        for code, st in parse_provinces(r[PROVINCE_COL]).items():
            provs.append({"taxon_id": t.taxon_id, "region_code": code, "checklist_status": st,
                          "checklist_category": CATEGORY.get(st, "other")})

    mt = pd.DataFrame(matches)
    at = pd.DataFrame(attrs)
    pv = pd.DataFrame(provs)
    n_prov = pv[pv["checklist_category"].isin(PRESENT)].groupby("taxon_id").size()
    at["n_provinces_checklist"] = at["taxon_id"].map(n_prov).fillna(0).astype(int)
    args.out.mkdir(parents=True, exist_ok=True)
    mt.to_csv(args.out / "finbif_matches.csv", index=False)
    at.to_csv(args.out / "finbif_taxa.csv", index=False)
    pv.sort_values(["taxon_id", "region_code"]).to_csv(
        args.out / "finbif_taxon_provinces.csv", index=False)

    print(mt["finbif_match"].value_counts().to_string())
    show = mt[mt["finbif_match"] != "scientific"]
    print("\nNot an exact scientific-name match:")
    print(show[["taxon_id", "name", "group", "finbif_match", "finbif_name", "check", "note"]]
          .to_string(index=False))
    print(f"\n{(mt['check'] == 'REVIEW').sum()} REVIEW; "
          f"{len(pv)} taxon x province rows; outputs in {args.out}")


if __name__ == "__main__":
    main()
