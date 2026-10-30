"""
Parse a saved copy of the BioColour swatch page
(https://library-biocolour.2.rahtiapp.fi/library/color/swatches/)
into tidy tables for the Natural Dye Atlas.

Usage (from anywhere; paths default to the project's data/ folder):
    python scripts/parse_swatches.py
    python scripts/parse_swatches.py --html path/to/pull.html --out other_dir/

Defaults: reads data/pull.html, writes into data/.

Writes (in --out):
    samples.csv  one row per color measurement (664)
    taxa.csv     one row per scientific name, with organism group and name parts
"""
import argparse
import html
import re
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"

# ---------------------------------------------------------------- lookups

# Genera in the library that are not vascular plants (WCVP will not cover them).
NON_PLANT_GENERA = {
    "fungus": {
        "Cortinarius", "Pisolithus", "Phaeolus", "Mutinus", "Gomphidius",
        "Lactarius", "Sarcodon", "Suillus", "Hapalopilus", "Inonotus",
        "Paxillus", "Russula", "Phellinus", "Imleria", "Tapinella",
        "Xerocomus", "Hypomyces", "Fomes", "Hydnellum", "Chalciporus",
    },
    "lichen": {
        "Pseudevernia", "Platismatia", "Cetraria", "Parmelia",
        "Umbilicaria", "Hypogymnia", "Evernia",
    },
    "alga": {"Fucus"},
    "insect": {"Dactylopius"},
}
GENUS_TO_GROUP = {g: grp for grp, gs in NON_PLANT_GENERA.items() for g in gs}

# Plant-part keywords found in the sample labels, mapped to a small vocabulary.
PART_KEYWORDS = [
    ("whole plant", "whole plant"), ("leaves", "leaves"), ("leaf", "leaves"),
    ("needles", "leaves"), ("tops", "aerial parts"), ("top parts", "aerial parts"),
    ("flowers", "flowers"), ("flower", "flowers"), ("fllowers", "flowers"),
    ("inflorescence", "flowers"), ("panicles", "flowers"), ("catkins", "flowers"),
    ("berries", "fruit"), ("fruits", "fruit"), ("rose hip", "fruit"),
    ("seedpods", "seeds"), ("seeds", "seeds"), ("seed", "seeds"),
    ("cones", "cones"), ("bark", "bark"), ("wood", "wood"),
    ("branches", "stems/twigs"), ("branch", "stems/twigs"), ("twigs", "stems/twigs"),
    ("stems", "stems/twigs"), ("stem", "stems/twigs"), ("stipes", "stems/twigs"),
    ("roots", "roots"), ("root", "roots"),
    ("peels", "peel/skin"), ("peel", "peel/skin"), ("skins", "peel/skin"),
    ("skin", "peel/skin"), ("caps", "caps"), ("cores", "cores"),
    ("spruce tips", "shoots"),
]


def mordant_class(raw: str) -> str:
    """Collapse the free-text mordant field into one primary class."""
    if not raw:
        return "not recorded"
    s = raw.lower()
    if s in {"none", "no mordanting"}:
        return "none"
    if "vat" in s:
        return "vat (reduction)"
    metals = []
    if "alum" in s or "aluna" in s or "(al)" in s:
        metals.append("alum")
    if "feso" in s or "iron" in s or "rust" in s:
        metals.append("iron")
    if "cuso" in s:
        metals.append("copper")
    if "sncl" in s:
        metals.append("tin")
    if len(metals) > 1:
        return "+".join(metals)
    if metals:
        return metals[0]
    if "tannin" in s or "willow" in s:
        return "tannin"
    if "oxalic" in s:
        return "oxalic acid"
    if "cream of tartar" in s:
        return "cream of tartar only"
    if "salt" in s:
        return "salt"
    if "ammo" in s:
        return "ammonia"
    return "other"


def hue_family(C: float, h: float) -> str:
    """Rough named hue bands on the CIELAB hue angle; low chroma counts as neutral."""
    if C < 10:
        return "neutral"
    if h < 40 or h >= 345:
        return "red/pink"
    if h < 70:
        return "orange/brown"
    if h < 100:
        return "yellow"
    if h < 170:
        return "green"
    if h < 290:
        return "blue"
    return "purple"


def parse_label(title: str) -> dict:
    """'Chokeberry (fresh berries 10:1) (5)' -> common name, prep text, replicate."""
    m = re.match(r"^(.*?)\s*\((\d+)\)\s*$", title)
    label, rep = (m.group(1), int(m.group(2))) if m else (title, None)
    first_paren = label.find("(")
    common = label[:first_paren].strip() if first_paren >= 0 else label.strip()
    prep = " ".join(re.findall(r"\(([^)]*)\)", label))
    text = prep.lower() if prep else label.lower()

    parts = []
    for kw, cat in PART_KEYWORDS:
        if re.search(rf"\b{re.escape(kw)}\b", text) and cat not in parts:
            parts.append(cat)
    if "rose hip" in label.lower() and "fruit" not in parts:
        parts.append("fruit")

    state = None
    if re.search(r"\bhalf[- ]dry\b", text):
        state = "half-dry"
    elif re.search(r"\b(dried|dry)\b", text):
        state = "dried"
    elif re.search(r"\bfresh\b", text):
        state = "fresh"

    ratio = re.search(r"\b(\d+:\d+)\b", text)
    mixed = bool(re.search(r"\+| and (?:\w+ )?\(", label)) or "afterdye" in label.lower()
    return {
        "label": label,
        "replicate": rep,
        "common_name": common,
        "prep_raw": prep or None,
        "part": "; ".join(parts) if parts else None,
        "state": state,
        "plant_to_fiber_ratio": ratio.group(1) if ratio else None,
        "is_mixed_dye": mixed,
    }


def parse_sci_name(raw: str) -> dict:
    """Split a scientific name into genus, binomial and rank, dropping authors."""
    s = raw.strip()
    tokens = s.split()
    genus = tokens[0]
    cultivar = re.search(r"'([^']+)'", s)
    if len(tokens) > 1 and tokens[1].lower() in {"spp.", "sp."}:
        return {"genus": genus, "binomial": None, "rank": "genus", "cultivar": None}
    epithet = tokens[1].lower() if len(tokens) > 1 else None
    binomial = f"{genus} {epithet}" if epithet else None
    rank = "species"
    if re.search(r"\b(ssp\.|subsp\.|var\.)", s):
        rank = "infraspecific"
    if cultivar:
        rank = "cultivar"
    return {"genus": genus, "binomial": binomial, "rank": rank,
            "cultivar": cultivar.group(1) if cultivar else None}


# ---------------------------------------------------------------- parsing

CARD = re.compile(
    r'<div class="col color-swatch"\s+'
    r'data-cielab-l="(?P<L>-?[\d.]+)"\s+'
    r'data-cielab-a="(?P<a>-?[\d.]+)"\s+'
    r'data-cielab-b="(?P<b>-?[\d.]+)">'
    r'.*?fill:rgb\((?P<r>\d+),(?P<g>\d+),(?P<bl>\d+)\)'
    r'.*?<h5 class="card-title">(?P<title>.*?)</h5>'
    r'.*?<p class="card-text">\n(?P<text>.*?)\n\s*</p>'
    r'.*?href="/library/color/detail/(?P<id>\d+)/"',
    re.S,
)


def parse_card_text(text: str) -> tuple:
    """The template prints species, then four slots: ',' fiber ',' mordant.
    Missing values leave blank lines, so position tells fiber from mordant."""
    lines = [ln.strip() for ln in text.split("\n")]
    species, slots = lines[0], (lines[1:] + [""] * 4)[:4]
    fiber = slots[1] if slots[0] == "," else ""
    mordant = slots[3] if slots[2] == "," else ""
    return html.unescape(species), html.unescape(fiber), html.unescape(mordant)


def main(src: Path, out_dir: Path) -> None:
    raw = src.read_text(encoding="utf-8")
    rows = []
    for m in CARD.finditer(raw):
        species, fiber, mordant = parse_card_text(m.group("text"))
        r, g, b = (int(m.group(k)) for k in ("r", "g", "bl"))
        row = {
            "sample_id": int(m.group("id")),
            "L": float(m.group("L")),
            "a": float(m.group("a")),
            "b": float(m.group("b")),
            "hex": f"#{r:02x}{g:02x}{b:02x}",
            "sci_name_raw": species,
            "fiber_raw": fiber or None,
            "mordant_raw": mordant or None,
            "mordant_class": mordant_class(mordant),
            "fiber_class": "wool" if "wool" in fiber.lower() else ("not recorded" if not fiber else "other"),
            "detail_url": f"https://library-biocolour.2.rahtiapp.fi/library/color/detail/{m.group('id')}/",
        }
        row.update(parse_label(html.unescape(m.group("title"))))
        rows.append(row)

    samples = pd.DataFrame(rows).sort_values("sample_id")

    # Polar form of a*b* (chroma, hue angle) and a coarse hue family for filters.
    samples["C"] = (samples["a"] ** 2 + samples["b"] ** 2) ** 0.5
    samples["h"] = (np.degrees(np.arctan2(samples["b"], samples["a"])) + 360) % 360
    samples["hue_family"] = [hue_family(c, h) for c, h in zip(samples["C"], samples["h"])]
    samples[["C", "h"]] = samples[["C", "h"]].round(1)
    # Near-white readings on dyed wool are suspect (likely an undyed or blank measurement).
    samples["qc_flag"] = np.where(samples["L"] > 90, "near-white: check detail page", None)

    taxa = (samples.groupby("sci_name_raw")
            .agg(common_names=("common_name", lambda s: "; ".join(sorted(set(s)))),
                 n_samples=("sample_id", "size"))
            .reset_index())
    name_parts = taxa["sci_name_raw"].apply(parse_sci_name).apply(pd.Series)
    taxa = pd.concat([taxa, name_parts], axis=1)
    taxa["group"] = taxa["genus"].map(GENUS_TO_GROUP).fillna("vascular plant")
    taxa.insert(0, "taxon_id", range(1, len(taxa) + 1))

    samples = samples.merge(taxa[["sci_name_raw", "taxon_id", "group"]], on="sci_name_raw")
    cols = ["sample_id", "taxon_id", "group", "sci_name_raw", "common_name", "label",
            "replicate", "prep_raw", "part", "state", "plant_to_fiber_ratio",
            "is_mixed_dye", "fiber_raw", "fiber_class", "mordant_raw",
            "mordant_class", "L", "a", "b", "C", "h", "hue_family", "hex",
            "qc_flag", "detail_url"]
    samples = samples[cols].sort_values("sample_id")

    out_dir.mkdir(parents=True, exist_ok=True)
    samples.to_csv(out_dir / "samples.csv", index=False)
    taxa.to_csv(out_dir / "taxa.csv", index=False)
    print(f"{len(samples)} samples, {len(taxa)} taxa -> {out_dir}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--html", type=Path, default=DATA / "pull.html")
    p.add_argument("--out", type=Path, default=DATA)
    a = p.parse_args()
    if not a.html.exists():
        p.error(f"input not found: {a.html}")
    main(a.html, a.out)
