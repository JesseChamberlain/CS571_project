"""
Download Finland's 21 biogeographical provinces (eliömaakunnat) and save a clean copy.

Source: FinBIF's GeoServer layer LajiMapData:biogeographical_provinces (the same
layer laji.fi uses for its map), via the workspace WFS endpoint. The data is the
LUOMUS / Syke 1:1 000 000 province boundaries joined to Finland's borders,
including territorial waters (Syke metadata: https://ckan.ymparisto.fi/dataset/eliomaakunnat).
Syke does not offer it as an open download, and FinBIF's global WFS endpoint
is disabled, so the workspace URL below is the one that works.

Usage (from anywhere; paths default to the project's data/ folder):
    python scripts/get_provinces.py

Outputs (in data/boundaries/):
    finbif_biogeographical_provinces_raw.geojson   untouched WFS response (EPSG:3067)
    provinces.gpkg                                 cleaned, EPSG:3067 (ETRS-TM35FIN)
    provinces_source.json                          source URL and download time, for citation

Columns in provinces.gpkg:
    code       Finnish abbreviation (A, V, U, ... InL); used as region_code downstream
    code_lat   Latin abbreviation (Al, Ab, N, ... Li)
    name_fi    Finnish name (e.g. Etelä-Häme)
    name_lat   Latin name (e.g. Tavastia australis)
    area_km2   area as given in the source (includes territorial waters)
"""
import argparse
import json
import time
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"

WFS = "https://geoserver.laji.fi/geoserver/LajiMapData/wfs"
LAYER = "LajiMapData:biogeographical_provinces"
COLUMNS = {
    "maakuntalyh_fi": "code",
    "maakuntalyh_lat": "code_lat",
    "maakunta_fi": "name_fi",
    "maakunta_lat": "name_lat",
    "areabr": "area_km2",
}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", type=Path, default=DATA / "boundaries")
    args = p.parse_args()

    import geopandas as gpd
    import requests

    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    params = dict(service="WFS", version="2.0.0", request="GetFeature",
                  typeNames=LAYER, outputFormat="application/json")
    r = requests.get(WFS, params=params, timeout=300)
    r.raise_for_status()
    raw = out / "finbif_biogeographical_provinces_raw.geojson"
    raw.write_bytes(r.content)

    g = gpd.read_file(raw)
    if g.crs is None or g.crs.to_epsg() != 3067:
        # GeoJSON has no CRS member here in some GeoServer versions; the layer is ETRS-TM35FIN.
        g = g.set_crs(3067, allow_override=True)
    g = g[[*COLUMNS, "geometry"]].rename(columns=COLUMNS)
    g["geometry"] = g.geometry.make_valid()
    g = g.sort_values("code").reset_index(drop=True)
    if len(g) != 21 or g["code"].duplicated().any():
        raise SystemExit(f"Expected 21 unique provinces, got {len(g)}: {g['code'].tolist()}")
    g.to_file(out / "provinces.gpkg", layer="provinces", driver="GPKG")

    (out / "provinces_source.json").write_text(json.dumps({
        "url": r.url,
        "layer": LAYER,
        "downloaded": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "attribution": "Biogeographical provinces: LUOMUS (Finnish Museum of Natural "
                       "History) and Syke, via FinBIF (laji.fi)",
    }, indent=2), encoding="utf-8")

    print(g.drop(columns="geometry").to_string(index=False))
    print(f"\n{len(g)} provinces, {g['area_km2'].sum():,.0f} km2 -> {out / 'provinces.gpkg'}")


if __name__ == "__main__":
    main()
