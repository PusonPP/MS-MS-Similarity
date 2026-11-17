"""Download the fixed pilot sample. Python standard library only."""
import concurrent.futures
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

COMMIT = "705afb7bccc3b2c42410a744eef73674716a60ef"
ROOT = Path(__file__).resolve().parents[1]


def download(item, select=False):
    url = f"https://raw.githubusercontent.com/MassBank/MassBank-data/{COMMIT}/{item}"
    data = urlopen(url, timeout=60).read()
    if select and b"AC$MASS_SPECTROMETRY: COLLISION_ENERGY 20 eV\n" not in data.replace(b"\r\n", b"\n"):
        return None
    path = ROOT / "data/raw" / Path(item).name
    path.write_bytes(data)
    return {"file": path.name, "sha256": hashlib.sha256(data).hexdigest(), "url": url}


if __name__ == "__main__":
    manifest = ROOT / "data/manifest.json"
    if manifest.exists():
        records = json.loads(manifest.read_text())["records"]
        paths = ["Athens_Univ/" + r["file"] for r in records]
    else:
        base = "https://api.github.com/repos/MassBank/MassBank-data/git/trees/"
        root = json.load(urlopen(base + COMMIT, timeout=60))
        directory = next(x["sha"] for x in root["tree"] if x["path"] == "Athens_Univ")
        tree = json.load(urlopen(base + directory, timeout=120))
        if tree.get("truncated"):
            raise ValueError("Incomplete source tree")
        paths = sorted("Athens_Univ/" + x["path"] for x in tree["tree"] if x["path"].endswith(".txt"))
    (ROOT / "data/raw").mkdir(parents=True, exist_ok=True)
    if manifest.exists():
        with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
            records_new = list(pool.map(download, paths))
    else:
        records_new = []
        for item in paths:
            record = download(item, select=True)
            if record:
                records_new.append(record)
            if len(records_new) == 65:
                break
        if len(records_new) != 65:
            raise ValueError("Not enough 20 eV records")
    if manifest.exists():
        if records_new != records:
            raise ValueError("Source checksum mismatch")
    else:
        manifest.write_text(json.dumps({"release": "2026.03", "commit": COMMIT,
            "selection": "First 65 sorted Athens_Univ filenames with exact COLLISION_ENERGY 20 eV; no score-based selection.",
            "records": records_new}, indent=2) + "\n")
    print(f"Verified {len(records_new)} raw records")
