import argparse
import hashlib
from pathlib import Path

import numpy as np
from rdkit import Chem

from common import read_json, write_csv, write_json


def parse_record(path):
    tags, peaks, in_peaks = {}, [], False
    for line in path.read_text().splitlines():
        if line == "//":
            break
        if line.startswith("PK$PEAK:"):
            in_peaks = True
        elif in_peaks:
            values = line.split()
            peaks.append([float(values[0]), float(values[1])])
        elif ": " in line:
            key, value = line.split(": ", 1)
            if key in {"AC$MASS_SPECTROMETRY", "MS$FOCUSED_ION"}:
                subkey, value = value.split(" ", 1)
                key += ":" + subkey
            tags.setdefault(key, value)
    mapping = {"id": "ACCESSION", "name": "CH$NAME", "smiles": "CH$SMILES",
        "instrument": "AC$INSTRUMENT", "instrument_type": "AC$INSTRUMENT_TYPE",
        "ms_type": "AC$MASS_SPECTROMETRY:MS_TYPE", "ion_mode": "AC$MASS_SPECTROMETRY:ION_MODE",
        "collision_energy": "AC$MASS_SPECTROMETRY:COLLISION_ENERGY",
        "precursor_type": "MS$FOCUSED_ION:PRECURSOR_TYPE"}
    record = {key: tags.get(tag, "") for key, tag in mapping.items()}
    record["precursor_mz"] = float(tags.get("MS$FOCUSED_ION:PRECURSOR_M/Z", "nan"))
    record["peaks"] = peaks
    return record


def preprocess(raw, manifest, config, out):
    out.mkdir(parents=True, exist_ok=True)
    kept, qc, seen = [], [], set()
    entries = read_json(manifest)["records"]
    if {p.name for p in raw.glob("*.txt")} != {e["file"] for e in entries}:
        raise ValueError("Raw files do not match manifest")
    for entry in entries:
        path = raw / entry["file"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            raise ValueError(f"Checksum mismatch: {path.name}")
        r = parse_record(path)
        row = {k: v for k, v in r.items() if k != "peaks"}
        row.update(raw_peaks=len(r["peaks"]), clean_peaks=0, conditions_pass=False, reason="kept")
        mismatches = [k for k, v in config["conditions"].items() if r[k] != v]
        if mismatches:
            row["reason"] = "condition:" + ",".join(mismatches)
        else:
            row["conditions_pass"] = True
            mol = Chem.MolFromSmiles(r["smiles"])
            if mol is None or mol.GetNumAtoms() == 0:
                row["reason"] = "invalid_structure"
            elif not np.isfinite(r["precursor_mz"]) or r["precursor_mz"] <= 0:
                row["reason"] = "invalid_precursor"
            else:
                r["canonical_smiles"] = Chem.MolToSmiles(mol, isomericSmiles=False)
                peaks = np.asarray(r["peaks"], dtype=float).reshape(-1, 2)
                valid = np.isfinite(peaks).all(axis=1) & (peaks[:, 1] > 0)
                valid &= (peaks[:, 0] >= config["min_mz"])
                valid &= (peaks[:, 0] < r["precursor_mz"] - config["precursor_exclusion_da"])
                peaks = peaks[valid]
                if len(peaks):
                    # Merge exact duplicate masses.
                    masses, inverse = np.unique(peaks[:, 0], return_inverse=True)
                    intensity = np.bincount(inverse, weights=peaks[:, 1])
                    intensity /= intensity.max()
                    mask = intensity >= config["min_relative_intensity"]
                    peaks = np.column_stack([masses[mask], intensity[mask]])
                row["clean_peaks"] = len(peaks)
                if len(peaks) < config["min_peaks"]:
                    row["reason"] = "too_few_fragments"
                elif r["canonical_smiles"] in seen:
                    row["reason"] = "duplicate_structure"
                else:
                    seen.add(r["canonical_smiles"])
                    r["peaks"] = peaks.tolist()
                    kept.append(r)
        qc.append(row)
    if len(kept) < 2:
        raise ValueError("Need at least two valid molecules")
    write_csv(out / "qc.csv", qc)
    write_json(out / "spectra.json", kept)
    write_csv(out / "molecules.csv", [{k: v for k, v in r.items() if k != "peaks"} for r in kept])
    counts = {"raw_records": len(qc), "conditions_pass": sum(r["conditions_pass"] for r in qc),
              "molecules": len(kept), "pairs": len(kept) * (len(kept) - 1) // 2}
    counts["reasons"] = {reason: sum(r["reason"] == reason for r in qc) for reason in sorted({r["reason"] for r in qc})}
    write_json(out / "qc_summary.json", counts)
    write_json(out / "parameters.json", config)
    print(counts)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("raw", type=Path)
    p.add_argument("manifest")
    p.add_argument("params")
    p.add_argument("out", type=Path)
    a = p.parse_args()
    preprocess(a.raw, a.manifest, read_json(a.params), a.out)
