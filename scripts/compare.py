import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from common import read_json, write_csv, write_json


def compare(spectra, scores, params, out):
    out.mkdir(parents=True, exist_ok=True)
    config = read_json(params)
    records = {r["id"]: r for r in read_json(spectra)}
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=config["morgan_radius"],
        fpSize=config["morgan_bits"], includeChirality=config["include_chirality"])
    fingerprints = {key: generator.GetFingerprint(Chem.MolFromSmiles(r["canonical_smiles"])) for key, r in records.items()}
    table = pd.read_csv(scores)
    table["name_a"] = table.id_a.map(lambda k: records[k]["name"])
    table["name_b"] = table.id_b.map(lambda k: records[k]["name"])
    table["tanimoto"] = [DataStructs.TanimotoSimilarity(fingerprints[a], fingerprints[b]) for a, b in zip(table.id_a, table.id_b)]
    table["delta"] = table.modified_cosine - table.cosine
    write_csv(out / "pairs.csv", table)
    summary = {"pairs": len(table), "scores": {}, "mean_delta": float(table.delta.mean()),
               "increased_pairs": int((table.delta > 1e-8).sum()), "decreased_pairs": int((table.delta < -1e-8).sum())}
    for column in ["cosine", "modified_cosine"]:
        summary["scores"][column] = {"mean": float(table[column].mean()), "median": float(table[column].median()),
            "pearson_r": float(table[column].corr(table.tanimoto)),
            "spearman_rho": float(table[column].corr(table.tanimoto, method="spearman")),
            "high_spectrum_low_structure": int(((table[column] >= 0.7) & (table.tanimoto <= 0.3)).sum()),
            "low_spectrum_high_structure": int(((table[column] <= 0.3) & (table.tanimoto >= 0.7)).sum())}
    write_json(out / "summary.json", summary)
    # Fixed, descriptive case rules.
    first = table.sort_values(["delta", "id_a", "id_b"], ascending=[False, True, True]).iloc[0]
    other = table.drop(first.name).copy()
    other["gap"] = other.modified_cosine - other.tanimoto
    second = other.sort_values(["gap", "id_a", "id_b"], ascending=[False, True, True]).iloc[0]
    cases = []
    for rule, row in [("largest_modified_gain", first), ("largest_spectral_structure_gap", second)]:
        cases.append({"rule": rule, **{key: row[key].item() if hasattr(row[key], "item") else row[key] for key in table.columns}})
    write_json(out / "cases.json", cases)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for name in ["spectra", "scores", "params"]:
        p.add_argument(name)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    compare(a.spectra, a.scores, a.params, a.out)
