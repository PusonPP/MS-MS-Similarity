import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
from matchms import Spectrum, set_matchms_logger_level
from matchms.similarity import CosineGreedy, ModifiedCosine

from common import read_json, write_csv


def make_spectrum(record):
    peaks = np.asarray(record["peaks"])
    return Spectrum(mz=peaks[:, 0], intensities=peaks[:, 1],
                    metadata={"precursor_mz": record["precursor_mz"]}, metadata_harmonization=False)


def score(spectra_path, params, out):
    set_matchms_logger_level("ERROR")
    config = read_json(params)
    records = read_json(spectra_path)
    spectra = [make_spectrum(r) for r in records]
    settings = {"tolerance": config["tolerance_da"], "mz_power": config["mz_power"],
                "intensity_power": config["intensity_power"]}
    cosine, modified = CosineGreedy(**settings), ModifiedCosine(**settings)
    rows = []
    for i, j in combinations(range(len(records)), 2):
        c, m = cosine.pair(spectra[i], spectra[j]), modified.pair(spectra[i], spectra[j])
        rows.append({"id_a": records[i]["id"], "id_b": records[j]["id"],
                     "cosine": float(c["score"]), "modified_cosine": float(m["score"]),
                     "cosine_matches": int(c["matches"]), "modified_matches": int(m["matches"])})
    write_csv(out, rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("spectra")
    p.add_argument("params")
    p.add_argument("out", type=Path)
    a = p.parse_args()
    score(a.spectra, a.params, a.out)
