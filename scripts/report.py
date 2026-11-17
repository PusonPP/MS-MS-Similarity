import argparse
import base64
from html import escape
from importlib.metadata import version
from pathlib import Path
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

from common import read_json, write_csv, write_json

BLUE, ORANGE = "#226487", "#b45a32"
plt.rcParams.update({"font.size": 12, "axes.spines.top": False, "axes.spines.right": False})


def embed(path):
    return f'<img src="data:image/png;base64,{base64.b64encode(path.read_bytes()).decode()}" alt="{escape(path.stem)}">'


def candidate_peaks(a, b, tolerance):
    shift = a["precursor_mz"] - b["precursor_mz"]
    candidates = []
    for ma, ia in a["peaks"]:
        for mb, ib in b["peaks"]:
            direct, shifted = abs(ma - mb), abs(ma - mb - shift)
            if min(direct, shifted) <= tolerance:
                candidates.append({"mz_a": ma, "mz_b": mb,
                    "type": "direct" if direct <= tolerance else "shifted",
                    "intensity_product": ia * ib,
                    "loss_a": a["precursor_mz"] - ma, "loss_b": b["precursor_mz"] - mb})
    return sorted(candidates, key=lambda r: -r["intensity_product"])


def report(clean, comparison, out):
    out.mkdir(parents=True, exist_ok=True)
    records = {r["id"]: r for r in read_json(clean / "spectra.json")}
    qc = read_json(clean / "qc_summary.json")
    summary = read_json(comparison / "summary.json")
    cases = read_json(comparison / "cases.json")
    config = read_json(clean / "parameters.json")
    pairs = pd.read_csv(comparison / "pairs.csv")
    for file in ["pairs.csv", "summary.json", "cases.json"]:
        shutil.copyfile(comparison / file, out / file)
    for file in ["qc.csv", "qc_summary.json", "parameters.json", "molecules.csv"]:
        shutil.copyfile(clean / file, out / file)
    versions = {package: version(package) for package in ["matchms", "rdkit", "numpy", "pandas", "matplotlib", "scipy"]}
    write_json(out / "versions.json", versions)
    figures = []
    for key, title, color in [("cosine", "Cosine", BLUE), ("modified_cosine", "Modified cosine", ORANGE)]:
        fig, ax = plt.subplots(figsize=(6, 4.8), layout="constrained")
        ax.scatter(pairs.tanimoto, pairs[key], s=19, alpha=0.55, color=color, edgecolors="none")
        ax.set(xlabel="Structure similarity (Tanimoto)", ylabel=f"{title} similarity", xlim=(-0.02, 1.02), ylim=(-0.02, 1.02))
        stats = summary["scores"][key]
        ax.set_title(f"{title}: r = {stats['pearson_r']:.3f}, rho = {stats['spearman_rho']:.3f}")
        ax.grid(alpha=0.15)
        path = out / f"{key}_scatter.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        figures.append(path)
    blocks, descriptions = [], []
    for number, case in enumerate(cases, 1):
        a, b = records[case["id_a"]], records[case["id_b"]]
        structure_path = out / f"case_{number}_structures.png"
        Draw.MolsToGridImage([Chem.MolFromSmiles(a["canonical_smiles"]), Chem.MolFromSmiles(b["canonical_smiles"])],
            molsPerRow=2, subImgSize=(500, 260), legends=[a["name"], b["name"]]).save(structure_path)
        fig, ax = plt.subplots(figsize=(10, 3.5), layout="constrained")
        for r, sign, color in [(a, 1, BLUE), (b, -1, ORANGE)]:
            peaks = np.asarray(r["peaks"])
            ax.vlines(peaks[:, 0], 0, sign * peaks[:, 1], color=color, linewidth=1.4, label=r["name"])
            for index in np.argsort(peaks[:, 1])[-2:]:
                ax.text(peaks[index, 0], sign * (peaks[index, 1] + 0.06), f"{peaks[index, 0]:.2f}",
                        fontsize=9, ha="center", va="bottom" if sign == 1 else "top")
        ax.axhline(0, color="black", linewidth=0.6)
        ax.set(xlabel="Fragment m/z", ylabel="Relative intensity", ylim=(-1.35, 1.35))
        ax.legend(loc="upper left", fontsize=10, frameon=False)
        mirror_path = out / f"case_{number}_mirror.png"
        fig.savefig(mirror_path, dpi=160)
        plt.close(fig)
        candidates = candidate_peaks(a, b, config["tolerance_da"])
        write_csv(out / f"case_{number}_candidate_peaks.csv", candidates)
        shifted = [r for r in candidates if r["type"] == "shifted"]
        text = (f"{a['name']} / {b['name']}: Cosine {case['cosine']:.3f}, Modified cosine "
                f"{case['modified_cosine']:.3f}, Tanimoto {case['tanimoto']:.3f}. "
                f"Matched peaks: {case['cosine_matches']} direct; {case['modified_matches']} modified. ")
        if shifted:
            p = shifted[0]
            text += (f"A strong shifted candidate is m/z {p['mz_a']:.4f} versus {p['mz_b']:.4f}. "
                     f"Their precursor losses are {p['loss_a']:.4f} and {p['loss_b']:.4f} Da. "
                     "Similar losses can raise the modified score even when whole structures differ. ")
        text += "Peak matches describe mass patterns; they do not prove a shared fragment structure."
        descriptions.append({"case": number, "rule": case["rule"], "text": text})
        blocks.append(f"<h2>Case {number}: {escape(a['name'])} / {escape(b['name'])}</h2><p>{escape(text)}</p>"
                      + embed(structure_path) + embed(mirror_path))
    write_json(out / "case_explanations.json", descriptions)
    c, m = summary["scores"]["cosine"], summary["scores"]["modified_cosine"]
    findings = (f"Cosine: Pearson r = {c['pearson_r']:.3f}; Spearman rho = {c['spearman_rho']:.3f}. "
                f"Modified cosine: r = {m['pearson_r']:.3f}; rho = {m['spearman_rho']:.3f}. "
                f"Modified cosine raised {summary['increased_pairs']} of {summary['pairs']} scores "
                f"and lowered {summary['decreased_pairs']}. Mean change: {summary['mean_delta']:.3f}.")
    limitation = (f"This is a small convenience sample at {config['conditions']['collision_energy']}. Molecules occur in many pairs, so pair scores are not independent. "
                  "Correlations are descriptive; no p-values or causal claims are used. Morgan fingerprints describe 2D substructures, "
                  "not full chemical identity. Chirality is ignored. Results depend on peak filters, tolerance, collision energy, "
                  "fingerprint settings, and greedy matching. A higher spectral score alone cannot prove a similar structure.")
    html = f"""<!doctype html><html lang="en"><meta charset="utf-8"><title>MS/MS similarity report</title>
<style>body{{max-width:1000px;margin:40px auto;padding:0 24px;font:17px/1.6 Arial;color:#20313b;background:white}}img{{max-width:100%}}.plots{{display:flex}}.plots img{{width:50%;object-fit:contain}}table{{border-collapse:collapse}}td,th{{padding:6px 14px;border-bottom:1px solid #ddd}}h1,h2{{line-height:1.2}}</style>
<h1>MS/MS spectra and molecular structure</h1><p>Puzheng Wang | BIOF 501</p>
<p>Question: Do similar MS/MS spectra imply similar molecular structures?</p>
<h2>Data and quality control</h2><p>MassBank 2026.03, University of Athens. {qc['raw_records']} raw records; {qc['conditions_pass']} with shared conditions; {qc['molecules']} unique structures; {qc['pairs']} pairs.</p>
<p>Sample: first 65 sorted Athens_Univ filenames with exact collision energy 20 eV. This is a pilot sample, not a random sample of MassBank. An initial 10 eV pilot left only nine valid molecules; the design was changed to 20 eV before the final comparison. The proposed 61 / 48 counts were not assumed.</p>
<p>Shared conditions: {escape(', '.join(config['conditions'].values()))}. Retain m/z &gt;= {config['min_mz']} and &lt; precursor minus {config['precursor_exclusion_da']} Da. Normalize remaining peaks to their maximum; keep intensity &gt;= {100 * config['min_relative_intensity']:.1f}% and at least {config['min_peaks']} peaks.</p>
{pd.DataFrame(list(qc['reasons'].items()), columns=['QC result','Records']).to_html(index=False)}
<h2>Workflow</h2><p>PREPROCESS_SPECTRA → SCORE_SPECTRA → COMPARE_STRUCTURES → BUILD_REPORT</p>
<p>matchms: greedy one-to-one matching, tolerance {config['tolerance_da']} Da, intensity power {config['intensity_power']}, m/z power {config['mz_power']}. Modified cosine also allows the precursor mass shift. RDKit: Morgan radius {config['morgan_radius']}, {config['morgan_bits']} bits, chirality {config['include_chirality']}; Tanimoto.</p>
<h2>Overall results</h2><p>{findings}</p><div class="plots">{''.join(embed(p) for p in figures)}</div>
<p>High spectrum / low structure (score &gt;= 0.7, Tanimoto &lt;= 0.3): Cosine {c['high_spectrum_low_structure']}; Modified cosine {m['high_spectrum_low_structure']}. Low spectrum / high structure (score &lt;= 0.3, Tanimoto &gt;= 0.7): Cosine {c['low_spectrum_high_structure']}; Modified cosine {m['low_spectrum_high_structure']}.</p>
<p>These cutoffs describe examples; they are not validated identification rules. Cases below use the largest modified gain and then the largest modified-score minus Tanimoto gap, excluding the first pair.</p>
{''.join(blocks)}<h2>Limits and interpretation</h2><p>{limitation}</p>
<h2>Reproduce</h2><p>See README.md. All four steps run in Nextflow with a fixed Docker environment. See pairs.csv, qc.csv, summary.json, parameters.json, versions.json and the run trace.</p>
<p>Sources: <a href="https://github.com/MassBank/MassBank-data/releases/tag/2026.03">MassBank 2026.03</a>; <a href="https://github.com/matchms/matchms/tree/0.28.2">matchms 0.28.2</a>; <a href="https://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints">RDKit fingerprints</a>. Raw records retain their authors and CC BY license.</p></html>"""
    (out / "report.html").write_text(html, encoding="utf-8")
    (out / "findings.txt").write_text(findings + "\n\n" + limitation + "\n", encoding="utf-8")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for name in ["clean", "comparison", "out"]:
        p.add_argument(name, type=Path)
    a = p.parse_args()
    report(a.clean, a.comparison, a.out)
