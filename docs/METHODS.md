# Methods

## Data

Use the first 65 sorted filenames with exact collision energy `20 eV` in Athens_Univ, MassBank 2026.03.
This rule chooses a small pilot sample without looking at similarity scores.
It is a convenience sample, not a random sample of all compounds.
The exact commit, URLs, and SHA-256 hashes are in `data/manifest.json`.
Raw records are included. No download is needed for analysis.

Keep Bruker maXis Impact, LC-ESI-QTOF, MS2, positive mode, [M+H]+, and 20 eV.
Do not mix collision energies.
An initial filename-only pilot at 10 eV left just nine valid molecules after QC.
We changed to fixed 20 eV to retain more fragment-rich spectra. The final sample rule uses metadata only.
This was an exploratory design change, not a preregistered analysis.
The proposed 65 / 61 / 48 counts were planning estimates, not observed results.

## Processing

1. Verify every raw file against the manifest.
2. Check conditions, SMILES, and precursor m/z.
3. Remove non-finite values, non-positive intensities, m/z below 10, and peaks at or above precursor minus 1.5 Da.
4. Merge exact duplicate m/z values. Normalize to the largest remaining fragment.
5. Keep peaks at least 1% of this maximum. Require at least five fragments.
6. Keep the first valid record per canonical SMILES, in manifest order. Ignore chirality for this 2D comparison.

No salt stripping, charge neutralization, or tautomer merging is used.
Every record has one final QC reason. `qc.csv` also keeps the measured conditions and peak counts.

## Scores

Use matchms 0.28.2 `CosineGreedy` and `ModifiedCosine`.
Tolerance: 0.02 Da. Intensity power: 1. m/z power: 0.
Both use greedy one-to-one peak matching. Modified cosine also permits a precursor mass shift.
Greedy assignment is not guaranteed to find the global optimum.

Use RDKit Morgan bit fingerprints: radius 2, 2048 bits, no chirality.
Tanimoto is the number of shared set bits divided by the number of set bits in either fingerprint.
Score every unordered pair of different retained structures, including zero scores.

## Interpretation

Report Pearson and Spearman correlations, score means and medians, and modified minus ordinary cosine.
Use descriptive cutoffs: high >= 0.7; low <= 0.3. These are not validated thresholds.
Do not use p-values: pairs share molecules and are not independent.

Case 1 has the largest modified score gain. Case 2 has the largest modified score minus Tanimoto gap, excluding case 1.
Break ties by accession. Candidate peak tables show possible mass matches, not confirmed fragment identities or the full accepted greedy assignment.

This sample and one low collision energy cannot represent all chemistry or instruments.
Fingerprint size, peak filters, and mass tolerance can change the result.
The workflow tests an association, not a causal relation or an identification rule.

## Sources

- [MassBank 2026.03](https://github.com/MassBank/MassBank-data/releases/tag/2026.03)
- [Fixed Athens data](https://github.com/MassBank/MassBank-data/tree/705afb7bccc3b2c42410a744eef73674716a60ef/Athens_Univ)
- [matchms ModifiedCosine source](https://github.com/matchms/matchms/blob/0.28.2/matchms/similarity/ModifiedCosine.py)
- [RDKit fingerprints](https://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints)
- [Nextflow containers](https://www.nextflow.io/docs/stable/container.html)
