# Validation

Tested on 14 November 2026: Ubuntu in WSL2, Linux x86_64, Docker 28.4.0, Java 17, Nextflow 24.10.4, Python 3.11.11.

## Run checks

- All four Nextflow processes completed with exit code 0. See `expected/trace.tsv`.
- A fresh folder with no Nextflow cache ran all four processes again.
- All 19 report files matched the reference exactly by SHA-256, including CSV, JSON, PNG, and HTML files.
- A resume test reused all four processes. All report files still matched.
- Three tests passed: known shifted spectra, rejection of a wrong raw-file hash, and pair/QC integrity.

The final counts are 65 raw records, 58 condition matches, 38 unique structures, and 703 pairs.
Seven records fail ion-mode or precursor-type checks, 14 have too few fragments, and six duplicate a retained structure.

## Resources

One task runs at a time. Docker enforces one CPU and 2 GB RAM per task.
The local trace reports about 4 seconds of total process time, excluding setup and orchestration.
The largest sampled RSS was about 316 MB. Very short processes may be under-sampled.
Docker, WSL, and the Nextflow Java process need additional host memory.
The course's exact numeric resource limits were not supplied; these are the measured settings, not a claim against unknown limits.

## Reproducibility

The base image digest, Python version, all Python package versions, and direct system package versions are pinned.
`expected/image-id.txt` records the tested local image. This is not a published registry image.
OS dependency packages still come from Debian repositories. Keep the built image for the closest repeat.
Fresh source builds were tested on the same machine and architecture, not on all computers.
The analysis containers have no network access. Initial image building and optional data restoration need internet.

Expected files are a regression reference, not independent proof of scientific truth.
Pair correlations are descriptive because pairs share molecules.
