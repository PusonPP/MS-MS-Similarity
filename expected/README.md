# Expected output

Reference results for the bundled 65 records and default `params.json`.

- `preprocess/`: cleaned spectra and QC.
- `report/`: 703 pair scores, summaries, figures, cases, and offline HTML report.
- `trace.tsv`: successful local run.
- `fresh-trace.tsv`: successful run in a new folder.
- `image-id.txt`: tested Docker image ID.

Check a new run with:

```bash
python3 scripts/check_results.py results/report expected/report
```

The trace changes between runs and is not compared byte for byte.
