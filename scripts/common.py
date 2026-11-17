import json
from pathlib import Path

import pandas as pd


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def write_csv(path, rows):
    pd.DataFrame(rows).to_csv(path, index=False, float_format="%.10g")
