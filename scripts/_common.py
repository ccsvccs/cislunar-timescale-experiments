"""Shared CLI mechanics only; core physics lives in src/tt_tdb."""
from pathlib import Path
import sys
# Let a clean checkout run the scripts without an editable install.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import json
import os
import platform
import importlib.metadata
import numpy as np

# A read-only home directory should not cause configuration writes or warnings.
os.environ.setdefault('MPLCONFIGDIR', '/tmp/tt-tdb-matplotlib')


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n')


def versions():
    packages = ('numpy','scipy','jplephem','spiceypy','pyerfa','mpmath','matplotlib','pytest')
    return {'python': platform.python_version(), 'platform': platform.platform(),
            'packages': {name:importlib.metadata.version(name) for name in packages}}


def make_plots_directory():
    import matplotlib
    matplotlib.use('Agg')
    p = Path('results/figures')
    p.mkdir(parents=True, exist_ok=True)
    return p


def csv_table(path, headers, values):
    import csv
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(np.asarray(values).tolist())
