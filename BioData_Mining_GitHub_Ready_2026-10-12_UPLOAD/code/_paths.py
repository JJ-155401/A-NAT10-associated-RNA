"""Shared repository paths for the reproducibility package.

Large public input files are intentionally kept outside version control.
Set AD_NAT10_DATA_DIR to the directory containing the downloaded inputs, or
use the default ``data/external`` directory inside this repository.
"""
from pathlib import Path
import os


REPO_ROOT = Path(__file__).resolve().parents[1]
EXTERNAL_DATA = Path(
    os.environ.get("AD_NAT10_DATA_DIR", str(REPO_ROOT / "data" / "external"))
).expanduser().resolve()
CURATED_TARGETS = REPO_ROOT / "data" / "human_ac4c_targets_curated.csv"
RESULTS_DIR = REPO_ROOT / "results"
FIGURES_DIR = REPO_ROOT / "figures"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
