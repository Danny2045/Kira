"""Run the archived Script 20 LODO model from the repo root.

Per the Disposition Plan, archived/scripts/20_selectivity_model_v2.py is
preserved unchanged as historical evidence. This wrapper executes the
archived script in-process with a virtual __file__ pointing to the repo
root, so that Script 20's ``os.path.dirname(__file__)``-rooted data
lookups resolve to the repo's actual ``data/`` directory.

The Script 20 output JSON (data/models/model_v2_results.json) is
regenerated against the corrected ESM-2 cache from PR #17.
"""
from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "archived" / "scripts" / "20_selectivity_model_v2.py"

if not SCRIPT.exists():
    raise FileNotFoundError(f"Script not found: {SCRIPT}")

os.chdir(REPO_ROOT)
virtual_file = str(REPO_ROOT / "20_selectivity_model_v2.py")
source = SCRIPT.read_text()
exec(
    compile(source, virtual_file, "exec"),
    {"__name__": "__main__", "__file__": virtual_file},
)
