"""Guards that production modules do not perform network I/O at import time.

The `chembl_webresource_client` package fetches its API schema synchronously at
import time (see `new_client.py:65` in upstream). On 2026-05-20 an EBI ChEMBL
outage returned HTTP 500 for `/api/data/spore`, which crashed pytest
*collection* on `tests/test_selectivity_v4_data.py` because
`kira.experiments.selectivity_v4_data` imported the chembl client at module
top level. CI on `main` had been green only because EBI happened to be healthy
at every prior merge — a latent reproducibility hole.

The fix deferred the chembl import inside `_fetch_chembl_record`. This file
asserts the invariant so that any future regression moving the import back to
module level fails loudly here.
"""

from __future__ import annotations

import importlib
import sys


def test_selectivity_v4_data_imports_without_chembl_client_loaded() -> None:
    """Importing the module must not pull in chembl_webresource_client.

    If `chembl_webresource_client` appears in `sys.modules` immediately after
    importing `kira.experiments.selectivity_v4_data`, someone has moved the
    chembl import back to module top level, which makes import-time network
    success a CI prerequisite. Reject that.
    """
    for mod in [
        "kira.experiments.selectivity_v4_data",
        "chembl_webresource_client",
        "chembl_webresource_client.new_client",
    ]:
        sys.modules.pop(mod, None)

    importlib.import_module("kira.experiments.selectivity_v4_data")

    assert "chembl_webresource_client.new_client" not in sys.modules, (
        "kira.experiments.selectivity_v4_data must not import "
        "chembl_webresource_client at module level — that triggers a "
        "synchronous EBI schema fetch and makes pytest collection depend "
        "on an external service being up. Defer the import into "
        "_fetch_chembl_record (after the cache short-circuit, inside the "
        "existing try/except)."
    )
