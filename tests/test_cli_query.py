"""CLI tests for kira query."""

from __future__ import annotations

import json

from typer.testing import CliRunner

from kira.cli import app


runner = CliRunner()


def test_query_json_does_not_crash():
    result = runner.invoke(app, ["query", "--disease", "schistosomiasis", "--format", "json"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert isinstance(payload, list)
    assert len(payload) > 0


def test_query_known_target_uses_real_orthologue_keys():
    result = runner.invoke(
        app,
        [
            "query",
            "--disease",
            "schistosomiasis",
            "--target",
            "Thioredoxin glutathione reductase",
            "--format",
            "json",
        ],
    )
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert len(payload) == 1
    row = payload[0]
    assert row["target"] == "Thioredoxin glutathione reductase"
    assert row["essentiality"] == 1.0
    assert row["human_orthologue"] == "Thioredoxin reductase 1 (human)"
    assert row["human_id"] == "CHEMBL3952"


def test_query_target_without_orthologue_uses_safe_fallbacks():
    result = runner.invoke(
        app,
        [
            "query",
            "--disease",
            "schistosomiasis",
            "--target",
            "Venus kinase receptor 2",
            "--format",
            "json",
        ],
    )
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert len(payload) == 1
    row = payload[0]
    assert row["human_orthologue"] == "none"
    assert row["human_id"] == "N/A"
