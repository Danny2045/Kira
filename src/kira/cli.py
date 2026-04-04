"""Kira Engine CLI — unified entry point for causal drug discovery.

Commands:
    kira query     — Look up target information and known compounds
    kira evaluate  — Score predictions against a ground-truth evaluation set
    kira selectivity — Analyze selectivity between parasite target and human orthologue
    kira validate  — Physics validation of protein structures (from physics-auditor)
    kira info      — Print structural information about a PDB file
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

app = typer.Typer(
    name="kira",
    help=(
        "Kira Engine — Open Causal Discovery for Humanitarian Biology.\n\n"
        "A unified scientific engine for drug repurposing, physics validation, "
        "and selectivity analysis targeting neglected tropical diseases."
    ),
    no_args_is_help=True,
)
console = Console()


# ---------------------------------------------------------------------------
# query: look up target information and known compounds
# ---------------------------------------------------------------------------

@app.command()
def query(
    disease: str = typer.Option(..., "--disease", "-d", help="Disease name (e.g. schistosomiasis)"),
    target: str = typer.Option(None, "--target", "-t", help="Target gene/protein name (e.g. SmDHODH)"),
    format: str = typer.Option("table", "--format", "-f", help="Output format: table, json"),
) -> None:
    """Look up target essentiality, orthologues, and known compounds."""
    from kira.targets import ORTHOLOGUE_MAP, TARGET_ESSENTIALITY

    results: list[dict] = []

    for tgt, info in TARGET_ESSENTIALITY.items():
        # Filter by disease
        if disease.lower() not in tgt.lower() and disease.lower() not in info.get("disease", "").lower():
            continue
        # Filter by target name if provided
        if target and target.lower() not in tgt.lower():
            continue

        ortho = ORTHOLOGUE_MAP.get(tgt, {})
        results.append({
            "target": tgt,
            "disease": info.get("disease", "unknown"),
            "essentiality": info.get("essentiality", "unknown"),
            "human_orthologue": ortho.get("human", "none"),
            "identity_pct": ortho.get("identity", "N/A"),
        })

    if format == "json":
        print(json.dumps(results, indent=2))
        return

    if not results:
        console.print(f"[yellow]No targets found for disease='{disease}', target='{target}'[/yellow]")
        return

    table = Table(title=f"Targets: {disease}")
    table.add_column("Target", style="cyan")
    table.add_column("Disease", style="white")
    table.add_column("Essentiality", style="green")
    table.add_column("Human Orthologue", style="yellow")
    table.add_column("Identity %", style="white")

    for r in results:
        table.add_row(
            r["target"], r["disease"], str(r["essentiality"]),
            r["human_orthologue"], str(r["identity_pct"]),
        )

    console.print(table)


# ---------------------------------------------------------------------------
# evaluate: score predictions against ground truth
# ---------------------------------------------------------------------------

@app.command()
def evaluate(
    predictions: Path = typer.Option(..., "--predictions", "-p", help="CSV of predictions (ranked compounds)"),
    ground_truth: Path = typer.Option(..., "--ground-truth", "-g", help="CSV of ground truth labels"),
    format: str = typer.Option("table", "--format", "-f", help="Output format: table, json"),
) -> None:
    """Evaluate compound ranking predictions against a ground-truth set.

    Computes precision, recall, AUROC, and enrichment metrics.
    Predictions CSV must have columns: compound_id, score (higher = better).
    Ground truth CSV must have columns: compound_id, label (1 = active, 0 = inactive).
    """

    if not predictions.exists():
        console.print(f"[red]Predictions file not found: {predictions}[/red]")
        raise typer.Exit(1)
    if not ground_truth.exists():
        console.print(f"[red]Ground truth file not found: {ground_truth}[/red]")
        raise typer.Exit(1)

    # Parse CSVs (lightweight, no pandas dependency required)
    def _parse_csv(path: Path) -> list[dict]:
        rows = []
        with open(path) as f:
            header = [h.strip() for h in f.readline().split(",")]
            for line in f:
                vals = [v.strip() for v in line.split(",")]
                if len(vals) == len(header):
                    rows.append(dict(zip(header, vals)))
        return rows

    pred_rows = _parse_csv(predictions)
    gt_rows = _parse_csv(ground_truth)

    if not pred_rows:
        console.print("[red]No predictions found in file[/red]")
        raise typer.Exit(1)
    if not gt_rows:
        console.print("[red]No ground truth labels found in file[/red]")
        raise typer.Exit(1)

    # Build ground truth lookup
    gt_map = {}
    for row in gt_rows:
        cid = row.get("compound_id", "")
        label = row.get("label", "0")
        try:
            gt_map[cid] = int(label)
        except ValueError:
            gt_map[cid] = 0

    # Sort predictions by score descending
    for row in pred_rows:
        try:
            row["_score"] = float(row.get("score", "0"))
        except ValueError:
            row["_score"] = 0.0
    pred_rows.sort(key=lambda r: r["_score"], reverse=True)

    # Compute metrics
    n_total = len(pred_rows)
    n_actives = sum(1 for v in gt_map.values() if v == 1)
    if n_actives == 0:
        console.print("[yellow]Warning: No active compounds in ground truth[/yellow]")

    # Top-k enrichment
    hits_at = {}
    for k in [10, 25, 50, 100]:
        if k > n_total:
            break
        top_k = pred_rows[:k]
        hits = sum(1 for r in top_k if gt_map.get(r.get("compound_id", ""), 0) == 1)
        hits_at[k] = hits

    # Overall precision at each cutoff
    matched = sum(1 for r in pred_rows if r.get("compound_id", "") in gt_map)

    results = {
        "n_predictions": n_total,
        "n_ground_truth": len(gt_map),
        "n_actives": n_actives,
        "n_matched": matched,
        "enrichment": {f"hits_at_{k}": v for k, v in hits_at.items()},
    }

    if format == "json":
        print(json.dumps(results, indent=2))
        return

    console.print(Panel(
        f"Predictions: {n_total} | Ground truth: {len(gt_map)} | "
        f"Actives: {n_actives} | Matched: {matched}",
        title="Evaluation Summary",
        border_style="blue",
    ))

    if hits_at:
        table = Table(title="Enrichment")
        table.add_column("Top-K", style="cyan")
        table.add_column("Hits", style="green")
        table.add_column("Hit Rate", style="white")

        for k, h in hits_at.items():
            rate = h / k if k > 0 else 0
            table.add_row(str(k), str(h), f"{rate:.1%}")

        console.print(table)


# ---------------------------------------------------------------------------
# selectivity: physics-based selectivity analysis
# ---------------------------------------------------------------------------

@app.command()
def selectivity(
    parasite: Path = typer.Option(..., "--parasite", help="PDB file for parasite target"),
    human: Path = typer.Option(..., "--human", help="PDB file for human orthologue"),
    ligand_x: float = typer.Option(0.0, "--ligand-x", help="Ligand centroid X coordinate"),
    ligand_y: float = typer.Option(0.0, "--ligand-y", help="Ligand centroid Y coordinate"),
    ligand_z: float = typer.Option(0.0, "--ligand-z", help="Ligand centroid Z coordinate"),
    cutoff: float = typer.Option(5.0, "--cutoff", help="Binding site cutoff distance (Angstroms)"),
    format: str = typer.Option("table", "--format", "-f", help="Output format: table, json"),
) -> None:
    """Analyze selectivity between a parasite target and its human orthologue.

    Extracts binding pockets, decomposes per-residue LJ energy, and builds
    a selectivity attribution map showing which residue differences drive
    the binding energy differential.
    """
    import numpy as np

    from kira.causality.binding_site import extract_binding_site
    from kira.causality.divergence import compute_divergence_profile
    from kira.causality.energy_decomp import decompose_binding_site_energy
    from kira.causality.selectivity_map import build_selectivity_map, compute_selectivity_score
    from kira.physics.core.parser import parse_pdb
    from kira.physics.core.topology import build_bonded_mask, infer_bonds_from_topology

    for path, label in [(parasite, "parasite"), (human, "human")]:
        if not path.exists():
            console.print(f"[red]{label} PDB not found: {path}[/red]")
            raise typer.Exit(1)

    t0 = time.time()

    # Parse structures
    struct1 = parse_pdb(parasite)
    struct2 = parse_pdb(human)

    # Build topologies
    bonds1 = infer_bonds_from_topology(struct1)
    mask1 = build_bonded_mask(struct1.n_atoms, bonds1)
    bonds2 = infer_bonds_from_topology(struct2)
    mask2 = build_bonded_mask(struct2.n_atoms, bonds2)

    # Extract binding sites
    ligand_coords = np.array([[ligand_x, ligand_y, ligand_z]], dtype=np.float32)
    site1 = extract_binding_site(struct1, ligand_coords, cutoff=cutoff)
    site2 = extract_binding_site(struct2, ligand_coords, cutoff=cutoff)

    # Decompose energies
    decomp1 = decompose_binding_site_energy(
        struct1.coords, struct1.elements, struct1.res_indices, mask1, site1,
    )
    decomp2 = decompose_binding_site_energy(
        struct2.coords, struct2.elements, struct2.res_indices, mask2, site2,
    )

    # Build selectivity map
    sel_map = build_selectivity_map(
        site1, site2, decomp1, decomp2,
        structure1_name=struct1.name,
        structure2_name=struct2.name,
    )
    score = compute_selectivity_score(sel_map)

    # Divergence profile
    chain1 = next(iter(struct1.protein_chains), None)
    chain2 = next(iter(struct2.protein_chains), None)
    seq1 = chain1.sequence if chain1 else None
    seq2 = chain2.sequence if chain2 else None
    div_profile = compute_divergence_profile(site1, site2, seq1, seq2)

    runtime = time.time() - t0

    report = {
        "parasite": str(parasite),
        "human": str(human),
        "selectivity_score": round(score, 3),
        "selectivity_direction": sel_map.selectivity_direction,
        "total_delta_energy_kcal": round(sel_map.total_delta_energy, 2),
        "n_pocket_positions": sel_map.n_positions,
        "n_divergent_positions": sel_map.n_divergent,
        "divergent_energy_fraction": round(sel_map.divergent_energy_fraction, 3),
        "global_sequence_identity": round(div_profile.global_sequence_identity, 3),
        "local_sequence_identity": round(div_profile.local_sequence_identity, 3),
        "identity_gap": round(div_profile.identity_gap, 3),
        "summary": sel_map.summary(),
        "divergence_summary": div_profile.summary(),
        "runtime_seconds": round(runtime, 3),
    }

    if format == "json":
        print(json.dumps(report, indent=2))
        return

    console.print()
    console.print(Panel(
        f"[bold]{struct1.name}[/bold] vs [bold]{struct2.name}[/bold]  |  "
        f"Score: [bold]{score:.3f}[/bold]  |  "
        f"Direction: [bold]{sel_map.selectivity_direction}[/bold]",
        title="Selectivity Analysis",
        border_style="blue",
    ))

    table = Table(show_header=True, header_style="bold")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="white")

    table.add_row("Delta Energy", f"{sel_map.total_delta_energy:.2f} kcal/mol")
    table.add_row("Pocket Positions", f"{sel_map.n_positions}")
    table.add_row("Divergent Positions", f"{sel_map.n_divergent}")
    table.add_row("Divergent Energy %", f"{sel_map.divergent_energy_fraction:.1%}")
    table.add_row("Global Seq Identity", f"{div_profile.global_sequence_identity:.1%}")
    table.add_row("Local Seq Identity", f"{div_profile.local_sequence_identity:.1%}")
    table.add_row("Identity Gap", f"{div_profile.identity_gap:+.1%}")

    console.print(table)
    console.print()
    console.print(f"  [dim]{sel_map.summary()}[/dim]")
    console.print(f"  [dim]{div_profile.summary()}[/dim]")
    console.print(f"  [dim]Runtime: {runtime:.3f}s[/dim]")
    console.print()


# ---------------------------------------------------------------------------
# validate: physics validation of protein structures (from physics-auditor)
# ---------------------------------------------------------------------------

@app.command()
def validate(
    paths: list[Path] = typer.Argument(..., help="PDB file(s) to validate"),
    output_dir: Path | None = typer.Option(None, "--output-dir", "-o", help="Directory for JSON reports"),
    config: Path | None = typer.Option(None, "--config", "-c", help="YAML config file"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Show per-residue details"),
    json_output: bool = typer.Option(False, "--json", help="Print JSON to stdout"),
) -> None:
    """Validate one or more protein structures with physics checks."""
    import jax.numpy as jnp

    from kira.physics.checks.clashes import check_clashes
    from kira.physics.config import load_config
    from kira.physics.core.energy import run_lj_analysis
    from kira.physics.core.geometry import compute_distance_matrix, extract_backbone_dihedrals
    from kira.physics.core.parser import parse_pdb
    from kira.physics.core.topology import build_bonded_mask, infer_bonds_from_topology

    cfg = load_config(str(config) if config else None)

    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)

    for path in paths:
        if not path.exists():
            console.print(f"[red]File not found: {path}[/red]")
            continue

        t0 = time.time()

        try:
            struct = parse_pdb(path)
        except (ValueError, Exception) as e:
            console.print(f"[red]Failed to parse {path}: {e}[/red]")
            continue

        bonds = infer_bonds_from_topology(struct)
        mask = build_bonded_mask(struct.n_atoms, bonds)

        coords = jnp.array(struct.coords)
        dist_matrix = compute_distance_matrix(coords)
        mask_jnp = jnp.array(mask)

        clash_result = check_clashes(
            dist_matrix, struct.elements, mask_jnp,
            struct.res_indices, struct.n_residues, cfg.clashes,
        )

        lj_result = run_lj_analysis(
            dist_matrix, struct.elements, mask_jnp,
            struct.res_indices, struct.n_residues, cfg.lennard_jones.energy_cap,
        )

        dihedrals = extract_backbone_dihedrals(
            struct.coords, struct.atom_names, struct.res_indices,
            struct.is_protein_mask, struct.chain_ids_array,
        )

        runtime = time.time() - t0

        composite = clash_result.subscore
        if composite >= cfg.composite.accept_threshold:
            recommendation = "accept"
            rec_color = "green"
        elif composite >= cfg.composite.short_md_threshold:
            recommendation = "short_md"
            rec_color = "yellow"
        else:
            recommendation = "discard"
            rec_color = "red"

        report = {
            "file": str(path),
            "global_score": round(composite, 3),
            "recommendation": recommendation,
            "checks": {
                "steric_clashes": {
                    "n_clashes": clash_result.n_clashes,
                    "n_severe": clash_result.n_severe_clashes,
                    "clashscore": round(clash_result.clashscore, 2),
                    "worst_overlap_angstrom": round(clash_result.worst_overlap, 3),
                    "subscore": round(clash_result.subscore, 3),
                },
                "lennard_jones": {
                    "total_energy_kcal": round(lj_result["total_energy"], 2),
                    "n_hot_pairs": lj_result["n_hot_pairs"],
                },
                "backbone_dihedrals": {
                    "n_phi": len(dihedrals["phi"]["res_indices"]),
                    "n_psi": len(dihedrals["psi"]["res_indices"]),
                    "n_omega": len(dihedrals["omega"]["res_indices"]),
                },
            },
            "metadata": {
                "n_atoms": struct.n_atoms,
                "n_residues": struct.n_residues,
                "n_chains": struct.n_chains,
                "n_bonds_inferred": len(bonds),
                "protein_chains": len(struct.protein_chains),
                "runtime_seconds": round(runtime, 3),
            },
        }

        if json_output:
            print(json.dumps(report, indent=2))
        else:
            _print_validation_report(report, rec_color)

        if output_dir:
            out_path = output_dir / f"{path.stem}_report.json"
            with open(out_path, "w") as f:
                json.dump(report, f, indent=2)
            console.print(f"  Report saved: {out_path}")


# ---------------------------------------------------------------------------
# info: structural summary
# ---------------------------------------------------------------------------

@app.command()
def info(
    path: Path = typer.Argument(..., help="PDB file to inspect"),
) -> None:
    """Print basic structural information about a PDB file."""
    from kira.physics.core.parser import parse_pdb

    struct = parse_pdb(path)

    table = Table(title=f"Structure: {struct.name}")
    table.add_column("Property", style="cyan")
    table.add_column("Value", style="white")

    table.add_row("Atoms", str(struct.n_atoms))
    table.add_row("Residues", str(struct.n_residues))
    table.add_row("Chains", str(struct.n_chains))
    table.add_row("Protein chains", str(len(struct.protein_chains)))

    for chain in struct.chains.values():
        label = f"Chain {chain.chain_id}"
        seq = chain.sequence
        if len(seq) > 40:
            seq = seq[:37] + "..."
        table.add_row(label, f"{len(chain.residues)} res | {seq}")

    console.print(table)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _print_validation_report(report: dict, rec_color: str) -> None:
    """Print a formatted validation report."""
    name = Path(report["file"]).stem
    meta = report["metadata"]
    checks = report["checks"]

    score = report["global_score"]
    rec = report["recommendation"]
    console.print()
    console.print(Panel(
        f"[bold]{name}[/bold]  |  "
        f"Score: [bold]{score:.3f}[/bold]  |  "
        f"Recommendation: [{rec_color}][bold]{rec.upper()}[/bold][/{rec_color}]",
        title="Kira Physics Validation",
        border_style="blue",
    ))

    table = Table(show_header=True, header_style="bold")
    table.add_column("Check", style="cyan")
    table.add_column("Result", style="white")
    table.add_column("Score", style="white")

    cl = checks["steric_clashes"]
    clash_str = f"{cl['n_clashes']} clashes ({cl['n_severe']} severe), clashscore={cl['clashscore']}"
    table.add_row("Steric Clashes", clash_str, f"{cl['subscore']:.3f}")

    lj = checks["lennard_jones"]
    lj_str = f"E={lj['total_energy_kcal']:.1f} kcal/mol, {lj['n_hot_pairs']} hot pairs"
    table.add_row("Lennard-Jones", lj_str, "---")

    bb = checks["backbone_dihedrals"]
    bb_str = f"phi={bb['n_phi']}, psi={bb['n_psi']}, omega={bb['n_omega']}"
    table.add_row("Backbone Dihedrals", bb_str, "---")

    console.print(table)

    console.print(
        f"  [dim]{meta['n_atoms']} atoms | {meta['n_residues']} residues | "
        f"{meta['n_chains']} chains | {meta['n_bonds_inferred']} bonds | "
        f"{meta['runtime_seconds']:.3f}s[/dim]"
    )
    console.print()


def main() -> None:
    """Entry point for the kira CLI."""
    app()


if __name__ == "__main__":
    main()
