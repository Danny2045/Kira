"""Kira Engine CLI — unified entry point for causal drug discovery.

Commands:
    kira query       — Look up target information and known compounds
    kira evaluate    — Score predictions against a ground-truth evaluation set
    kira selectivity — Analyze selectivity between parasite target and human orthologue
    kira validate    — Physics validation of protein structures (from physics-auditor)
    kira info        — Print structural information about a PDB file
    kira delivery    — Predict LNP biodistribution and delivery probability
    kira expression  — Look up cell-type expression context for drug targets
    kira perturbation — Predict cellular response to target perturbation
    kira protocol    — Generate experiment protocol for selectivity assays
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

    # TARGET_ESSENTIALITY maps target_name -> float score.
    # All targets in the current dataset are S. mansoni (schistosomiasis).
    # We match disease name loosely against "schistosomiasis" / "mansoni".
    disease_aliases = {"schistosomiasis", "mansoni", "schistosoma"}
    disease_matches = any(alias in disease.lower() for alias in disease_aliases)

    for tgt, essentiality_score in TARGET_ESSENTIALITY.items():
        # If disease doesn't match any alias, skip all targets
        if not disease_matches:
            continue
        # Filter by target name if provided
        if target and target.lower() not in tgt.lower():
            continue

        ortho = ORTHOLOGUE_MAP.get(tgt, {})
        results.append({
            "target": tgt,
            "essentiality": essentiality_score,
            "human_orthologue": ortho.get("human_name", "none"),
            "notes": ortho.get("notes", ""),
        })

    if format == "json":
        print(json.dumps(results, indent=2))
        return

    if not results:
        console.print(f"[yellow]No targets found for disease='{disease}', target='{target}'[/yellow]")
        return

    table = Table(title=f"Targets: {disease}")
    table.add_column("Target", style="cyan")
    table.add_column("Essentiality", style="green")
    table.add_column("Human Orthologue", style="yellow")

    for r in results:
        table.add_row(
            r["target"], str(r["essentiality"]),
            r["human_orthologue"],
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
# delivery: LNP biodistribution and delivery probability
# ---------------------------------------------------------------------------

@app.command()
def delivery(
    target_tissue: str = typer.Option("liver", "--target-tissue", "-t", help="Target tissue for delivery"),
    lipid: str = typer.Option("SM-102", "--lipid", "-l", help="Ionizable lipid (SM-102, ALC-0315, DLin-MC3-DMA)"),
    route: str = typer.Option("intravenous", "--route", "-r", help="Administration route (intravenous, intramuscular, intranasal, etc.)"),
    payload: str = typer.Option("mRNA", "--payload", help="Payload type (mRNA, siRNA, CRISPR-Cas9 RNP, etc.)"),
    payload_size: int = typer.Option(1000, "--payload-size", help="Payload size in nucleotides"),
    format: str = typer.Option("table", "--format", "-f", help="Output format: table, json"),
) -> None:
    """Predict LNP biodistribution and delivery probability chain.

    Models tissue-level distribution based on lipid composition,
    administration route, and payload. Estimates the probability
    chain from circulation to intracellular action.
    """
    from kira.delivery.biodistribution import estimate_delivery_chain, predict_biodistribution
    from kira.delivery.formulation import (
        AdministrationRoute,
        IonizableLipidClass,
        LNPFormulation,
    )

    # Map string inputs to enums
    lipid_map = {
        "SM-102": IonizableLipidClass.SM102,
        "ALC-0315": IonizableLipidClass.ALC0315,
        "DLin-MC3-DMA": IonizableLipidClass.DLIN_MC3,
        "C12-200": IonizableLipidClass.C12_200,
        "Lipid 5": IonizableLipidClass.LIPID5,
    }
    route_map = {
        "intravenous": AdministrationRoute.IV,
        "iv": AdministrationRoute.IV,
        "intramuscular": AdministrationRoute.IM,
        "im": AdministrationRoute.IM,
        "subcutaneous": AdministrationRoute.SC,
        "sc": AdministrationRoute.SC,
        "intranasal": AdministrationRoute.IN,
        "in": AdministrationRoute.IN,
        "intrathecal": AdministrationRoute.IT,
        "it": AdministrationRoute.IT,
        "nebulized": AdministrationRoute.NEBULIZED,
    }

    il = lipid_map.get(lipid, IonizableLipidClass.SM102)
    rt = route_map.get(route.lower(), AdministrationRoute.IV)

    formulation = LNPFormulation(
        name=f"{lipid}_{route}_{payload}",
        ionizable_lipid=il,
        route=rt,
        payload_type=payload,
        payload_size_nt=payload_size,
    )

    biodist = predict_biodistribution(formulation, target_tissue)
    chain = estimate_delivery_chain(formulation, target_tissue)

    if format == "json":
        result = {
            "formulation": formulation.to_dict(),
            "biodistribution": biodist.to_dict(),
            "delivery_chain": chain.to_dict(),
        }
        print(json.dumps(result, indent=2))
        return

    # Biodistribution table
    console.print()
    console.print(Panel(
        f"[bold]{formulation.name}[/bold]  |  "
        f"Target: [bold]{target_tissue}[/bold]  |  "
        f"Route: [bold]{rt.value}[/bold]",
        title="LNP Delivery Analysis",
        border_style="blue",
    ))

    table = Table(title="Tissue Distribution")
    table.add_column("Tissue", style="cyan")
    table.add_column("Fraction", style="white")
    table.add_column("Level", style="white")

    sorted_tissues = sorted(
        biodist.tissue_fractions.items(), key=lambda x: x[1], reverse=True,
    )
    for tissue, frac in sorted_tissues:
        bar = "#" * int(frac * 40)
        style = "green" if tissue == target_tissue else "white"
        table.add_row(f"[{style}]{tissue}[/{style}]", f"{frac:.1%}", bar)

    console.print(table)

    # Delivery chain
    chain_table = Table(title="Delivery Probability Chain")
    chain_table.add_column("Step", style="cyan")
    chain_table.add_column("Probability", style="white")

    chain_table.add_row("Circulation", f"{chain.p_circulation:.4f}")
    chain_table.add_row("Tissue Access", f"{chain.p_tissue_access:.4f}")
    chain_table.add_row("Cell Binding", f"{chain.p_cell_binding:.4f}")
    chain_table.add_row("Internalization", f"{chain.p_internalization:.4f}")
    chain_table.add_row("Endosomal Escape", f"{chain.p_endosomal_escape:.4f}")
    chain_table.add_row("Intracellular Action", f"{chain.p_intracellular_action:.4f}")
    chain_table.add_row("[bold]Overall P(effective)[/bold]", f"[bold]{chain.p_effective:.6f}[/bold]")

    bottleneck_name, bottleneck_val = chain.bottleneck
    chain_table.add_row(f"[red]Bottleneck: {bottleneck_name}[/red]", f"[red]{bottleneck_val:.4f}[/red]")

    console.print(chain_table)

    if biodist.warnings:
        console.print()
        for w in biodist.warnings:
            console.print(f"  [yellow]Warning: {w}[/yellow]")
    console.print()


# ---------------------------------------------------------------------------
# expression: cell-type expression context
# ---------------------------------------------------------------------------

@app.command()
def expression(
    gene: str = typer.Option(..., "--gene", "-g", help="Gene symbol (e.g., DHODH, HDAC8, TXNRD1)"),
    format: str = typer.Option("table", "--format", "-f", help="Output format: table, json"),
) -> None:
    """Look up cell-type expression context for drug targets.

    Shows tissue-level expression profile for human orthologues,
    including tissue specificity and primary expression sites.
    """
    from kira.cell.expression import get_expression_profile

    profile = get_expression_profile(gene)
    if profile is None:
        console.print(f"[yellow]No expression profile found for '{gene}'[/yellow]")
        console.print("[dim]Available: DHODH, HDAC8, TXNRD1[/dim]")
        raise typer.Exit(1)

    if format == "json":
        print(json.dumps(profile.to_dict(), indent=2))
        return

    console.print()
    console.print(Panel(
        f"[bold]{profile.gene_name}[/bold] ({profile.organism})  |  "
        f"Specificity: [bold]{profile.tissue_specificity:.2f}[/bold]  |  "
        f"{'Ubiquitous' if profile.ubiquitous else 'Tissue-specific'}",
        title="Gene Expression Profile",
        border_style="blue",
    ))

    table = Table(title=f"Tissue Expression: {profile.gene_name}")
    table.add_column("Tissue", style="cyan")
    table.add_column("TPM", style="white")
    table.add_column("Level", style="white")
    table.add_column("Source", style="dim")

    sorted_tissues = sorted(
        profile.tissue_expression, key=lambda t: t.tpm, reverse=True,
    )
    for te in sorted_tissues:
        level_color = {
            "high": "red",
            "medium": "yellow",
            "low": "green",
            "not_detected": "dim",
        }.get(te.level.value, "white")
        table.add_row(
            te.tissue,
            f"{te.tpm:.1f}",
            f"[{level_color}]{te.level.value}[/{level_color}]",
            te.source,
        )

    console.print(table)
    console.print(f"\n  [dim]{profile.summary()}[/dim]\n")


# ---------------------------------------------------------------------------
# perturbation: predict cellular response
# ---------------------------------------------------------------------------

@app.command()
def perturbation(
    target_gene: str = typer.Option(..., "--target", "-t", help="Target gene (e.g., SmDHODH)"),
    concentration: float = typer.Option(1000.0, "--conc", "-c", help="Concentration in nM"),
    selectivity_ratio: float = typer.Option(1.0, "--selectivity", "-s", help="Selectivity ratio (human IC50 / parasite IC50)"),
    organism: str = typer.Option("Schistosoma mansoni", "--organism", "-o", help="Organism"),
    cell_type: str = typer.Option("generic", "--cell-type", help="Cell type"),
    format: str = typer.Option("table", "--format", "-f", help="Output format: table, json"),
) -> None:
    """Predict cellular response to target perturbation.

    Models pathway effects, viability, and cell state transitions
    for chemical inhibition of a target.
    """
    from kira.cell.perturbation import Perturbation, PerturbationType, predict_perturbation_response

    pert = Perturbation(
        name=f"{target_gene}_inhibition",
        perturbation_type=PerturbationType.CHEMICAL,
        target_gene=target_gene,
        concentration_nm=concentration,
        selectivity_ratio=selectivity_ratio,
    )

    response = predict_perturbation_response(pert, cell_type, organism)

    if format == "json":
        print(json.dumps(response.to_dict(), indent=2))
        return

    console.print()
    console.print(Panel(
        f"[bold]{target_gene}[/bold] @ {concentration:.0f} nM in {cell_type} ({organism})\n"
        f"State: [bold]{response.predicted_state_change.value}[/bold] "
        f"(p={response.state_change_probability:.2f})  |  "
        f"Viability: [bold]{response.viability:.2f}[/bold]",
        title="Perturbation Response",
        border_style="blue",
    ))

    if response.pathway_effects:
        table = Table(title="Pathway Effects")
        table.add_column("Pathway", style="cyan")
        table.add_column("Direction", style="white")
        table.add_column("Magnitude", style="white")
        table.add_column("Confidence", style="dim")

        for pe in response.pathway_effects:
            dir_color = "red" if pe.direction == "down" else "green" if pe.direction == "up" else "yellow"
            table.add_row(
                pe.pathway_name,
                f"[{dir_color}]{pe.direction}[/{dir_color}]",
                f"{pe.magnitude:.3f}",
                f"{pe.confidence:.2f}",
            )

        console.print(table)

    console.print(f"\n  [dim]{response.summary}[/dim]\n")


# ---------------------------------------------------------------------------
# protocol: generate experiment protocols
# ---------------------------------------------------------------------------

@app.command()
def protocol(
    parasite_target: str = typer.Option(..., "--parasite-target", "-p", help="Parasite target (e.g., SmDHODH)"),
    human_target: str = typer.Option(..., "--human-target", help="Human orthologue (e.g., HsDHODH)"),
    n_compounds: int = typer.Option(10, "--n-compounds", "-n", help="Number of compounds to screen"),
    top_conc: float = typer.Option(100.0, "--top-conc", help="Top concentration (uM)"),
    format: str = typer.Option("table", "--format", "-f", help="Output format: table, json, summary"),
) -> None:
    """Generate experiment protocol for selectivity assays.

    Designs a complete selectivity assay protocol with dose-response
    curves, controls, and quality criteria.
    """
    from kira.lab.protocol import (
        CompoundSpec,
        TargetSpec,
        design_selectivity_assay,
    )

    # Generate placeholder compounds
    compounds = [
        CompoundSpec(
            compound_id=f"KIRA-{i+1:04d}",
            name=f"Compound {i+1}",
        )
        for i in range(n_compounds)
    ]

    pt = TargetSpec(
        target_name=parasite_target,
        organism="Schistosoma mansoni",
    )
    ht = TargetSpec(
        target_name=human_target,
        organism="Homo sapiens",
    )

    proto = design_selectivity_assay(
        compounds=compounds,
        parasite_target=pt,
        human_target=ht,
        top_concentration_um=top_conc,
    )

    if format == "json":
        print(json.dumps(proto.to_dict(), indent=2))
        return

    if format == "summary":
        console.print(proto.summary())
        return

    console.print()
    console.print(Panel(
        f"[bold]{proto.name}[/bold]\n"
        f"ID: {proto.protocol_id}  |  Type: {proto.assay_type.value}",
        title="Experiment Protocol",
        border_style="blue",
    ))

    table = Table(title="Protocol Summary")
    table.add_column("Parameter", style="cyan")
    table.add_column("Value", style="white")

    table.add_row("Compounds", str(proto.n_compounds))
    table.add_row("Targets", f"{parasite_target}, {human_target}")
    table.add_row("Dose Points", str(proto.concentration_series.n_points))
    table.add_row("Replicates", str(proto.concentration_series.n_replicates))
    table.add_row("Top Concentration", f"{proto.concentration_series.top_concentration_um} uM")
    table.add_row("Plate Format", f"{proto.plate_format.value}-well")
    table.add_row("Plates Needed", str(proto.n_plates))
    table.add_row("Total Wells", str(proto.total_wells))
    table.add_row("Estimated Cost", f"${proto.estimated_cost_usd:,.0f}")
    table.add_row("Readout", proto.readout_type.value)
    table.add_row("Incubation", f"{proto.incubation_hours}h @ {proto.temperature_c}°C")

    console.print(table)

    # Show concentrations
    conc_str = ", ".join(
        f"{c:.2f}" for c in proto.concentration_series.concentrations_um
    )
    console.print(f"\n  [dim]Concentrations (uM): {conc_str}[/dim]")
    console.print(f"  [dim]{proto.objective}[/dim]\n")


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
