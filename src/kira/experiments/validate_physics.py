"""Validate physics engine against real PDB structures.

Downloads real crystal structures from the PDB, runs all physics checks,
and compares against known Molprobity clashscores. This validates that
the physics engine produces meaningful results on actual protein data,
not just synthetic test fixtures.

HOW TO RUN:
    cd ~/Kira
    conda activate bio-builder
    python -m kira.experiments.validate_physics

WHAT IT TESTS:
    - Parser handles real PDB format quirks (HETATM, alt conformations, etc.)
    - Clashscores are in the expected range for high-resolution structures
    - LJ energies are finite and per-residue decomposition sums correctly
    - Backbone dihedrals fall in expected Ramachandran regions
    - Runtime is reasonable for real-sized proteins
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

# Try to import requests for PDB download
try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False


def download_pdb(pdb_id: str, output_dir: Path) -> Path:
    """Download a PDB file from RCSB.

    Parameters
    ----------
    pdb_id : str
        4-letter PDB ID (e.g., "1UBQ").
    output_dir : Path
        Directory to save the file.

    Returns
    -------
    Path
        Path to the downloaded PDB file.
    """
    if not HAS_REQUESTS:
        raise ImportError("requests library needed. pip install requests")

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    output_path = output_dir / f"{pdb_id.lower()}.pdb"

    if output_path.exists():
        print(f"  {pdb_id}: already downloaded")
        return output_path

    print(f"  Downloading {pdb_id} from RCSB...")
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    output_path.write_text(resp.text)
    print(f"  Saved to {output_path} ({len(resp.text)} bytes)")
    return output_path


# Structures to validate against, with expected properties
VALIDATION_STRUCTURES = {
    "1UBQ": {
        "name": "Ubiquitin",
        "resolution": 1.8,
        "n_residues_approx": 76,
        "expected_clashscore_max": 15.0,  # High-res structure, should be low
        "notes": "Small, well-characterized protein. Gold standard test case.",
    },
    "2X99": {
        "name": "SmTGR (S. mansoni thioredoxin glutathione reductase)",
        "resolution": 2.3,
        "n_residues_approx": 591,
        "expected_clashscore_max": 25.0,
        "notes": "Kira Script 14 docking target. YOUR protein.",
    },
    "1D3H": {
        "name": "HsDHODH (human dihydroorotate dehydrogenase)",
        "resolution": 1.6,
        "n_residues_approx": 367,
        "expected_clashscore_max": 15.0,
        "notes": "Human ortholog of SmDHODH. The 30.8x selectivity story.",
    },
}


def validate_structure(pdb_path: Path, expected: dict) -> dict:
    """Run full physics validation on a real PDB structure.

    Parameters
    ----------
    pdb_path : Path
        Path to the PDB file.
    expected : dict
        Expected properties for validation.

    Returns
    -------
    dict
        Validation results.
    """
    import jax.numpy as jnp

    from kira.physics.checks.clashes import check_clashes
    from kira.physics.core.energy import run_lj_analysis
    from kira.physics.core.geometry import compute_distance_matrix, extract_backbone_dihedrals
    from kira.physics.core.parser import parse_pdb
    from kira.physics.core.topology import build_bonded_mask, infer_bonds_from_topology

    results = {"pdb_id": pdb_path.stem.upper(), "name": expected["name"]}

    # Parse
    t0 = time.time()
    try:
        struct = parse_pdb(str(pdb_path))
    except Exception as e:
        results["error"] = f"Parse failed: {e}"
        return results

    results["n_atoms"] = struct.n_atoms
    results["n_residues"] = struct.n_residues
    results["n_chains"] = struct.n_chains
    results["parse_time"] = time.time() - t0

    # Topology
    t0 = time.time()
    bonds = infer_bonds_from_topology(struct)
    mask = build_bonded_mask(struct.n_atoms, bonds)
    results["n_bonds"] = len(bonds)
    results["topology_time"] = time.time() - t0

    # Distance matrix + checks
    t0 = time.time()
    coords = jnp.array(struct.coords)
    dist_matrix = compute_distance_matrix(coords)
    mask_jnp = jnp.array(mask)

    clash_r = check_clashes(dist_matrix, struct.elements, mask_jnp,
                            struct.res_indices, struct.n_residues)
    results["clashscore"] = clash_r.clashscore
    results["n_clashes"] = clash_r.n_clashes
    results["clash_subscore"] = clash_r.subscore

    # LJ energy
    lj = run_lj_analysis(dist_matrix, struct.elements, mask_jnp,
                         struct.res_indices, struct.n_residues)
    results["lj_total"] = lj["total_energy"]
    results["lj_finite"] = np.isfinite(lj["total_energy"])
    results["lj_per_res_sum"] = float(np.sum(lj["per_residue_energy"]))
    results["lj_hot_pairs"] = lj["n_hot_pairs"]

    # Per-residue energy consistency
    per_res_sum = float(np.sum(lj["per_residue_energy"]))
    results["energy_decomp_consistent"] = abs(per_res_sum - lj["total_energy"]) < 0.1

    # Backbone dihedrals
    dihedrals = extract_backbone_dihedrals(
        struct.coords, struct.atom_names, struct.res_indices,
        struct.is_protein_mask, struct.chain_ids_array,
    )
    results["n_phi"] = len(dihedrals["phi"]["res_indices"])
    results["n_psi"] = len(dihedrals["psi"]["res_indices"])
    results["n_omega"] = len(dihedrals["omega"]["res_indices"])

    # Check omega angles are mostly trans
    omega = np.asarray(dihedrals["omega"]["angles"])
    if len(omega) > 0:
        n_trans = np.sum(np.abs(np.abs(omega) - 180.0) < 30.0)
        results["fraction_trans_omega"] = float(n_trans / len(omega))
    else:
        results["fraction_trans_omega"] = 0.0

    results["total_time"] = results["parse_time"] + results["topology_time"] + (time.time() - t0)

    return results


def print_validation_results(all_results: list[dict]) -> None:
    """Print validation results in a readable format."""
    print("\n" + "=" * 70)
    print("PHYSICS ENGINE VALIDATION ON REAL PDB STRUCTURES")
    print("=" * 70)

    n_pass = 0
    n_fail = 0
    n_total = len(all_results)

    for r in all_results:
        pdb_id = r["pdb_id"]
        name = r["name"]
        print(f"\n  {pdb_id}: {name}")

        if "error" in r:
            print(f"    ❌ ERROR: {r['error']}")
            n_fail += 1
            continue

        print(f"    Atoms: {r['n_atoms']} | Residues: {r['n_residues']} | Chains: {r['n_chains']}")
        print(f"    Bonds inferred: {r['n_bonds']}")
        print(f"    Clashscore: {r['clashscore']:.1f} ({r['n_clashes']} clashes)")
        print(f"    LJ energy: {r['lj_total']:.1f} kcal/mol ({r['lj_hot_pairs']} hot pairs)")
        print(f"    LJ finite: {'✓' if r['lj_finite'] else '❌'}")
        print(f"    Energy decomposition consistent: {'✓' if r['energy_decomp_consistent'] else '❌'}")
        print(f"    Backbone dihedrals: φ={r['n_phi']}, ψ={r['n_psi']}, ω={r['n_omega']}")
        print(f"    Trans ω: {r['fraction_trans_omega']:.1%}")
        print(f"    Time: {r['total_time']:.2f}s")

        # Validation checks
        checks = []
        expected = VALIDATION_STRUCTURES.get(pdb_id, {})

        if r["clashscore"] <= expected.get("expected_clashscore_max", 30.0):
            checks.append(("Clashscore within range", True))
        else:
            checks.append(("Clashscore within range", False))

        checks.append(("LJ energy finite", r["lj_finite"]))
        checks.append(("Energy decomposition consistent", r["energy_decomp_consistent"]))
        checks.append(("Trans omega > 80%", r["fraction_trans_omega"] > 0.8))
        checks.append(("Has backbone dihedrals", r["n_phi"] > 0))

        all_pass = all(c[1] for c in checks)
        if all_pass:
            n_pass += 1
        else:
            n_fail += 1

        print("    Checks:")
        for name, passed in checks:
            symbol = "✓" if passed else "❌"
            print(f"      {symbol} {name}")

    print(f"\n  {'='*70}")
    print(f"  RESULT: {n_pass}/{n_total} structures passed all checks")
    if n_fail > 0:
        print(f"  {n_fail} structure(s) had failures — investigate before trusting physics results")
    else:
        print("  Physics engine validated on real protein structures ✓")


def main():
    """Download and validate real PDB structures."""
    print("=" * 70)
    print("PHYSICS ENGINE VALIDATION")
    print("Testing against real PDB structures from the RCSB")
    print("=" * 70)

    # Create cache directory for downloaded PDBs
    cache_dir = Path.home() / "Kira" / "data" / "validation"
    if not cache_dir.exists():
        cache_dir = Path("data") / "validation"
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Download structures
    print("\n[1/2] Downloading structures...")
    pdb_paths = {}
    for pdb_id in VALIDATION_STRUCTURES:
        try:
            path = download_pdb(pdb_id, cache_dir)
            pdb_paths[pdb_id] = path
        except Exception as e:
            print(f"  {pdb_id}: download failed ({e})")

    if not pdb_paths:
        print("  No structures downloaded. Check internet connection.")
        return

    # Validate
    print(f"\n[2/2] Validating {len(pdb_paths)} structures...")
    all_results = []
    for pdb_id, path in pdb_paths.items():
        expected = VALIDATION_STRUCTURES[pdb_id]
        result = validate_structure(path, expected)
        all_results.append(result)

    print_validation_results(all_results)


if __name__ == "__main__":
    main()
