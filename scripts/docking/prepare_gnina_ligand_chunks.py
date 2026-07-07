#!/usr/bin/env python3
"""Prepare 3D ligand SDF chunks for GNINA docking."""

from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
import os
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--chunk-size", default=500, type=int)
    parser.add_argument("--max-iters", default=120, type=int)
    parser.add_argument("--workers", default=max(1, min(8, os.cpu_count() or 1)), type=int)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def prepare_chunk(payload: tuple[int, int, list[dict[str, str]], Path, int]) -> dict[str, object]:
    chunk_idx, start_idx, rows, output_dir, max_iters = payload
    chunk_sdf = output_dir / f"fxia_filtered_all_chunk_{chunk_idx:04d}.sdf"
    writer = Chem.SDWriter(str(chunk_sdf))
    manifest_rows: list[dict[str, str]] = []
    failed: list[dict[str, str]] = []

    for offset, row in enumerate(rows, start=1):
        idx = start_idx + offset
        cid = row.get("candidate_id") or f"FXIA-GEN-{idx:06d}"
        smiles = row.get("parent_smiles") or row.get("canonical_smiles") or row.get("generated_smiles")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            failed.append({"candidate_id": cid, "reason": "parse_failed", "smiles": smiles or ""})
            continue

        mol = Chem.AddHs(mol)
        seed = (0xF11A + idx) % 2_147_483_647
        status = AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=True, maxAttempts=200)
        if status != 0:
            status = AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=True, maxAttempts=1000)
        if status != 0:
            failed.append({"candidate_id": cid, "reason": "embed_failed", "smiles": smiles or ""})
            continue

        try:
            if AllChem.MMFFHasAllMoleculeParams(mol):
                AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iters)
            else:
                AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
        except Exception:
            pass

        mol.SetProp("_Name", cid)
        for key in (
            "candidate_id",
            "source_models",
            "parent_inchi_key",
            "mw",
            "tpsa",
            "clogp",
            "hbd",
            "hba",
            "rotatable_bonds",
            "max_clinical_similarity",
            "max_patent_proxy_similarity",
        ):
            mol.SetProp(key, row.get(key, ""))

        writer.write(mol)
        manifest_rows.append(
            {
                "candidate_id": cid,
                "chunk_id": f"{chunk_idx:04d}",
                "chunk_sdf": chunk_sdf.name,
                "source_models": row.get("source_models", ""),
                "parent_smiles": row.get("parent_smiles", ""),
                "parent_inchi_key": row.get("parent_inchi_key", ""),
                "mw": row.get("mw", ""),
                "tpsa": row.get("tpsa", ""),
                "clogp": row.get("clogp", ""),
                "max_clinical_similarity": row.get("max_clinical_similarity", ""),
                "max_patent_proxy_similarity": row.get("max_patent_proxy_similarity", ""),
            }
        )

    writer.close()
    return {
        "chunk_id": f"{chunk_idx:04d}",
        "chunk_sdf": chunk_sdf.name,
        "ligands_written": len(manifest_rows),
        "ligands_failed": len(failed),
        "manifest_rows": manifest_rows,
        "failed": failed,
    }


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_csv = args.output_dir / "all_filtered_ligand_manifest.csv"
    summary_json = args.output_dir / "all_filtered_ligand_prep_summary.json"

    with args.input_csv.open() as handle:
        rows = list(csv.DictReader(handle))

    tasks = [
        (idx + 1, start, rows[start : start + args.chunk_size], args.output_dir, args.max_iters)
        for idx, start in enumerate(range(0, len(rows), args.chunk_size))
    ]

    if args.overwrite:
        for path in args.output_dir.glob("fxia_filtered_all_chunk_*.sdf"):
            path.unlink()

    results: list[dict[str, object]] = []
    if args.workers == 1:
        for task in tasks:
            results.append(prepare_chunk(task))
    else:
        with mp.Pool(processes=args.workers) as pool:
            for result in pool.imap_unordered(prepare_chunk, tasks):
                results.append(result)

    results.sort(key=lambda item: item["chunk_id"])
    manifest_rows: list[dict[str, str]] = []
    failed: list[dict[str, str]] = []
    for result in results:
        manifest_rows.extend(result["manifest_rows"])
        failed.extend(result["failed"])

    with manifest_csv.open("w", newline="") as handle:
        fieldnames = list(manifest_rows[0].keys()) if manifest_rows else []
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest_rows)

    summary = {
        "input_csv": str(args.input_csv),
        "input_rows": len(rows),
        "ligands_written": len(manifest_rows),
        "ligands_failed": len(failed),
        "chunk_size": args.chunk_size,
        "chunk_count": len(results),
        "workers": args.workers,
        "output_dir": str(args.output_dir),
        "manifest_csv": str(manifest_csv),
        "failed_examples": failed[:20],
    }
    summary_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
