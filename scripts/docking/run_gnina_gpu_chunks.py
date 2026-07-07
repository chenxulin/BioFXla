#!/usr/bin/env python3
"""Run GNINA docking chunks across one or more GPU devices."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunks-dir", required=True, type=Path)
    parser.add_argument("--receptor", required=True, type=Path)
    parser.add_argument("--autobox-ligand", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--devices", default="0,1,2,3,4,5,6,7")
    parser.add_argument("--gnina-bin", default="gnina")
    parser.add_argument("--exhaustiveness", default=2, type=int)
    parser.add_argument("--num-modes", default=1, type=int)
    parser.add_argument("--cpu", default=4, type=int)
    parser.add_argument("--cnn", default="fast")
    parser.add_argument("--autobox-add", default=6.0, type=float)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def run_chunk(task: tuple[Path, int, argparse.Namespace]) -> dict[str, object]:
    chunk_sdf, device, args = task
    args.output_dir.mkdir(parents=True, exist_ok=True)
    poses_dir = args.output_dir / "poses"
    logs_dir = args.output_dir / "logs"
    poses_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    stem = chunk_sdf.stem
    out_sdf = poses_dir / f"{stem}_gnina_poses.sdf"
    log_path = logs_dir / f"{stem}_gnina.log"
    meta_path = logs_dir / f"{stem}_gnina.meta.json"

    if meta_path.exists() and not args.overwrite:
        try:
            meta = json.loads(meta_path.read_text())
        except json.JSONDecodeError:
            meta = {}
        completed = (
            meta.get("returncode") == 0
            and out_sdf.exists()
            and out_sdf.stat().st_size > 0
        )
    else:
        completed = False

    if completed:
        return {
            "chunk": chunk_sdf.name,
            "device": device,
            "status": "skipped_existing",
            "output_sdf": str(out_sdf),
            "log": str(log_path),
            "returncode": 0,
            "elapsed_sec": 0.0,
        }

    command = [
        args.gnina_bin,
        "-r",
        str(args.receptor),
        "-l",
        str(chunk_sdf),
        "--autobox_ligand",
        str(args.autobox_ligand),
        "--autobox_add",
        str(args.autobox_add),
        "--exhaustiveness",
        str(args.exhaustiveness),
        "--num_modes",
        str(args.num_modes),
        "--cpu",
        str(args.cpu),
        "--device",
        str(device),
        "--cnn",
        args.cnn,
        "-o",
        str(out_sdf),
        "--log",
        str(log_path),
    ]

    started = time.time()
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(device)
    # GNINA sees the selected GPU as device 0 after CUDA_VISIBLE_DEVICES.
    command[command.index("--device") + 1] = "0"
    proc = subprocess.run(command, cwd=Path.cwd(), env=env, text=True)
    elapsed = time.time() - started
    result = {
        "chunk": chunk_sdf.name,
        "device": device,
        "status": "ok" if proc.returncode == 0 else "failed",
        "output_sdf": str(out_sdf),
        "log": str(log_path),
        "returncode": proc.returncode,
        "elapsed_sec": round(elapsed, 3),
        "command": command,
    }
    meta_path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    return result


def main() -> int:
    args = parse_args()
    args.chunks_dir = args.chunks_dir.resolve()
    args.receptor = args.receptor.resolve()
    args.autobox_ligand = args.autobox_ligand.resolve()
    args.output_dir = args.output_dir.resolve()
    devices = [int(item) for item in args.devices.split(",") if item.strip()]
    chunks = sorted(args.chunks_dir.glob("fxia_filtered_all_chunk_*.sdf"))
    if not chunks:
        raise SystemExit(f"No chunks found in {args.chunks_dir}")

    tasks = [(chunk, devices[idx % len(devices)], args) for idx, chunk in enumerate(chunks)]
    run_manifest = args.output_dir / "gnina_chunk_run_manifest.csv"
    results: list[dict[str, object]] = []

    with ThreadPoolExecutor(max_workers=len(devices)) as executor:
        futures = [executor.submit(run_chunk, task) for task in tasks]
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(json.dumps({k: result[k] for k in ("chunk", "device", "status", "returncode", "elapsed_sec")}, ensure_ascii=False), flush=True)

    results.sort(key=lambda item: item["chunk"])
    with run_manifest.open("w", newline="") as handle:
        fieldnames = ["chunk", "device", "status", "returncode", "elapsed_sec", "output_sdf", "log"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: result.get(key, "") for key in fieldnames} for result in results)

    failed = [result for result in results if result["returncode"] != 0]
    summary = {
        "chunk_count": len(chunks),
        "devices": devices,
        "failed_count": len(failed),
        "manifest": str(run_manifest),
    }
    (args.output_dir / "gnina_chunk_run_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
