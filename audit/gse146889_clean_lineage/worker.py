#!/usr/bin/env python3
"""Sharded Monte Carlo worker for the GSE146889 clean-lineage audit."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from ripser import ripser
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

N_PCS = 50
PCA_SEED = 42

TASK_BASE_SEEDS: dict[str, int] = {
    "mixed_gauss": 42,
    "mixed_perm": 42,
    "tumor_gauss": 14688901,
    "normal_gauss": 14688902,
    "tumor_perm": 14688911,
    "normal_perm": 14688912,
    "q1_gauss": 14689001,
    "q2_gauss": 14689002,
    "q3_gauss": 14689003,
    "q4_gauss": 14689004,
    "q1_perm": 14689101,
    "q2_perm": 14689102,
    "q3_perm": 14689103,
    "q4_perm": 14689104,
    "bootstrap_mixed": 14689200,
    "bootstrap_tumor": 14689201,
    "bootstrap_normal": 14689202,
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def max_h1(point_cloud: np.ndarray) -> float:
    h1 = np.asarray(ripser(point_cloud, maxdim=1)["dgms"][1], dtype=float)
    if h1.size == 0:
        return 0.0
    finite = h1[np.isfinite(h1[:, 1])]
    if finite.size == 0:
        return 0.0
    return float(np.max(finite[:, 1] - finite[:, 0]))


def pca_max_h1(matrix: np.ndarray, *, standardize: bool) -> float:
    X = np.asarray(matrix, dtype=np.float64)
    if standardize:
        X = StandardScaler().fit_transform(X)
    n_components = min(N_PCS, X.shape[0] - 1, X.shape[1])
    projected = PCA(n_components=n_components, random_state=PCA_SEED).fit_transform(X)
    return max_h1(projected)


def rng_for(task: str, draw_index: int) -> np.random.Generator:
    return np.random.default_rng(
        np.random.SeedSequence([TASK_BASE_SEEDS[task], int(draw_index)])
    )


def compute_draw(task: str, draw_index: int, prepared: Any) -> float:
    if task.endswith("_gauss"):
        subset = task.removesuffix("_gauss")
        source = np.asarray(prepared[f"{subset}_std"], dtype=np.float64)
        rng = rng_for(task, draw_index)
        gaussian = rng.normal(size=source.shape)
        return pca_max_h1(gaussian, standardize=True)

    if task.endswith("_perm"):
        subset = task.removesuffix("_perm")
        source = np.asarray(prepared[f"{subset}_std"], dtype=np.float64)
        rng = rng_for(task, draw_index)
        permuted = rng.permuted(source, axis=0)
        return pca_max_h1(permuted, standardize=False)

    if task.startswith("bootstrap_"):
        subset = task.removeprefix("bootstrap_")
        source = np.asarray(prepared[f"{subset}_pca"], dtype=np.float64)
        rng = rng_for(task, draw_index)
        indices = rng.integers(0, source.shape[0], size=source.shape[0])
        return max_h1(source[indices])

    raise ValueError(f"Unknown task: {task}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared", required=True)
    parser.add_argument("--task", required=True, choices=sorted(TASK_BASE_SEEDS))
    parser.add_argument("--start", required=True, type=int)
    parser.add_argument("--count", required=True, type=int)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    prepared_path = Path(args.prepared)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prepared = np.load(prepared_path, allow_pickle=False)

    indices = np.arange(args.start, args.start + args.count, dtype=np.int64)
    values = np.empty(args.count, dtype=np.float64)
    for local_index, draw_index in enumerate(indices):
        values[local_index] = compute_draw(args.task, int(draw_index), prepared)
        if (local_index + 1) % max(1, min(50, args.count // 10 or 1)) == 0:
            print(
                f"{args.task}: {local_index + 1}/{args.count} "
                f"(global draw {int(draw_index)})",
                flush=True,
            )

    output_path = output_dir / f"{args.task}_{args.start:05d}_{args.count:05d}.npz"
    np.savez_compressed(
        output_path,
        task=np.asarray(args.task),
        start=np.asarray(args.start, dtype=np.int64),
        count=np.asarray(args.count, dtype=np.int64),
        indices=indices,
        values=values,
        base_seed=np.asarray(TASK_BASE_SEEDS[args.task], dtype=np.int64),
        prepared_sha256=np.asarray(sha256_file(prepared_path)),
        worker_sha256=np.asarray(sha256_file(Path(__file__))),
    )
    summary = {
        "task": args.task,
        "start": args.start,
        "count": args.count,
        "mean": float(values.mean()),
        "sd_population": float(values.std(ddof=0)),
        "min": float(values.min()),
        "max": float(values.max()),
        "output": str(output_path),
    }
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
