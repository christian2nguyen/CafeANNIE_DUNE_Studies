#!/usr/bin/env python3
"""
plot_matrices.py
----------------
Reads all TMatrixT<double> objects from a ROOT file and produces:
  - A heatmap of the raw covariance matrix
  - A heatmap of the derived correlation matrix

Usage:
    python plot_matrices.py <input.root> [--outdir <output_dir>]

Requirements:
    pip install uproot numpy matplotlib
"""

from __future__ import annotations

import argparse
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import TwoSlopeNorm
from typing import Dict, Optional
import uproot


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_matrices(root_path: str) -> Dict[str, np.ndarray]:
    """Return all TMatrixT<double> objects found anywhere in the ROOT file."""
    matrices = {}

    def _visit(name, obj):
        classname = getattr(obj, "classname", None) or type(obj).__name__
        # uproot represents TMatrixT<double> as TMatrixD or similar
        if "TMatrixT" in classname or "TMatrixD" in classname:
            try:
                arr = obj.to_numpy()          # shape (nrows, ncols)
                matrices[name] = arr
            except Exception as exc:
                print(f"[WARN] Could not convert '{name}': {exc}")

    with uproot.open(root_path) as f:
        # iterate over all keys, including those inside directories/cycles
        for key in f.keys(recursive=True, cycle=False):
            obj = f[key]
            _visit(key, obj)

    if not matrices:
        print("[ERROR] No TMatrixT<double> objects found in the file.")
        sys.exit(1)

    print(f"Found {len(matrices)} matrix(es): {list(matrices.keys())}")
    return matrices


def cov_to_corr(cov: np.ndarray) -> np.ndarray:
    """Convert a covariance matrix to a correlation matrix."""
    std = np.sqrt(np.diag(cov))
    # Guard against zero-variance rows/columns
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = cov / np.outer(std, std)
    corr = np.where(np.isfinite(corr), corr, 0.0)
    np.fill_diagonal(corr, 1.0)
    return corr


def make_label(name: str) -> str:
    """Strip ROOT path separators and cycle numbers for display."""
    label = name.split(";")[0].replace("/", " / ")
    return label


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_matrix(
    ax: plt.Axes,
    mat: np.ndarray,
    title: str,
    cmap: str,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    center: Optional[float] = None,
    fmt: str = ".2g",
    annot_threshold: int = 20,
):
    """Draw a single heatmap on *ax*."""
    n = mat.shape[0]

    # Normalisation
    if center is not None:
        _vmin = vmin if vmin is not None else mat.min()
        _vmax = vmax if vmax is not None else mat.max()
        norm = TwoSlopeNorm(vmin=_vmin, vcenter=center, vmax=_vmax)
    else:
        norm = None

    im = ax.imshow(
        mat,
        cmap=cmap,
        aspect="auto",
        norm=norm,
        vmin=vmin if center is None else None,
        vmax=vmax if center is None else None,
        interpolation="nearest",
    )
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    # Annotate cells only when matrix is small enough to be readable
    if n <= annot_threshold:
        for i in range(n):
            for j in range(n):
                val = mat[i, j]
                text = format(val, fmt)
                # Pick contrasting text colour
                bg = im.norm(val)
                lum = plt.get_cmap(cmap)(bg)[0] * 0.299 + plt.get_cmap(cmap)(bg)[1] * 0.587 + plt.get_cmap(cmap)(bg)[2] * 0.114
                color = "white" if lum < 0.5 else "black"
                ax.text(j, i, text, ha="center", va="center", fontsize=7, color=color)

    ax.set_title(title, fontsize=11, fontweight="bold", pad=8)
    ax.set_xlabel("Column index")
    ax.set_ylabel("Row index")

    ticks = np.arange(n)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    if n <= 20:
        ax.set_xticklabels(ticks, fontsize=7)
        ax.set_yticklabels(ticks, fontsize=7)
    else:
        # Only show every Nth tick to avoid clutter
        step = max(1, n // 10)
        ax.xaxis.set_major_locator(mticker.MultipleLocator(step))
        ax.yaxis.set_major_locator(mticker.MultipleLocator(step))


def plot_pair(name: str, cov: np.ndarray, outdir: str):
    """Save a figure with the covariance and correlation matrices side by side."""
    corr = cov_to_corr(cov)
    label = make_label(name)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(label, fontsize=13, y=1.01)

    # --- Covariance (symmetric diverging palette around 0) ---
    abs_max = np.max(np.abs(cov))
    plot_matrix(
        axes[0],
        cov,
        title="Covariance matrix",
        cmap="RdBu_r",
        vmin=-abs_max,
        vmax=abs_max,
        center=0.0,
    )

    # --- Correlation (fixed [-1, 1] range) ---
    plot_matrix(
        axes[1],
        corr,
        title="Correlation matrix",
        cmap="RdBu_r",
        vmin=-1.0,
        vmax=1.0,
        center=0.0,
        fmt=".2f",
    )

    plt.tight_layout()

    safe_name = name.replace("/", "__").replace(";", "_").replace(" ", "_")
    out_path = os.path.join(outdir, f"{safe_name}.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"  Saved → {out_path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Plot covariance & correlation matrices from a ROOT file."
    )
    parser.add_argument("root_file", help="Path to the input ROOT file")
    parser.add_argument(
        "--outdir",
        default="matrix_plots",
        help="Directory for output PNG files (default: ./matrix_plots)",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.root_file):
        print(f"[ERROR] File not found: {args.root_file}")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    matrices = load_matrices(args.root_file)

    for name, cov in matrices.items():
        print(f"Processing '{name}' ({cov.shape[0]}×{cov.shape[1]}) …")
        if cov.shape[0] != cov.shape[1]:
            print(f"  [WARN] Not square — skipping correlation conversion.")
            fig, ax = plt.subplots(figsize=(7, 6))
            plot_matrix(ax, cov, title=make_label(name), cmap="viridis")
            plt.tight_layout()
            safe = name.replace("/", "__").replace(";", "_")
            fig.savefig(os.path.join(args.outdir, f"{safe}.png"), dpi=150, bbox_inches="tight")
            plt.close(fig)
        else:
            plot_pair(name, cov, args.outdir)

    print(f"\nDone! All plots saved to '{args.outdir}/'")


if __name__ == "__main__":
    main()
