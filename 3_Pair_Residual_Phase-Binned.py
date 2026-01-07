# MIT License
# Copyright (c) 2026 Christopher Dean White

"""
Generate:
1) Pair residual phase-binned two ways (±Y ns)
2) Bin counts (nbins)
3) π-inversion score (best circular phase shift + correlation), written to CSV

Input: *_pair.csv like /mnt/data/MKEA_RIO2_pair.csv

Expected columns:
  resid_<SITEA>, resid_<SITEB>
  <SITEA>_alt_deg, <SITEA>_slope_deg_per_s
  <SITEB>_alt_deg, <SITEB>_slope_deg_per_s
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    mu = np.nanmean(x)
    sd = np.nanstd(x)
    if not np.isfinite(sd) or sd == 0:
        raise ValueError("zscore: std is zero / non-finite")
    return (x - mu) / sd

def phase_from_alt_slope(alt_deg: np.ndarray, slope_deg_per_s: np.ndarray) -> np.ndarray:
    """φ = atan2(z(alt), z(slope)) mapped to [0, 2π)."""
    za = zscore(alt_deg)
    zs = zscore(slope_deg_per_s)
    phi = np.arctan2(za, zs)
    return np.mod(phi, 2 * np.pi)

def set_phase_ticks(ax):
    """Format x-axis as phase rotation with canonical π ticks."""
    ticks = [0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]
    labels = [r"$0$", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"]
    ax.set_xlim(0, 2*np.pi)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)

def peak_to_peak(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.nanmax(x) - np.nanmin(x))

def bin_means(phi: np.ndarray, y: np.ndarray, nbins: int = 24):
    """Bin y by phi across [0, 2π). Return centers, means, counts."""
    edges = np.linspace(0, 2 * np.pi, nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    means = np.full(nbins, np.nan, dtype=float)
    counts = np.zeros(nbins, dtype=int)

    idx = np.digitize(phi, edges, right=False) - 1
    idx = np.clip(idx, 0, nbins - 1)

    for b in range(nbins):
        m = idx == b
        counts[b] = int(np.sum(m))
        if counts[b] > 0:
            means[b] = float(np.nanmean(y[m]))

    return centers, means, counts

def infer_sites_from_filename(path: Path) -> tuple[str, str] | None:
    """Infer SITEA/SITEB from filenames like MKEA_RIO2_pair.csv."""
    stem = path.stem  # e.g. "MKEA_RIO2_pair"
    parts = stem.split("_")
    if len(parts) >= 3 and parts[-1].lower() == "pair":
        return parts[0], parts[1]
    return None

def nan_corr(a: np.ndarray, b: np.ndarray) -> float:
    """Pearson correlation ignoring NaNs."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3:
        return np.nan
    aa = a[m]
    bb = b[m]
    sa = np.std(aa)
    sb = np.std(bb)
    if sa == 0 or sb == 0:
        return np.nan
    return float(np.corrcoef(aa, bb)[0, 1])

def best_circular_shift_score(a: np.ndarray, b: np.ndarray, prefer_pi: bool = True):
    """
    Find best circular shift k (0..N-1) for b such that:
      - either a ~ shift(b, k)  (direct)
      - or a ~ -shift(b, k)     (inverted)
    Returns dict with best mode, shift, corr, etc.

    If prefer_pi=True, ties are broken in favor of shift closest to π (N/2).
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n = len(a)
    if len(b) != n:
        raise ValueError("a and b must have same length")

    target_k_pi = n / 2.0

    def tie_key(k: int) -> float:
        # smaller is better (distance to π-shift)
        d = abs(k - target_k_pi)
        return min(d, n - d)  # wraparound distance

    best = {
        "mode": None,          # "direct" or "inverted"
        "k_bins": None,        # integer bin shift
        "corr": -np.inf,       # maximize absolute corr, store signed corr separately
        "signed_corr": None,   # signed correlation for the best transform
        "abs_corr": None,
    }

    for k in range(n):
        b_shift = np.roll(b, k)

        c_direct = nan_corr(a, b_shift)
        c_invert = nan_corr(a, -b_shift)

        # Choose best between direct/invert at this k by abs corr
        candidates = [
            ("direct", c_direct),
            ("inverted", c_invert),
        ]

        for mode, c in candidates:
            if not np.isfinite(c):
                continue
            abs_c = abs(c)

            better = False
            if abs_c > best["corr"]:
                better = True
            elif abs_c == best["corr"] and prefer_pi and best["k_bins"] is not None:
                # tie-break: closer to π shift
                if tie_key(k) < tie_key(best["k_bins"]):
                    better = True

            if better:
                best["mode"] = mode
                best["k_bins"] = k
                best["corr"] = abs_c
                best["signed_corr"] = float(c)
                best["abs_corr"] = float(abs_c)

    if best["k_bins"] is None:
        return {
            "mode": None,
            "k_bins": None,
            "shift_rad": np.nan,
            "shift_deg": np.nan,
            "signed_corr": np.nan,
            "abs_corr": np.nan,
        }

    shift_rad = best["k_bins"] * (2 * np.pi / n)
    shift_deg = shift_rad * (180.0 / np.pi)

    return {
        "mode": best["mode"],
        "k_bins": int(best["k_bins"]),
        "shift_rad": float(shift_rad),
        "shift_deg": float(shift_deg),
        "signed_corr": float(best["signed_corr"]),
        "abs_corr": float(best["abs_corr"]),
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pair_csv", help="Path to *_pair.csv (e.g., MKEA_RIO2_pair.csv)")
    ap.add_argument("--siteA", default=None, help="Override SITEA code (e.g., MKEA)")
    ap.add_argument("--siteB", default=None, help="Override SITEB code (e.g., RIO2)")
    ap.add_argument("--nbins", type=int, default=24)
    ap.add_argument("--ylim", type=float, default=100.0, help="Y-limit in ns for residual plot")
    ap.add_argument("--outdir", default=".", help="Output directory")
    args = ap.parse_args()

    pair_path = Path(args.pair_csv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    inferred = infer_sites_from_filename(pair_path)
    siteA = args.siteA or (inferred[0] if inferred else None)
    siteB = args.siteB or (inferred[1] if inferred else None)
    if siteA is None or siteB is None:
        raise SystemExit("Could not infer site codes; pass --siteA and --siteB")

    df = pd.read_csv(pair_path)

    residA_col = f"resid_{siteA}"
    residB_col = f"resid_{siteB}"
    altA_col = f"{siteA}_alt_deg"
    slopeA_col = f"{siteA}_slope_deg_per_s"
    altB_col = f"{siteB}_alt_deg"
    slopeB_col = f"{siteB}_slope_deg_per_s"

    for col in [residA_col, residB_col, altA_col, slopeA_col, altB_col, slopeB_col]:
        if col not in df.columns:
            raise SystemExit(f"Missing required column '{col}'. Found: {list(df.columns)}")

    # φ per site
    phiA = phase_from_alt_slope(df[altA_col].to_numpy(), df[slopeA_col].to_numpy())
    phiB = phase_from_alt_slope(df[altB_col].to_numpy(), df[slopeB_col].to_numpy())

    # Pair residual (ns)
    pair_ns = (df[residA_col].to_numpy() - df[residB_col].to_numpy()) * 1e9

    # Bin two ways
    centers, mean_by_A, count_by_A = bin_means(phiA, pair_ns, nbins=args.nbins)
    _,       mean_by_B, count_by_B = bin_means(phiB, pair_ns, nbins=args.nbins)

    # Summary CSV
    summary = pd.DataFrame({
        "phi_center_rad": centers,
        f"pair_mean_by_phi_{siteA}_ns": mean_by_A,
        f"count_by_phi_{siteA}": count_by_A,
        f"pair_mean_by_phi_{siteB}_ns": mean_by_B,
        f"count_by_phi_{siteB}": count_by_B,
    })
    summary_csv = outdir / f"{siteA}_{siteB}_pair_binned_summary.csv"
    summary.to_csv(summary_csv, index=False)

    # p2p for legend
    p2p_A = peak_to_peak(mean_by_A)
    p2p_B = peak_to_peak(mean_by_B)

    # π-inversion score between the two binned curves
    # We score mean_by_A vs mean_by_B with optimal circular shift, allowing sign inversion.
    inv_score = best_circular_shift_score(mean_by_A, mean_by_B, prefer_pi=True)

    # Write score CSV (one row)
    score_csv = outdir / f"{siteA}_{siteB}_pi_inversion_score.csv"
    score_row = {
        "siteA": siteA,
        "siteB": siteB,
        "nbins": args.nbins,
        "mode": inv_score["mode"],          # "direct" or "inverted"
        "k_bins": inv_score["k_bins"],      # shift in bins
        "shift_rad": inv_score["shift_rad"],
        "shift_deg": inv_score["shift_deg"],
        "signed_corr": inv_score["signed_corr"],
        "abs_corr": inv_score["abs_corr"],
        "p2p_by_phi_siteA_ns": p2p_A,
        "p2p_by_phi_siteB_ns": p2p_B,
        "n_total_rows": int(len(df)),
    }
    pd.DataFrame([score_row]).to_csv(score_csv, index=False)

    # -------- Plot 1: residuals --------
    resid_png = outdir / f"{siteA}_{siteB}_pair_residual_binned_pm{int(args.ylim)}.png"
    plt.figure(figsize=(14, 8), dpi=120)
    plt.plot(
        centers, mean_by_A, marker="o",
        label=f"Pair ({siteA} − {siteB}) binned by φ_{siteA}  (p2p ≈ {p2p_A:.1f} ns)"
    )
    plt.plot(
        centers, mean_by_B, marker="o",
        label=f"Pair ({siteA} − {siteB}) binned by φ_{siteB}  (p2p ≈ {p2p_B:.1f} ns)"
    )

    # Annotate π-inversion score on plot
    ann = (
        f"Best shift: {inv_score['shift_deg']:.1f}° "
        f"({inv_score['k_bins']} bins), mode={inv_score['mode']}, "
        f"corr={inv_score['signed_corr']:.3f} (|corr|={inv_score['abs_corr']:.3f})"
    )

    ax = plt.gca()
    ax.text(0.02, 0.03, ann, transform=ax.transAxes, fontsize=10, alpha=0.9)

    ax.axhline(0, linewidth=1)
    ax.grid(True, alpha=0.25)
    ax.set_title(f"{siteA} − {siteB} pair residual, phase-binned two ways (±{int(args.ylim)} ns)")

    set_phase_ticks(ax)
    for x in (np.pi/2, 3*np.pi/2):
        ax.axvline(x, alpha=0.15, linewidth=1)
    ax.set_xlabel(r"Lunar-slope phase rotation $\phi$ (rad)")
    ax.set_ylabel("Pair residual (ns) (24-bin mean)")
    ax.set_ylim(-args.ylim, args.ylim)
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(resid_png)
    plt.close()

    # -------- Plot 2: bin counts --------
    counts_png = outdir / f"{siteA}_{siteB}_pair_bin_counts.png"
    plt.figure(figsize=(14, 6), dpi=120)
    plt.plot(centers, count_by_A, marker="o", label=f"Counts per bin (φ_{siteA})")
    plt.plot(centers, count_by_B, marker="o", label=f"Counts per bin (φ_{siteB})")
    plt.grid(True, alpha=0.25)
    plt.title(f"{siteA} − {siteB} bin counts ({args.nbins} bins)")
    ax = plt.gca()
    set_phase_ticks(ax)
    ax.set_xlabel(r"Lunar-slope phase rotation $\phi$ (rad)")
    plt.ylabel("Count")
    plt.legend(loc="lower center")
    plt.tight_layout()
    plt.savefig(counts_png)
    plt.close()

    print("Wrote:")
    print(f"  {summary_csv}")
    print(f"  {score_csv}")
    print(f"  {resid_png}")
    print(f"  {counts_png}")


if __name__ == "__main__":
    main()
