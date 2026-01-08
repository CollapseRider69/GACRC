
"""
ld_plots.py

Compute a data-based "Lagrangian descriptor" style finite-time trajectory integral
from atomic clock residuals and plot against lunar altitude.

M(t0) = ∫_{t0-τ}^{t0+τ} |dr/dt| dt

- dr/dt via numpy.gradient using median dt
- finite-time integral via trapezoid rule over a symmetric window
- no lunar variables used in constructing M(t0); altitude only used for plotting
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------
# Core LD computation
# -----------------------------

def compute_timebase_seconds(df: pd.DataFrame, time_col: str = "datetime") -> tuple[np.ndarray, float]:
    """Return time array in seconds from start and median dt (seconds)."""
    dt_series = pd.to_datetime(df[time_col])
    t = (dt_series - dt_series.iloc[0]).dt.total_seconds().to_numpy()
    diffs = np.diff(t)
    # Guard: handle any duplicates / weirdness
    diffs = diffs[np.isfinite(diffs) & (diffs > 0)]
    if len(diffs) == 0:
        raise ValueError("Could not compute dt from timestamps (need increasing datetime).")
    dt = float(np.median(diffs))
    return t, dt


def data_ld(signal: np.ndarray, t: np.ndarray, dt: float, tau_seconds: float) -> np.ndarray:
    """
    Compute M(t0)=∫|d(signal)/dt| dt over a symmetric ±tau window.
    Edges are NaN where the window cannot be centered.
    """
    if len(signal) != len(t):
        raise ValueError("signal and t must have same length")

    dsignal_dt = np.gradient(signal, dt)

    half = int(round(tau_seconds / dt))
    if half < 1:
        raise ValueError("tau_seconds too small relative to dt")

    M = np.full(len(signal), np.nan, dtype=float)

    # Trapezoid integral of |ds/dt| over [i-half, i+half)
    for i in range(half, len(signal) - half):
        sl = slice(i - half, i + half)
        M[i] = np.trapz(np.abs(dsignal_dt[sl]), t[sl])

    return M


# -----------------------------
# Plotting helpers
# -----------------------------

def plot_ld_vs_alt(
    df: pd.DataFrame,
    resid_col: str,
    alt_col: str,
    tau_hours: float = 3.0,
    mirror_alt: bool = False,
    logy: bool = False,
    title: str | None = None,
    outpath: str | None = None,
):
    df = df.copy()
    df["datetime"] = pd.to_datetime(df["datetime"])
    df = df.sort_values("datetime").reset_index(drop=True)

    t, dt = compute_timebase_seconds(df, "datetime")

    r = df[resid_col].to_numpy(dtype=float)
    alt = df[alt_col].to_numpy(dtype=float)
    if mirror_alt:
        alt = -alt

    M = data_ld(r, t, dt, tau_seconds=tau_hours * 3600.0)

    plt.figure(figsize=(10, 4.2))
    plt.scatter(alt, M, s=2, alpha=0.6)
    plt.axvline(0, linestyle="--", linewidth=1)

    plt.xlabel(f"{'Mirrored ' if mirror_alt else ''}Lunar altitude (deg)")
    plt.ylabel(rf"$M=\int|\dot r|dt$ ({'log scale' if logy else 'linear'}), $\tau={tau_hours:g}$ h")
    if logy:
        plt.yscale("log")

    if title is None:
        title = f"Data-based Lagrangian Descriptor for {resid_col} vs {alt_col}"
    plt.title(title)

    plt.tight_layout()
    if outpath:
        plt.savefig(outpath, dpi=200)
        plt.close()
    else:
        plt.show()


def plot_stacked_two_panel(
    top_csv: str,
    top_resid: str,
    top_alt: str,
    bottom_csv: str,
    bottom_resid: str,
    bottom_alt: str,
    tau_hours: float = 3.0,
    top_mirror: bool = False,
    bottom_mirror: bool = False,
    bottom_logy: bool = True,
    title: str = "Finite-time trajectory integral (LD view) across stations",
    outpath: str = "ld_two_panel_stacked.png",
):
    # Load + prepare top
    df1 = pd.read_csv(top_csv)
    df1["datetime"] = pd.to_datetime(df1["datetime"])
    df1 = df1.sort_values("datetime").reset_index(drop=True)
    t1, dt1 = compute_timebase_seconds(df1)

    r1 = df1[top_resid].to_numpy(dtype=float)
    x1 = df1[top_alt].to_numpy(dtype=float)
    if top_mirror:
        x1 = -x1
    M1 = data_ld(r1, t1, dt1, tau_seconds=tau_hours * 3600.0)

    # Load + prepare bottom
    df2 = pd.read_csv(bottom_csv)
    df2["datetime"] = pd.to_datetime(df2["datetime"])
    df2 = df2.sort_values("datetime").reset_index(drop=True)
    t2, dt2 = compute_timebase_seconds(df2)

    r2 = df2[bottom_resid].to_numpy(dtype=float)
    x2 = df2[bottom_alt].to_numpy(dtype=float)
    if bottom_mirror:
        x2 = -x2
    M2 = data_ld(r2, t2, dt2, tau_seconds=tau_hours * 3600.0)

    # Plot stacked
    plt.figure(figsize=(10, 8.6))

    ax1 = plt.subplot(2, 1, 1)
    ax1.scatter(x1, M1, s=2, alpha=0.6)
    ax1.axvline(0, linestyle="--", linewidth=1)
    ax1.set_xlabel(f"{'Mirrored ' if top_mirror else ''}lunar altitude (deg)")
    ax1.set_ylabel(rf"$M=\int|\dot r|dt$ (linear), $\tau={tau_hours:g}$ h")
    ax1.set_title(f"Top: {top_resid} vs {top_al_
