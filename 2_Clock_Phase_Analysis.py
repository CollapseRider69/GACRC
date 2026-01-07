# MIT License
# Copyright (c) 2026 Christopher Dean White

import argparse, os, zipfile, gzip, math
from datetime import datetime, timezone
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def parse_args():
    p = argparse.ArgumentParser(description="Analyze JPL FIN clock time drift vs Moon phase (Skyfield CSVs).")
    p.add_argument("--zip", required=True, help="Path to JPL.zip containing FIN *_CLK.CLK.gz files")
    p.add_argument("--stations", nargs="+", required=True, help="Station codes to extract (e.g., FAIR MAW1)")
    p.add_argument("--moon", nargs="+", required=True, help="Mappings like FAIR=/path/FAIR.csv MAW1=/path/MAW1.csv")
    p.add_argument("--out", default="out", help="Output directory")
    p.add_argument("--tol", type=int, default=60, help="As-of join tolerance in seconds for phase labels (default 60s)")
    p.add_argument("--rolling", type=str, default="3H", help="Rolling window for correlation (default 3H)")
    return p.parse_args()

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)
    return p

def parse_line_fast(ln):
    parts = ln.split()
    if len(parts) < 11:
        return None
    try:
        dt = datetime(
            int(parts[2]), int(parts[3]), int(parts[4]),
            int(parts[5]), int(parts[6]), int(float(parts[7]))
        ).replace(tzinfo=timezone.utc)
        off = float(parts[9])   # FIN "offset" = correction applied
        dr  = float(parts[10])  # FIN "drift" rate (may be unused here)
        return dt, off, dr
    except Exception:
        return None

def extract_station_series(fin_files, station):
    rows = []
    for p in fin_files:
        with gzip.open(p, 'rt', errors='ignore') as f:
            for ln in f:
                if not ln.startswith("AR "):
                    continue
                code = ln[3:8].strip()
                if code != station:
                    continue
                parsed = parse_line_fast(ln)
                if parsed is None:
                    continue
                rows.append(parsed)
    if not rows:
        return pd.DataFrame(columns=["datetime","offset","drift","tdrift"])
    df = pd.DataFrame(rows, columns=["datetime","offset","drift"]).sort_values("datetime")
    df = df.drop_duplicates(subset=["datetime"]).reset_index(drop=True)
    # Define observed time drift as the negative of the correction (offset)
    df["tdrift"] = -df["offset"]
    return df

def load_moon_csv(path):
    df = pd.read_csv(path)
    if "timestamp_utc" not in df.columns or "alt_deg" not in df.columns:
        raise ValueError(f"Moon CSV {path} must contain 'timestamp_utc' and 'alt_deg'.")
    ts = df["timestamp_utc"].astype(str).str.replace("Z","+00:00", regex=False)
    df["datetime"] = pd.to_datetime(ts, utc=True, errors="coerce")
    df = df.sort_values("datetime").dropna(subset=["datetime","alt_deg"]).reset_index(drop=True)
    # Compute slope d(alt)/dt with *aligned* time deltas.
    # Use a forward difference so the slope at time t reflects (alt[t+1] - alt[t]) / (t[t+1] - t[t]).
    dt_f = (df["datetime"].shift(-1) - df["datetime"]).dt.total_seconds()
    dt_f = dt_f.replace(0, np.nan)
    df["slope_deg_per_s"] = (df["alt_deg"].shift(-1) - df["alt_deg"]) / dt_f
    df["phase"] = np.where(df["slope_deg_per_s"] > 0, "ascending", "descending")
    return df[["datetime","alt_deg","slope_deg_per_s","phase"]].dropna()

def fisher_ci(r, n):
    if n is None or pd.isna(r) or n <= 3:
        return (np.nan, np.nan)
    z = np.arctanh(r)
    se = 1/np.sqrt(n-3)
    lo = np.tanh(z - 1.96*se)
    hi = np.tanh(z + 1.96*se)
    return lo, hi

def P2_from_alt_deg(alt_deg):
    a = np.deg2rad(alt_deg.astype(float))
    return 0.5*(3.0*np.sin(a)**2 - 1.0)

def z(x):
    x = np.asarray(x, float)
    return (x - np.nanmean(x)) / np.nanstd(x)

def attach(pair, moon_df, label, tol_seconds):
    """Attach nearest moon altitude/slope/phase to a pair DataFrame by datetime."""
    moon_sorted = moon_df.sort_values("datetime")
    pair_sorted = pair.sort_values("datetime")
    merged = pd.merge_asof(
        pair_sorted,
        moon_sorted,
        on="datetime",
        direction="nearest",
        tolerance=pd.Timedelta(seconds=tol_seconds),
    )
    return pair_sorted.assign(
        **{
            f"{label}_alt_deg": merged["alt_deg"].to_numpy(),
            f"{label}_slope_deg_per_s": merged["slope_deg_per_s"].to_numpy(),
            f"{label}_phase": merged["phase"].to_numpy(),
        }
    )

def main():
    args = parse_args()
    out = ensure_dir(args.out)

    # Unzip FIN files
    with zipfile.ZipFile(args.zip, 'r') as zf:
        zf.extractall(out)

    # Collect FIN *_CLK.CLK.gz files
    fin_files = []
    for root, dirs, files in os.walk(out):
        for f in files:
            if f.startswith("JPL0OPSFIN") and f.endswith("_CLK.CLK.gz"):
                fin_files.append(os.path.join(root, f))
    fin_files = sorted(fin_files)
    if not fin_files:
        raise SystemExit("No FIN *_CLK.CLK.gz files found after unzip.")

    # Extract series per station (includes tdrift = -offset)
    series = {}
    for st in args.stations:
        df = extract_station_series(fin_files, st)
        if df.empty:
            print(f"[WARN] No data for station {st}")
        series[st] = df
        df.to_csv(os.path.join(out, f"{st}_offsets_drifts.csv"), index=False)

    if len(args.stations) < 2:
        print("Only one station provided; skipping pair analyses.")
        return

    A_code, B_code = args.stations[0], args.stations[1]
    A, B = series[A_code].copy(), series[B_code].copy()

    # Merge pair on exact timestamps (UTC)
    pair = pd.merge(
        A[["datetime","tdrift","offset","drift"]],
        B[["datetime","tdrift","offset","drift"]],
        on="datetime",
        how="inner",
        suffixes=(f"_{A_code}", f"_{B_code}")
    ).sort_values("datetime").reset_index(drop=True)

    # Build moon_map from CLI
    moon_map = {}
    for m in args.moon:
        if "=" not in m:
            raise SystemExit(f"--moon entries must be CODE=/path/file.csv; got: {m}")
        code, path = m.split("=", 1)
        moon_map[code.strip()] = load_moon_csv(path.strip())

    # Safety: ensure Moon CSVs provided for both stations
    if A_code not in moon_map or B_code not in moon_map:
        raise SystemExit(
            f"Missing Moon CSV for {A_code} or {B_code}. Have: {list(moon_map.keys())}"
        )

    # Attach Moon altitude/phase to the merged pair
    pair = attach(pair, moon_map[A_code], A_code, args.tol)
    pair = attach(pair, moon_map[B_code], B_code, args.tol)

    # Sanity-check required columns exist and drop NaNs
    req_cols = [
        f"{A_code}_alt_deg", f"{B_code}_alt_deg",
        f"{A_code}_slope_deg_per_s", f"{B_code}_slope_deg_per_s",
    ]
    missing = [c for c in req_cols if c not in pair.columns]
    if missing:
        raise RuntimeError(
            f"Missing expected columns after moon attach: {missing}. "
            f"Stations: {[A_code, B_code]}  Moon keys: {list(moon_map.keys())}"
        )
    pair = pair.dropna(subset=req_cols).reset_index(drop=True)

    # --- split analysis (operate on time drift, not correction) ---
    if pair.empty:
        print("[WARN] No rows after moon alignment; skipping residual analysis.")
    else:
        # 1) Remove GR-like tidal term per site via P2(alt) regression on time drift
        for code in (A_code, B_code):
            y = pair[f"tdrift_{code}"].astype(float).values
            x = P2_from_alt_deg(pair[f"{code}_alt_deg"]).values
            X = np.column_stack([np.ones_like(x), x])  # intercept + P2
            beta = np.linalg.lstsq(X, y, rcond=None)[0]
            fit  = X @ beta
            pair[f"resid_{code}"] = y - fit

        # 2) Δslope (normalized)
        pair["slope_A"] = z(pair[f"{A_code}_slope_deg_per_s"].astype(float))
        pair["slope_B"] = z(pair[f"{B_code}_slope_deg_per_s"].astype(float))
        pair["delta_slope"] = z(pair["slope_A"] - pair["slope_B"])

        # 3) Rolling residual correlation (between residual time drifts)
        pair_dt = pair.set_index("datetime").sort_index()
        rwin = args.rolling
        roll_resid = pair_dt[f"resid_{A_code}"].rolling(rwin).corr(
                        pair_dt[f"resid_{B_code}"]).dropna().reset_index()
        roll_resid.columns = ["datetime", f"rolling_corr_resid_{rwin}"]

        # 4) Join Δslope and test (no extra minus; dependent var already drift)
        joined = pd.merge_asof(
            roll_resid.sort_values("datetime"),
            pair[["datetime", "delta_slope"]].sort_values("datetime"),
            on="datetime",
            direction="nearest",
            tolerance=pd.Timedelta(seconds=args.tol)
        ).dropna()

        if len(joined) >= 3:
            rname = f"rolling_corr_resid_{rwin}"
            r0 = np.corrcoef(joined[rname], joined["delta_slope"])[0,1]

            # lag scan ±12h (if enough points and regular-ish spacing)
            if len(joined) >= 10:
                dt_s = (joined["datetime"].iloc[1] - joined["datetime"].iloc[0]).total_seconds()
                if dt_s > 0:
                    max_h = 12
                    steps = max(1, int(max_h*3600/dt_s))
                    lags = np.arange(-steps, steps+1)
                    cors = []
                    for k in lags:
                        if k > 0:
                            cors.append(np.corrcoef(joined[rname][k:], joined["delta_slope"][:-k])[0,1])
                        elif k < 0:
                            cors.append(np.corrcoef(joined[rname][:k], joined["delta_slope"][-k:])[0,1])
                        else:
                            cors.append(r0)
                    best_idx = int(np.nanargmax(cors))
                    best_r   = cors[best_idx]
                    best_lag_h = lags[best_idx]*dt_s/3600.0
                else:
                    best_r = np.nan; best_lag_h = np.nan
            else:
                best_r = np.nan; best_lag_h = np.nan

            asc_mask  = joined["delta_slope"] > 0
            desc_mask = joined["delta_slope"] < 0
            r_asc  = np.corrcoef(joined.loc[asc_mask,  rname], joined.loc[asc_mask,  "delta_slope"])[0,1] if asc_mask.any() else np.nan
            r_desc = np.corrcoef(joined.loc[desc_mask, rname], joined.loc[desc_mask, "delta_slope"])[0,1] if desc_mask.any() else np.nan
            bias_amp = 0.5*(r_asc - r_desc) if np.isfinite(r_asc) and np.isfinite(r_desc) else np.nan

            print("\n=== gPCT residual tests (time drift, after GR-like removal) ===")
            print(f"Pearson r(resid-corr vs Δslope) = {r0:.3f}")
            print(f"Best lag alignment: {best_lag_h:+.2f} h (r={best_r:.3f})")
            print(f"Ascending r : {r_asc:.3f}  |  Descending r : {r_desc:.3f}  |  Bias amplitude : {bias_amp:.3f}")

            # Save joined + residual-rolling-corr CSVs for inspection
            joined_out = os.path.join(out, f"{A_code}_{B_code}_joined_RESID_TDRIFT_{rwin}.csv")
            roll_resid_out = os.path.join(out, f"{A_code}_{B_code}_rolling_corr_RESID_TDRIFT_{rwin}.csv")
            joined.to_csv(joined_out, index=False)
            roll_resid.to_csv(roll_resid_out, index=False)

            # Phase-fold plot of residual corr vs Δslope (tidal semidiurnal ~12.42h)
            T = 12.42*3600.0
            t0 = joined["datetime"].iloc[0]
            phase = ((joined["datetime"] - t0).dt.total_seconds() % T)/T * 2*np.pi
            bins = np.linspace(0, 2*np.pi, 73)
            cent = 0.5*(bins[1:] + bins[:-1])
            b_r, b_ds = [], []
            for i in range(len(bins)-1):
                m = (phase >= bins[i]) & (phase < bins[i+1])
                if m.any():
                    b_r.append(np.nanmean(joined.loc[m, rname]))
                    b_ds.append(np.nanmean(joined.loc[m, "delta_slope"]))
                else:
                    b_r.append(np.nan); b_ds.append(np.nan)

            plt.figure(figsize=(8,4))
            plt.plot(np.degrees(cent), b_r, label="Residual rolling corr (tdrift)", lw=1)
            plt.plot(np.degrees(cent), np.array(b_ds), label="Δslope (norm)", lw=1)
            plt.axhline(0, color="k", lw=0.5)
            plt.xlabel("Lunar phase (deg, 0–360 ≈ 12.42h)")
            plt.ylabel("Normalized")
            plt.title(f"{A_code}-{B_code}: Phase-folded residual corr vs Δslope (tdrift)")
            plt.legend(fontsize=8)
            plt.tight_layout()
            plt.savefig(os.path.join(out, f"{A_code}_{B_code}_phase_folded_RESID_TDRIFT_{rwin}.png"), dpi=200)
            plt.close()
        else:
            print("[INFO] Not enough joined points for residual tests; skipping.")

    # --- Save the merged pair (with moon columns & tdrift) ---
    pair_path = os.path.join(out, f"{A_code}_{B_code}_pair.csv")
    pair.to_csv(pair_path, index=False)

    # Quadrant stats on time drift
    rows = []
    combos = [("ascending","ascending"),("ascending","descending"),
              ("descending","ascending"),("descending","descending")]
    for pa, pb in combos:
        sub = pair[(pair[f'{A_code}_phase']==pa) & (pair[f'{B_code}_phase']==pb)]
        n = len(sub)
        r = sub[f'tdrift_{A_code}'].corr(sub[f'tdrift_{B_code}']) if n >= 3 else np.nan
        lo, hi = fisher_ci(r, n) if n >= 4 else (np.nan, np.nan)
        rows.append({
            f'{A_code}_phase': pa, f'{B_code}_phase': pb,
            "N": n, "r_tdrift": r, "95%CI_low": lo, "95%CI_high": hi
        })
    quad = pd.DataFrame(rows)
    quad_path = os.path.join(out, f"{A_code}_{B_code}_quadrants_TDRIFT.csv")
    quad.to_csv(quad_path, index=False)

    # Rolling correlation on time drift
    pair_dt = pair.set_index("datetime").sort_index()
    rwin = args.rolling
    roll = pair_dt[f'tdrift_{A_code}'].rolling(rwin).corr(
               pair_dt[f'tdrift_{B_code}']).dropna().reset_index()
    roll.columns = ["datetime", f"rolling_corr_TDRIFT_{rwin}"]
    roll_path = os.path.join(out, f"{A_code}_{B_code}_rolling_corr_TDRIFT_{rwin}.csv")
    roll.to_csv(roll_path, index=False)

    plt.figure(figsize=(12,5))
    plt.plot(roll["datetime"], roll[f"rolling_corr_TDRIFT_{rwin}"])
    plt.xlabel("Time (UTC)")
    plt.ylabel(f"Rolling correlation on time drift ({rwin})")
    plt.title(f"{A_code} vs {B_code} — Rolling correlation (time drift, {rwin})")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(out, f"{A_code}_{B_code}_rolling_corr_TDRIFT_{rwin}.png"))
    plt.close()

    # Daily by phase (using A station’s phase labels)
    pair_dt["date"] = pair_dt.index.date
    daily_rows = []
    for d, sub in pair_dt.groupby("date"):
        for phase in ["ascending","descending"]:
            chunk = sub[sub[f"{A_code}_phase"] == phase]
            r = chunk[f'tdrift_{A_code}'].corr(chunk[f'tdrift_{B_code}']) if len(chunk) >= 5 else np.nan
            daily_rows.append({"date": str(d), f"{A_code}_phase": phase, "r_tdrift": r, "N": len(chunk)})
    daily = pd.DataFrame(daily_rows)
    daily.to_csv(os.path.join(out, f"{A_code}_{B_code}_daily_by_{A_code}_phase_TDRIFT.csv"), index=False)

    print("DONE. Outputs in:", out)

if __name__ == "__main__":
    main()
