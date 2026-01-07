# GACRC
Phase-dependent Structure in Global Atomic Clock Residual Correlations

## License

The code in this repository is released under the MIT License.

The associated paper,
"Phase-dependent Structure in Global Atomic Clock Residual Correlations",
is licensed separately under CC BY 4.0.

# Phase-dependent Structure in Global Atomic Clock Residual Correlations (GACRC)

This repository contains analysis code used to generate the results and figures
in the paper:

**“Phase-dependent Structure in Global Atomic Clock Residual Correlations”**  
Christopher Dean White (2026), CC BY 4.0

The code implements a reproducible pipeline for analyzing inter-station atomic
clock residual correlations conditioned on lunar geometry (altitude, slope,
and phase). The repository is intentionally limited to analysis scripts and does
not redistribute raw GNSS clock data.

---

## Contents

- Lunar geometry calculation using Skyfield  
- Extraction and alignment of JPL final GNSS clock residuals  
- Phase- and slope-conditioned residual analysis  
- Quadrant-partitioned correlation statistics  
- Phase-folded rolling correlation analysis  
- Phase-binned residual plots and summary metrics  

---

## Data Availability

This repository **does not include raw GNSS clock data**.

All analyses are performed on publicly available Jet Propulsion Laboratory (JPL)
final GNSS RINEX clock products (CLK files), which must be obtained directly from
the authoritative source.

Required file naming pattern:

JPL0OPSFIN_YYYYDDD0000_01D_30S_CLK.CLK.gz

Data source:

https://sideshow.jpl.nasa.gov/pub/JPL_GNSS_Products/

Raw data are not redistributed here to ensure all analyses are performed on
canonical, unmodified source files.

---

## Requirements

- Python 3.x  
- numpy  
- pandas  
- matplotlib  
- skyfield  

---

## Data Preparation

1. Download the required daily JPL final GNSS clock files (`*.CLK.gz`) covering
   the analysis interval from the JPL GNSS Products repository.

2. Place all downloaded `.CLK.gz` files into a single directory.

3. Create a zip archive containing the clock files, for example:

   zip clocks.zip *.CLK.gz

4. Retain the original JPL filenames inside the archive.  
   The analysis code expects a single zip archive containing the raw clock files.

---

## Analysis Pipeline

The analysis is performed in three stages.

---

### 1. Lunar Geometry Calculation (per station)

**Script:**  
`1_Lunar_Calc_for AtomicClocks.py`

This interactive script computes lunar geometry for a given station using
Skyfield.

You will be prompted for:
- Station latitude and longitude  
- Start and end timestamps (UTC)  
- Time step in seconds (use 30 s to match JPL clock cadence)  

**Output:**
- `moon_altitude.csv` containing:
  - `timestamp_utc`
  - `alt_deg`
  - `az_deg`
  - `distance_km`

Run this script once per station.  
Rename each output file to the corresponding station code, for example:

FAIR.csv  
MAW1.csv  
MKEA.csv  

---

### 2. Pairwise Clock Phase Analysis

**Script:**  
`2_Clock_Phase_Analysis.py`

This is the main analysis script. It extracts clock residuals for a pair of
stations, aligns them in time, attaches lunar geometry, removes an
altitude-dependent term, and computes multiple correlation products.

**Required arguments:**
- `--zip` : path to the zip archive containing JPL CLK files  
- `--stations` : two station codes (e.g. `FAIR MAW1`)  
- `--moon` : mapping of station codes to lunar CSV files  
- `--out` : output directory  

**Optional arguments:**
- `--tol` : timestamp matching tolerance in seconds (default: 60)  
- `--rolling` : rolling correlation window (default: `3H`)  

**Example:**

python 2_Clock_Phase_Analysis.py  
--zip clocks.zip  
--stations FAIR MAW1  
--moon FAIR=FAIR.csv MAW1=MAW1.csv  
--out out/FAIR_MAW1  
--tol 60  
--rolling 3H  

**Key operations:**
- Defines observed time drift as `tdrift = -offset`  
- Aligns station time series on UTC timestamps  
- Attaches lunar altitude, slope, and ascending/descending labels  
- Removes a per-station altitude-dependent (P2) term  
- Computes:
  - Quadrant-partitioned correlations with 95% confidence intervals  
  - Phase-conditioned daily correlations  
  - Phase-folded rolling residual correlations  

**Outputs include:**
- `{A}_{B}_pair.csv`  
- `{A}_{B}_quadrants_TDRIFT.csv`  
- `{A}_{B}_rolling_corr_TDRIFT_3H.csv`  
- `{A}_{B}_rolling_corr_TDRIFT_3H.png`  
- `{A}_{B}_daily_by_{A}_phase_TDRIFT.csv`  

---

### 3. Phase-binned Residual Analysis

**Script:**  
`3_Pair_Residual_Phase-Binned.py`

This script consumes the `{A}_{B}_pair.csv` file generated in Step 2 and produces
phase-binned residual statistics and plots.

Phase is defined using a combined altitude–slope representation:

phi = atan2(z(altitude), z(slope)) ∈ [0, 2π)

**Example:**

python 3_Pair_Residual_Phase-Binned.py  
--in out/FAIR_MAW1/FAIR_MAW1_pair.csv  
--out out/FAIR_MAW1/  

**Outputs include:**
- Phase-binned residual plots  
- Bin counts  
- Summary metrics (including π-inversion score)  

---

## License

- **Code:** MIT License  
  (Provided “as is”, without warranty or liability.)

- **Paper:** CC BY 4.0  
  Licensed separately from the code.

See the `LICENSE` file for details.

---

## Notes

This repository is intended for reproducibility and inspection of the analysis
pipeline only. No causal or physical interpretation is implied by the code
itself; interpretation is discussed in the accompanying paper.
