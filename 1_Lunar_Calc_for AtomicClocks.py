# MIT License
# Copyright (c) 2026 Christopher Dean White

import csv
from skyfield.api import load, wgs84, utc
from datetime import datetime, timedelta

def parse_dt(s: str) -> datetime:
    """Parse ISO date or 'YYYY-MM-DD HH:MM' and default to UTC if naive."""
    dt = datetime.fromisoformat(s)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=utc)
    return dt

def main():
    # --- User input ---
    lat = float(input("Latitude (deg): ").strip())
    lon = float(input("Longitude (deg): ").strip())

    start_str = input("Start datetime (YYYY-MM-DD or YYYY-MM-DD HH:MM): ").strip()
    end_str   = input("End datetime   (YYYY-MM-DD or YYYY-MM-DD HH:MM): ").strip()

    start_dt = parse_dt(start_str)
    end_dt   = parse_dt(end_str)

    step_seconds = int(input("Step in seconds (default 30): ") or 30)

    # --- Skyfield setup ---
    ts = load.timescale()
    eph = load('de421.bsp')

    earth = eph['earth']
    moon  = eph['moon']

    observer = earth + wgs84.latlon(latitude_degrees=lat, longitude_degrees=lon)

    # --- CSV output file ---
    outname = "moon_altitude.csv"
    print(f"\nWriting CSV → {outname}")

    with open(outname, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["timestamp_utc", "alt_deg", "az_deg", "distance_km"])

        # --- Loop through time ---
        current = start_dt
        count = 0

        while current <= end_dt:
            t = ts.utc(current)
            apparent = observer.at(t).observe(moon).apparent()
            alt, az, distance = apparent.altaz()

            writer.writerow([
                current.isoformat(),
                f"{alt.degrees:.6f}",
                f"{az.degrees:.6f}",
                f"{distance.km:.3f}",
            ])

            count += 1
            current += timedelta(seconds=step_seconds)

    print(f"✅ Done. {count} rows written.")

if __name__ == "__main__":
    main()
