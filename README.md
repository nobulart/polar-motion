# Wobble-Suppressed Calibration of Mean-Pole Polar Motion

Paper, figures, dated data snapshot and reproduction script for separating the
Chandler/annual wobble from the secular drift of the Earth's mean pole.

## Contents

- `main.tex` / `main.pdf` --- paper source and built PDF (22 pages).
- `figures/` --- nine PDF figures (trajectory, epoch rates with bootstrap CI,
  rolling rates, cumulative displacement, baseline residual rates, residual-rate trajectory,
  residual displacement path, burst half-width sweep, Chandler-window raw rates).
- `data/` --- dated EOP snapshot (2026-11-24) for IERS finals, JPL EOP2 and C04,
  with SHA-256 checksums recorded in the paper.
- `scripts/polar_drift_repro.py` --- reproduces diagnostics, epoch rates with
  block-bootstrap intervals, net displacement, rolling rates, burst sensitivity,
  baseline residual subtraction and all figures.
## Build

```
cd latex/papers/polar_motion_wobble_calibration   # if inside the vault tree
tectonic main.tex
```

## Reproduce

```
/Users/craig/.venvs/figgen/bin/python3 scripts/polar_drift_repro.py
```

The script loads the archived snapshot by default and verifies its checksums;
`--live` fetches the DRIFT endpoint instead. Figures are written to `figures/`.

## Result summary

- Corrected one-basis co-fit (periods in years): 3.2 mas/yr.
- Raw daily-window slope 20.0 mas/yr is Chandler leakage: a 4.12-year window sits near a half-integer Chandler cycle; whole multiples give 20.3 (2.37 yr), 5.7 (3.55 yr) and 4.9 (4.74 yr).
- Centered +/-2 yr burst 5.9 mas/yr is unstable: half-width sweep 1.0--3.0 yr gives 8.1--2.7 mas/yr; harmonic fit over the last 4.1 yr gives 8.7; longer baselines give raw OLS 4.5 (2000--end), 5.2 (2010--end), and calendar-year means 3.7--4.7 (2020--2026). The data support roughly 4--5 mas/yr since about 2000, not a distinct 18 cm/yr episode.
- Middle common window 2010--2020: 4.6 mas/yr in all three products.
- Net displacement: 5.94 m (IERS finals), 7.18 m (JPL), 7.59 m (C04); recent window adds 0.74 m and is approximate.
- Baseline residual: raw OLS 2010--end minus Adhikari & Ivins 2016 climate+GIA budget leaves 0.8 mas/yr (CI 0.35--4.08 finals; 0.30--4.14 JPL; 0.28--4.03 C04); about 22% of replicates below 1 mas/yr, consistent with zero to about 1 mas/yr. The benchmark overlaps 2010--2020, so this is not an independent test; the unstable four-year residual is not the headline.
## Scope

Observed record calibration only. No mechanism, coupling or causal claim is made.
