# Wobble-Suppressed Calibration of Mean-Pole Polar Motion

Paper, figures, dated data snapshot and reproduction script for separating the
Chandler/annual wobble from the secular drift of the Earth's mean pole.

## Contents

- `main.tex` / `main.pdf` --- paper source and built PDF (23 pages).
- `figures/` --- nine PDF figures (trajectory, epoch rates with bootstrap CI,
  rolling rates, cumulative displacement, baseline residual vector ellipses, residual-rate trajectory,
  residual displacement path, burst half-width sweep, sliding-end-time raw rate spread).
- `data/` --- dated EOP snapshot (2026-11-24) for IERS finals, JPL EOP2 and C04,
  with SHA-256 checksums recorded in the paper.
- `scripts/polar_drift_repro.py` --- reproduces diagnostics, epoch rates with
  block-bootstrap intervals, net displacement, rolling rates, burst sensitivity,
  baseline residual vector ellipse and all figures.
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
- Raw daily-window slope 19.9 mas/yr is mostly annual leakage and phase: sliding end times over the last 12 years give ranges 2.0--22.9 (3.55 yr), 1.7--20.7 (4.12 yr), 0.5--16.9 (4.74 yr) and 0.3--11.3 (7.11 yr); short windows under about 8 years have end-time errors of order 10 mas/yr, and no length fixes that.
- Centered +/-2 yr burst 5.9 mas/yr is unstable: half-width sweep 1.0--3.0 yr gives 8.1--2.7 mas/yr; harmonic fits over 6--10 yr give 6.1--3.7; longer baselines give raw OLS 4.5 (2000--end), 5.2 (2010--end), 5.1 (2015--end), and calendar-year means 4.0--3.7 (2018--2020--2026). The data support roughly 4--5 mas/yr since about 2000, not a distinct 18 cm/yr episode.
- Middle common window 2010--2020: 4.6 mas/yr in all three products.
- Net displacement: 5.94 m (IERS finals), 7.18 m (JPL), 7.59 m (C04); recent window adds 0.74 m and is approximate.
- Baseline residual vector: raw OLS 2010--end minus Adhikari & Ivins 2016 climate+GIA budget is about (0.1, 0.8) mas/yr; component spreads ±1.4 and ±1.3 mas/yr; 95% ellipse axes 2.62 by 3.87 (finals), 2.63 by 3.79 (JPL) and 2.68 by 3.77 (C04); zero lies inside the ellipse (p≈0.95), consistent with zero to about 2.5--3 mas/yr. The benchmark overlaps 2010--2020 and is not an independent test; the unstable four-year residual is not the headline.
## Scope

Observed record calibration only. No mechanism, coupling or causal claim is made.
