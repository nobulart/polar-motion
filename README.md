# Wobble-Suppressed Calibration of Mean-Pole Polar Motion

Paper, figures, dated data snapshot and reproduction script for separating the
Chandler/annual wobble from the secular drift of the Earth's mean pole.

## Contents

- `main.tex` / `main.pdf` --- paper source and built PDF (20 pages).
- `figures/` --- seven PDF figures (trajectory, epoch rates with bootstrap CI,
  rolling rates, cumulative displacement, residual rates, residual-rate trajectory,
  residual displacement path).
- `data/` --- dated EOP snapshot (2026-11-24) for IERS finals, JPL EOP2 and C04,
  with SHA-256 checksums recorded in the paper.
- `scripts/polar_drift_repro.py` --- reproduces diagnostics, epoch rates with
  block-bootstrap intervals, net displacement, rolling rates, residual subtraction
  against the known climate+GIA budget and all figures.

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
- Raw daily-window slope: 20.0 mas/yr; sequential sinusoid: 10.7 mas/yr (wobble inflation).
- Centered +/-2 yr smoothing: burst 5.9 mas/yr (18.1 cm/yr), bootstrap CI 4.3--9.8;
  interior-only variant 3.8 mas/yr.
- Middle common window 2010--2020: 4.6 mas/yr in all three products.
- Net displacement: 5.94 m (IERS finals), 7.18 m (JPL), 7.58 m (C04); recent window adds 0.74 m.
- Residual benchmark after subtracting Adhikari & Ivins 2016 climate reconstruction (4.56 ±1.26 mas/yr at 27.7°E)
  and GIA trend (3.06 mas/yr at 75°W): early-era subtraction omitted as anachronistic; burst interior-only 3.2 mas/yr
  (CI not meaningful; 2.1-year window collapses to one 730-day block); burst full-window 3.8 mas/yr (CI 1.9--7.7) is edge-biased;
  middle 1.4 mas/yr (CI 0.4--5.0) is circular. Climate-error-propagated bootstrap: only ~12% of burst replicates below 3 mas/yr,
  fewer than 1% below 1 mas/yr; ellipse axes about 3.3 by 4.5 mas/yr; residual positive and not resolved to zero.
  Combined uncertainty ±3 mas/yr or more; IERS Conventions mean-pole model comparison is the next benchmark, not done here.
## Scope

Observed record calibration only. No mechanism, coupling or causal claim is made.
