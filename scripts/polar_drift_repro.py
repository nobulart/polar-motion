#!/usr/bin/env python3
"""Reproduce every claim in the DRIFT polar-motion calibration paper, with assertions.

Default: loads the dated snapshot in data/ (2026-11-24). Use --live to fetch
the DRIFT endpoint instead. Prints diagnostics, epoch rates with block-bootstrap
confidence intervals, net displacement and rolling rates; writes four PDF figures.

Run with the plotting interpreter:
    /Users/craig/.venvs/figgen/bin/python3 scripts/polar_drift_repro.py
"""
import argparse, hashlib, json, math
from datetime import date
import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt

API = "http://mac-studio:3111/api/eop?dataset="
R_KM = 6371.0
MAS_PER_RAD = 206264806.25
HALF = 2.0                                   # +/-2 yr => 4-year boxcar width
PERIOD_YEARS = [1.0, 0.5, 432.9 / 365.2425]   # annual, semiannual, Chandler
# Known contributions (Adhikari & Ivins 2016): climate reconstruction 4.56 +/- 1.26
# mas/yr at 27.7 +/- 7.6 deg E; GIA long-term 0.85 deg/Myr = 3.06 mas/yr at 75 deg W.
CLIMATE_MAG, CLIMATE_DIR, CLIMATE_UNC = 4.56, 27.7, 1.26
GIA_MAG, GIA_DIR = 3.06, 75.0
FIGDIR = "figures"
DATADIR = "data"
SNAPSHOT = {"finals": "eop_finals_2026-11-24.json",
            "jpl": "eop_jpl_2026-11-24.json",
            "c04": "eop_c04_2026-11-24.json"}
EXPECTED = {"finals": "b6baac64fd3472af806f91c1c67a834edee4e1a1ee96357515e0e86f8d2aa8c8",
            "jpl": "05c4546fec4f03dea8661199f8cd21aad81b4b7bc3acbd1f5167431f1edbbace",
            "c04": "f9e834b008fb696152f29d9d530332ebab5085e446d6ff954b4f51092db84b44"}
mas2cm = lambda m: R_KM * 1e5 * (m / MAS_PER_RAD)

def fetch(dataset):
    import urllib.request
    with urllib.request.urlopen(API + dataset, timeout=120) as r:
        return json.load(r)

def load_snapshot(dataset):
    path = DATADIR + "/" + SNAPSHOT[dataset]
    raw = open(path, "rb").read()
    sha = hashlib.sha256(raw).hexdigest()
    if sha != EXPECTED[dataset]:
        raise ValueError("checksum mismatch for " + dataset + ": " + sha)
    print(" ", dataset, "sha256 verified", sha)
    return json.loads(raw.decode())

def to_array(data):
    out = []
    for rec in data:
        y, m, d = map(int, rec["t"].split("-"))
        ty = y + (date(y, m, d) - date(y, 1, 1)).days / 365.0
        if math.isfinite(rec["xp"]) and math.isfinite(rec["yp"]):
            out.append((ty, rec["xp"], rec["yp"], y))
    return np.array(sorted(out))

def smooth(t, x, half=HALF):
    """Centered moving average over +/-half years; truncated at series edges."""
    out = np.full_like(x, np.nan)
    for i in range(len(t)):
        m = (t >= t[i] - half) & (t <= t[i] + half)
        if m.sum() > 0:
            out[i] = x[m].mean()
    return out

def slope(t, xx, yy, a, b):
    m = (t >= a) & (t < b)
    tt = t[m]
    A = np.vstack([np.ones_like(tt), tt]).T
    cx, *_ = np.linalg.lstsq(A, xx[m], rcond=None)
    cy, *_ = np.linalg.lstsq(A, yy[m], rcond=None)
    dx, dy = cx[1] * 1000, cy[1] * 1000
    return dx, dy, math.hypot(dx, dy), math.degrees(math.atan2(-dy, dx))

def one_basis(t, xx, yy):
    """Co-fit with periods converted to YEARS (days-as-years is a unit bug)."""
    B = [np.ones_like(t), t]
    for P in PERIOD_YEARS:
        w = 2 * math.pi / P
        B += [np.sin(w * t), np.cos(w * t)]
    A = np.vstack(B).T
    cx, *_ = np.linalg.lstsq(A, xx, rcond=None)
    cy, *_ = np.linalg.lstsq(A, yy, rcond=None)
    return math.hypot(cx[1] * 1000, cy[1] * 1000)

def sequential(t, xx, yy):
    """Linear trend first, then sinusoids on residual; corrected periods."""
    A0 = np.vstack([np.ones_like(t), t]).T
    cx0, *_ = np.linalg.lstsq(A0, xx, rcond=None)
    cy0, *_ = np.linalg.lstsq(A0, yy, rcond=None)
    rx = xx - A0 @ cx0
    ry = yy - A0 @ cy0
    B = []
    for P in PERIOD_YEARS:
        w = 2 * math.pi / P
        B += [np.sin(w * t), np.cos(w * t)]
    Pm = np.vstack(B).T
    cpx, *_ = np.linalg.lstsq(Pm, rx, rcond=None)
    cpy, *_ = np.linalg.lstsq(Pm, ry, rcond=None)
    return xx - Pm @ cpx, yy - Pm @ cpy

def boot(t, xx, yy, a, b, nboot=3000, block=730, seed=0):
    """Block bootstrap with 2-year blocks (~730 samples), seeded per call."""
    m = (t >= a) & (t < b)
    tt, ax, ay = t[m], xx[m], yy[m]
    L = len(tt)
    rng = np.random.default_rng(seed)
    rates = []
    for _ in range(nboot):
        starts = rng.integers(0, max(1, L - block), size=max(1, L // block))
        sel = []
        for s in starts:
            sel.extend(range(s, min(s + block, L)))
        if len(sel) < 60:
            continue
        A = np.vstack([np.ones(len(sel)), tt[sel]]).T
        cx, *_ = np.linalg.lstsq(A, ax[sel], rcond=None)
        cy, *_ = np.linalg.lstsq(A, ay[sel], rcond=None)
        rates.append(math.hypot(cx[1] * 1000, cy[1] * 1000))
    return np.array(rates)

def boot_vec(t, xx, yy, a, b, nboot=3000, block=730, seed=0):
    """Block bootstrap returning bx and by separately; seed per call."""
    m = (t >= a) & (t < b)
    tt, ax, ay = t[m], xx[m], yy[m]
    L = len(tt)
    rng = np.random.default_rng(seed)
    bx, by = [], []
    for _ in range(nboot):
        starts = rng.integers(0, max(1, L - block), size=max(1, L // block))
        sel = []
        for s in starts:
            sel.extend(range(s, min(s + block, L)))
        if len(sel) < 60:
            continue
        A = np.vstack([np.ones(len(sel)), tt[sel]]).T
        cx, *_ = np.linalg.lstsq(A, ax[sel], rcond=None)
        cy, *_ = np.linalg.lstsq(A, ay[sel], rcond=None)
        bx.append(cx[1] * 1000); by.append(cy[1] * 1000)
    return np.array(bx), np.array(by)

def residual_vector(dx, dy):
    """Subtract known climate+GIA vectors; frame y toward 90 W (west => +dy)."""
    cx = CLIMATE_MAG * math.cos(math.radians(CLIMATE_DIR))
    cy = -CLIMATE_MAG * math.sin(math.radians(CLIMATE_DIR))
    gx = GIA_MAG * math.cos(math.radians(GIA_DIR))
    gy = GIA_MAG * math.sin(math.radians(GIA_DIR))
    rx, ry = dx - cx - gx, dy - cy - gy
    return rx, ry, math.hypot(rx, ry), math.degrees(math.atan2(-ry, rx))

def known_components():
    cx = CLIMATE_MAG * math.cos(math.radians(CLIMATE_DIR))
    cy = -CLIMATE_MAG * math.sin(math.radians(CLIMATE_DIR))
    gx = GIA_MAG * math.cos(math.radians(GIA_DIR))
    gy = GIA_MAG * math.sin(math.radians(GIA_DIR))
    return cx + gx, cy + gy

def yearly_means(years, sx, sy):
    pos = []
    for y in np.unique(years):
        m = years == y
        if m.sum() > 200:
            pos.append((y, sx[m].mean(), sy[m].mean()))
    return np.array(pos)

def rolling(t, xx, yy, win=5.0):
    out = []
    for y in np.arange(t[0] + win, t[-1]):
        m = (t >= y - win) & (t < y)
        if m.sum() < 300:
            continue
        dx, dy, rm, lon = slope(t, xx, yy, y - win, y)
        out.append((y, rm))
    return np.array(out)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--live", action="store_true")
    args = ap.parse_args()

    print("=== DATA ===")
    store = {}
    for ds in SNAPSHOT:
        data = fetch(ds) if args.live else load_snapshot(ds)
        arr = to_array(data)
        t, xp, yp, years = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3].astype(int)
        sx, sy = smooth(t, xp), smooth(t, yp)
        store[ds] = dict(t=t, xp=xp, yp=yp, sx=sx, sy=sy, years=years)
        print(ds, "records", len(arr), "span", round(t[0], 2), "->", round(t[-1], 2))

    print("\n=== DIAGNOSTICS (IERS finals, burst window) ===")
    s = store["finals"]
    t, te = s["t"], s["t"][-1]
    raw = slope(t, s["xp"], s["yp"], te - 4.12, te)[2]
    cofit = one_basis(t, s["xp"], s["yp"])   # whole-record co-fit
    sx_seq, sy_seq = sequential(t, s["xp"], s["yp"])
    seq = slope(t, sx_seq, sy_seq, te - 4.12, te)[2]
    sm = slope(t, s["sx"], s["sy"], te - 4.12, te)[2]
    interior = slope(t, s["sx"], s["sy"], te - 4.12, te - 2.0)[2]
    print("raw burst slope          ", round(raw, 1), "mas/yr")
    print("corrected co-fit (whole record)", round(cofit, 1), "mas/yr")
    print("sequential sinusoid      ", round(seq, 1), "mas/yr")
    print("smoothed burst           ", round(sm, 1), "mas/yr")
    print("smoothed interior (excl last 2 yr)", round(interior, 1), "mas/yr")
    assert 18.0 <= raw <= 22.0, "raw burst example out of scatter range"
    assert abs(cofit - 3.2) <= 0.1, "whole-record co-fit mismatch"
    assert abs(seq - 10.7) <= 0.5, "sequential sinusoid mismatch"
    assert abs(sm - 5.9) <= 0.1, "smoothed burst mismatch"
    assert abs(interior - 3.8) <= 0.3, "interior variant mismatch"
    print("assertions passed: raw, co-fit, sequential, smoothed, interior")

    print("\n=== EPOCH RATES + BOOTSTRAP CI ===")
    for ds in SNAPSHOT:
        s = store[ds]
        t, te, t0 = s["t"], s["t"][-1], s["t"][0]
        windows = [("early", t0, t0 + 33),
                   ("middle 2010-2020", 2010, 2020),
                   ("burst", te - 4.12, te)]
        for lbl, a, b in windows:
            dx, dy, rm, lon = slope(t, s["sx"], s["sy"], a, b)
            seed = int(hashlib.sha256((ds + ":" + lbl).encode()).hexdigest()[:8], 16)
            bv = boot(t, s["sx"], s["sy"], a, b, seed=seed)
            lo, hi = np.percentile(bv, 2.5), np.percentile(bv, 97.5)
            print(ds, lbl, "bx", round(dx, 1), "by", round(dy, 1),
                  "rate", round(rm, 1), "mas/yr",
                  "cm/yr", round(mas2cm(rm), 1),
                  "CI", round(lo, 1), "-", round(hi, 1),
                  "az", round(lon, 1))

    print("\n=== RESIDUAL AFTER KNOWN SUBTRACTION ===")
    print("known climate vector 4.56 mas/yr at 27.7 deg E (+/-1.26)")
    print("known GIA vector 3.06 mas/yr at 75 deg W")
    print("early era subtraction omitted: modern budget vs 1973-2006 is not an unexplained-motion test")
    for ds in SNAPSHOT:
        s = store[ds]
        t, te = s["t"], s["t"][-1]
        windows = [("middle 2010-2020", 2010, 2020),
                   ("burst full", te - 4.12, te),
                   ("burst interior", te - 4.12, te - 2.0)]
        for lbl, a, b in windows:
            dx, dy, rm, lon = slope(t, s["sx"], s["sy"], a, b)
            rx, ry, rmag, raz = residual_vector(dx, dy)
            seed = int(hashlib.sha256((ds + ":residual:" + lbl).encode()).hexdigest()[:8], 16)
            bx, by = boot_vec(t, s["sx"], s["sy"], a, b, seed=seed)
            rb = np.array([residual_vector(x, y)[2] for x, y in zip(bx, by)])
            lo, hi = np.percentile(rb, 2.5), np.percentile(rb, 97.5)
            print(ds, lbl, "observed", round(rm, 1), "residual", round(rmag, 1),
                  "mas/yr", "CI", round(lo, 1), "-", round(hi, 1),
                  "az", round(raz, 1))
    print("combined uncertainty about +/-3 mas/yr or more; burst residual not distinguishable from zero")
    s = store["finals"]
    dx, dy, _, _ = slope(s["t"], s["sx"], s["sy"], s["t"][-1] - 4.12, s["t"][-1])
    assert abs(residual_vector(dx, dy)[2] - 3.82) <= 0.1, "computed burst residual mismatch"

    print("\n=== NET DISPLACEMENT ===")
    for ds in SNAPSHOT:
        s = store[ds]
        pos = yearly_means(s["years"], s["sx"], s["sy"])
        net = np.array([math.hypot(p[1] - pos[0, 1],
                                    p[2] - pos[0, 2]) * 1000 for p in pos])
        m = mas2cm(net[-1]) / 100
        seg = pos[(pos[:, 0] >= 2022)]
        burst_m = math.hypot(seg[-1, 1] - seg[0, 1],
                             seg[-1, 2] - seg[0, 2]) * 1000
        print(ds, "total", round(m, 2), "m", round(m * 3.281, 1), "ft",
              "| 2022-2026 adds", round(mas2cm(burst_m) / 100, 2), "m")

    print("\n=== ROLLING 5-YEAR RATES (IERS finals) ===")
    s = store["finals"]
    t = s["t"]
    smr = rolling(t, s["sx"], s["sy"])
    rawr = rolling(t, s["xp"], s["yp"])
    print("smoothed max", round(smr[:, 1].max(), 1), "mas/yr at",
          round(smr[smr[:, 1].argmax(), 0], 1))
    print("smoothed 2010-2020",
          [round(v, 1) for y, v in smr[(smr[:, 0] >= 2010) & (smr[:, 0] <= 2020)]])
    print("smoothed 2022+",
          [round(v, 1) for y, v in smr[(smr[:, 0] >= 2022)]])
    print("raw trailing range", round(rawr[:, 1].min(), 1), "-",
          round(rawr[:, 1].max(), 1), "mas/yr; latest",
          round(rawr[rawr[:, 0].argmax(), 1], 1))
    print("raw 2022+",
          [round(v, 1) for y, v in rawr[(rawr[:, 0] >= 2022)]])

    # ---------- FIGURES ----------
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Helvetica", "DejaVu Sans"],
                         "pdf.fonttype": 42, "axes.linewidth": 0.8})
    COL = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
           "grey": "#999999", "sky": "#56B4E9", "vermilion": "#D55E00"}

    pos = yearly_means(s["years"], s["sx"], s["sy"])
    pos_disp = pos.copy()
    pos_disp[:, 1] -= pos[0, 1]
    pos_disp[:, 2] -= pos[0, 2]
    fig, ax = plt.subplots(figsize=(7.4, 6.0))
    epochs = [("early era", 1973, 2006, COL["grey"]),
              ("2006-2020", 2006, 2020, COL["sky"]),
              ("Oct 2022-Nov 2026", 2023, 2027, COL["vermilion"])]
    for lbl, a, b, c in epochs:
        seg = pos_disp[(pos_disp[:, 0] >= a) & (pos_disp[:, 0] < b)]
        ax.plot(seg[:, 1], seg[:, 2], color=c, lw=2.0,
                marker="o", ms=3, label=lbl)
        if len(seg) >= 2:
            ax.annotate("", xy=(seg[-1, 1], seg[-1, 2]),
                        xytext=(seg[0, 1], seg[0, 2]),
                        arrowprops=dict(arrowstyle="-|>", color=c, lw=1.6))
    ax.set_xlabel("x_pole displacement from 1973 (arcsec)")
    ax.set_ylabel("y_pole displacement from 1973 (arcsec)")
    ax.set_title("Mean-pole trajectory, displacement from 1973 yearly mean\n"
                 "IERS finals 1973-2026", fontsize=11)
    ax.grid(alpha=0.3, lw=0.5)
    ax.legend(loc="upper left", fontsize=8)
    fig.savefig(FIGDIR + "/fig1_trajectory.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    labels = ["early", "middle 2010-2020", "burst"]
    xpos = np.arange(3)
    w = 0.26
    dslist = list(SNAPSHOT)
    cols = [COL["blue"], COL["orange"], COL["green"]]
    for i, ds in enumerate(dslist):
        s2 = store[ds]
        t2, te2, t02 = s2["t"], s2["t"][-1], s2["t"][0]
        rates, los, his = [], [], []
        for lbl, a, b in [("early", t02, t02 + 33),
                          ("middle", 2010, 2020),
                          ("burst", te2 - 4.12, te2)]:
            rm = slope(t2, s2["sx"], s2["sy"], a, b)[2]
            seed = int(hashlib.sha256((ds + ":" + lbl).encode()).hexdigest()[:8], 16)
            bv = boot(t2, s2["sx"], s2["sy"], a, b, seed=seed)
            rates.append(rm)
            los.append(np.percentile(bv, 2.5))
            his.append(np.percentile(bv, 97.5))
        ax.bar(xpos + (i - 1) * w, rates, width=w, label=ds,
               color=cols[i], edgecolor="white", lw=0.5)
        centers = xpos + (i - 1) * w
        ax.vlines(centers, los, his, color=COL["grey"], lw=1.0)
        ax.hlines(los, centers - 0.06, centers + 0.06, color=COL["grey"], lw=1.0)
        ax.hlines(his, centers - 0.06, centers + 0.06, color=COL["grey"], lw=1.0)
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels)
    ax.set_ylabel("secular drift rate (mas/yr)")
    ax.set_title("Secular drift rate per epoch with bootstrap CI\n"
                 "centered +/-2 yr smoothing", fontsize=11)
    ax.grid(axis="y", alpha=0.3, lw=0.5)
    ax.legend(loc="upper left", fontsize=8)
    ax.set_ylim(0, 12)
    fig.savefig(FIGDIR + "/fig2_epoch_rates.pdf")
    plt.close(fig)

    # Fig 5 residual rates with vector-bootstrap residual CI; early omitted
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    labels = ["middle 2010-2020", "burst full", "burst interior-only"]
    xpos = np.arange(3)
    w = 0.24
    dslist = list(SNAPSHOT)
    cols = [COL["blue"], COL["orange"], COL["green"]]
    for i, ds in enumerate(dslist):
        s2 = store[ds]
        t2, te2 = s2["t"], s2["t"][-1]
        rvals, los, his = [], [], []
        for lbl, a, b in [("middle", 2010, 2020),
                          ("burst full", te2 - 4.12, te2),
                          ("burst interior", te2 - 4.12, te2 - 2.0)]:
            dx, dy, rm, lon = slope(t2, s2["sx"], s2["sy"], a, b)
            rvals.append(residual_vector(dx, dy)[2])
            seed = int(hashlib.sha256((ds + ":fig5:" + lbl).encode()).hexdigest()[:8], 16)
            bx, by = boot_vec(t2, s2["sx"], s2["sy"], a, b, seed=seed)
            rb = np.array([residual_vector(x, y)[2] for x, y in zip(bx, by)])
            los.append(np.percentile(rb, 2.5)); his.append(np.percentile(rb, 97.5))
        centers = xpos + (i - 1) * w
        ax.bar(centers, rvals, width=w, label=ds,
               color=cols[i], edgecolor="white", lw=0.5)
        ax.errorbar(centers, rvals,
                    yerr=[np.array(rvals) - np.array(los), np.array(his) - np.array(rvals)],
                    fmt="none", ecolor=COL["grey"], capsize=3, lw=1.0)
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels)
    ax.set_ylabel("residual rate after known subtraction (mas/yr)")
    ax.set_title("Residual rates with vector-bootstrap CI\n"
                 "Adhikari & Ivins 2016 benchmark; early era omitted", fontsize=11)
    ax.grid(axis="y", alpha=0.3, lw=0.5)
    ax.legend(loc="upper left", fontsize=8)
    ax.set_ylim(0, 8)
    fig.savefig(FIGDIR + "/fig5_residual.pdf")
    plt.close(fig)

    # Fig 6 residual rate trajectory in rate space (trailing 5-year windows)
    rx_list, ry_list, yrs_list = [], [], []
    for y in np.arange(t[0] + 5, t[-1]):
        m = (t >= y - 5) & (t < y)
        if m.sum() < 300:
            continue
        dx, dy, rm, lon = slope(t, s["sx"], s["sy"], y - 5, y)
        rx, ry, _, _ = residual_vector(dx, dy)
        rx_list.append(rx); ry_list.append(ry); yrs_list.append(y)
    rx = np.array(rx_list); ry = np.array(ry_list); yrs = np.array(yrs_list)
    fig, ax = plt.subplots(figsize=(7.4, 6.0))
    ax.plot(rx, ry, color=COL["grey"], lw=0.8, alpha=0.6)
    epochs = [("early era", 1973, 2006, COL["grey"]),
              ("2006-2020", 2006, 2020, COL["sky"]),
              ("Oct 2022-Nov 2026", 2023, 2027, COL["vermilion"])]
    for lbl, a, b, c in epochs:
        mask = (yrs >= a) & (yrs < b)
        ax.scatter(rx[mask], ry[mask], color=c, s=12, label=lbl)
    ax.set_xlabel("residual x rate (mas/yr)")
    ax.set_ylabel("residual y rate (mas/yr)")
    ax.set_title("Residual rate trajectory after subtracting climate+GIA vectors\n"
                 "IERS finals, trailing 5-year windows", fontsize=11)
    ax.grid(alpha=0.3, lw=0.5)
    ax.legend(loc="upper left", fontsize=8)
    fig.savefig(FIGDIR + "/fig6_residual_trajectory.pdf")
    plt.close(fig)

    # Fig 7 residual displacement path (increments corrected, 2003-2026)
    pos = yearly_means(s["years"], s["sx"], s["sy"])
    seg = pos[pos[:, 0] >= 2003]
    inc_x = np.diff(seg[:, 1]); inc_y = np.diff(seg[:, 2])
    kx, ky = known_components()
    kx, ky = kx / 1000.0, ky / 1000.0   # mas/yr -> arcsec/yr
    cum_x = np.cumsum(inc_x - kx) * 1000.0
    cum_y = np.cumsum(inc_y - ky) * 1000.0
    yrs = seg[1:, 0]
    fig, ax = plt.subplots(figsize=(7.4, 6.0))
    epochs = [("2003-2006", 2003, 2006, COL["grey"]),
              ("2006-2020", 2006, 2020, COL["sky"]),
              ("Oct 2022-Nov 2026", 2023, 2027, COL["vermilion"])]
    ax.plot(cum_x, cum_y, color=COL["grey"], lw=1.0, alpha=0.7)
    for lbl, a, b, c in epochs:
        mask = (yrs >= a) & (yrs < b)
        ax.scatter(cum_x[mask], cum_y[mask], color=c, s=18, label=lbl)
    ax.set_xlabel("residual x displacement (mas)")
    ax.set_ylabel("residual y displacement (mas)")
    ax.set_title("Residual displacement path after subtracting climate+GIA\n"
                 "IERS finals, yearly increments, 2003-2026 benchmark", fontsize=11)
    ax.set_xlim(cum_x.min() - 2, cum_x.max() + 2)
    ax.set_ylim(cum_y.min() - 2, cum_y.max() + 3)
    ax.grid(alpha=0.3, lw=0.5)
    ax.legend(loc="lower right", fontsize=8)
    fig.savefig(FIGDIR + "/fig7_residual_path.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    ax.plot(rawr[:, 0], rawr[:, 1], color=COL["grey"], lw=1.0,
            label="raw series")
    ax.plot(smr[:, 0], smr[:, 1], color=COL["blue"], lw=2.0,
            label="smoothed (+/-2 yr)")
    ax.axvspan(2022.77, 2026.89, color=COL["vermilion"], alpha=0.10,
               label="burst window")
    ax.set_xlabel("year")
    ax.set_ylabel("trailing 5-year secular rate (mas/yr)")
    ax.set_title("Trailing 5-year secular rate, raw vs smoothed\n"
                 "IERS finals", fontsize=11)
    ax.grid(alpha=0.3, lw=0.5)
    ax.legend(loc="upper left", fontsize=8)
    fig.savefig(FIGDIR + "/fig3_rolling_rate.pdf")
    plt.close(fig)

    net = np.array([math.hypot(p[1] - pos[0, 1],
                                p[2] - pos[0, 2]) * 1000 for p in pos])
    net_m = np.array([mas2cm(v) / 100 for v in net])
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    ax.plot(pos[:, 0], net_m, color=COL["blue"], lw=2.0)
    ax.fill_between(pos[:, 0], net_m, color=COL["sky"], alpha=0.25)
    ax.axvspan(2022.77, 2026.89, color=COL["vermilion"], alpha=0.10,
               label="burst window")
    ax.set_xlabel("year")
    ax.set_ylabel("net displacement since 1973 (m)")
    ax.set_title("Cumulative net mean-pole displacement\n"
                 "IERS finals, smoothed yearly means", fontsize=11)
    ax.grid(alpha=0.3, lw=0.5)
    ax.legend(loc="upper left", fontsize=8)
    fig.savefig(FIGDIR + "/fig4_cumulative.pdf")
    plt.close(fig)

    print("\nFigures written to figures/: fig1_trajectory.pdf, "
          "fig2_epoch_rates.pdf, fig3_rolling_rate.pdf, fig4_cumulative.pdf")

if __name__ == "__main__":
    main()
