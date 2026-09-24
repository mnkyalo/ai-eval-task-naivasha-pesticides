"""Honest variants and plausible wrong answers for the Naivasha passive-sampling task.

Every function takes the data directory and returns (answers, F) where answers maps
(site, compound) -> ng/L and F maps site -> fitted sheet coefficient. None of them
reads a generating parameter or the answer key.

The reference solution itself is not re-implemented here; run_validation.py executes
task/solution/reference_solution.py unmodified in a scratch directory.
"""
import csv
import math
import os
from collections import defaultdict

import numpy as np
from scipy.optimize import least_squares

T = 30.0
WIN = {"A": 27.0, "B": 20.0, "C": 10.0, "D": 0.0}
AREA = 2 * 5.5 * 9.0                          # cm2, both faces
MASS = (5.5 * 9.0 * 0.05) * 1.15 / 1000.0     # kg
EXPN = 0.47


def _read(data, name):
    with open(os.path.join(data, name)) as f:
        return list(csv.DictReader(f))


def _mean_by(rows, key):
    acc = defaultdict(list)
    for r in rows:
        acc[r[key]].append(float(r["mass_ng"]))
    return {k: sum(v) / len(v) for k, v in acc.items()}


def load(data, blank_correct=True):
    props = {r["name"]: (float(r["log_Kpw"]), float(r["molar_mass_g_per_mol"]), r["role"])
             for r in _read(data, "compound_properties.csv")}
    blank = _mean_by(_read(data, "procedural_blanks.csv"), "analyte")
    sdb = _mean_by(_read(data, "speedisk_procedural_blanks.csv"), "analyte")
    if not blank_correct:
        blank = {k: 0.0 for k in blank}
        sdb = {k: 0.0 for k in sdb}
    n0 = {k: v - blank[k] for k, v in _mean_by(_read(data, "prc_reference_samplers.csv"), "prc").items()}
    prc = defaultdict(list)
    for r in _read(data, "prc_exposed_samplers.csv"):
        prc[(r["site"], r["prc"], float(r["deployment_days"]))].append(float(r["mass_ng"]) - blank[r["prc"]])
    tgt = defaultdict(list)
    for r in _read(data, "target_compound_masses.csv"):
        tgt[(r["site"], r["compound"], r["sampler_id"].split("-")[1][0])].append(
            float(r["mass_ng"]) - blank[r["compound"]])
    sd = defaultdict(list)
    for r in _read(data, "speedisk_compound_masses.csv"):
        sd[(r["site"], r["compound"])].append(float(r["mass_ng"]) - sdb[r["compound"]])
    sites = sorted({k[0] for k in tgt})
    compounds = [n for n, v in props.items() if v[2] == "target"]
    return props, n0, prc, tgt, sd, sites, compounds


def rs(F, M):
    return F * AREA / M ** EXPN


def fit_F(props, n0, prc, sites):
    F = {}
    for s in sites:
        pts = [(props[p][0], props[p][1], d, v / n0[p])
               for (site, p, d), vals in prc.items() if site == s for v in vals]
        res = least_squares(lambda q: [math.exp(-rs(q[0], M) * d / (10 ** lk * MASS)) - f
                                       for lk, M, d, f in pts], [1.0], bounds=(1e-4, 50))
        F[s] = float(res.x[0])
    return F


def _sheet_steady_state(tgt, s, c, K, ke, sset="D"):
    """Classic constant-concentration inversion of one sheet set."""
    dur = T - WIN[sset]
    return np.mean(tgt[(s, c, sset)]) / (K * (1 - math.exp(-ke * dur)))


def _speedisk_route(data, calib, with_phi, blank_correct=True):
    """Sheets for the compounds they still integrate, Speedisk for the rest.

    calib chooses which integrative compounds calibrate the disk:
      'all'      every compound whose 30-day sheet is below 30% of equilibrium
      'least'    only the least hydrophobic of those
    with_phi fits the DOC binding term G*(1+phi*Kpw); without it, one G per site.
    """
    props, n0, prc, tgt, sd, sites, compounds = load(data, blank_correct)
    F = fit_F(props, n0, prc, sites)
    out = {}
    for s in sites:
        info = {}
        for c in compounds:
            lk, M, _ = props[c]
            K = 10 ** lk * MASS
            ke = rs(F[s], M) / K
            info[c] = (K, ke, M, 10 ** lk)
        integ = sorted([c for c in compounds if 1 - math.exp(-info[c][1] * T) < 0.3],
                       key=lambda c: info[c][3])
        if calib == "least":
            integ = integ[:1]
        # sheet-derived 30-day TWA for the integrative compounds (set D, full equation)
        tw = {c: _sheet_steady_state(tgt, s, c, info[c][0], info[c][1]) for c in integ}
        ratio = {c: np.mean(sd[(s, c)]) * info[c][2] ** EXPN / (tw[c] * T) for c in integ}
        if with_phi and len(integ) >= 2:
            kp = np.array([info[c][3] for c in integ])
            r = np.array([ratio[c] for c in integ])
            fit = least_squares(lambda q: np.log(q[0] * (1 + q[1] * kp)) - np.log(r),
                                [np.median(r), 1e-7], bounds=([1e-3, 0], [1e5, 1e-3]))
            G, phi = fit.x
        else:
            G, phi = float(np.exp(np.mean(np.log(list(ratio.values()))))), 0.0
        for c in compounds:
            if c in tw and calib != "least":
                out[(s, c)] = tw[c]
            else:
                K, ke, M, kpw = info[c]
                out[(s, c)] = np.mean(sd[(s, c)]) * M ** EXPN / (G * (1 + phi * kpw) * T)
    return out, F


# ------------------------------------------------------------------ honest
def honest_disk_with_binding_term(data):
    """Sheets integrate the hydrophobic compounds; disk, corrected for the DOC term, the rest."""
    return _speedisk_route(data, "all", with_phi=True)


# ------------------------------------------------------------------ wrong
def wrong_textbook_one_G(data):
    """One Speedisk coefficient per site over the integrative compounds, no binding term."""
    return _speedisk_route(data, "all", with_phi=False)


def wrong_least_hydrophobic_calibrant(data):
    """G from the single least hydrophobic integrative compound."""
    return _speedisk_route(data, "least", with_phi=False)


def wrong_steady_state_30day_set(data):
    """Ignore the concentration decline: invert the 30-day set at constant Cw."""
    props, n0, prc, tgt, sd, sites, compounds = load(data)
    F = fit_F(props, n0, prc, sites)
    out = {}
    for s in sites:
        for c in compounds:
            lk, M, _ = props[c]
            K = 10 ** lk * MASS
            out[(s, c)] = _sheet_steady_state(tgt, s, c, K, rs(F[s], M) / K)
    return out, F


def wrong_no_blank_correction(data):
    """Honest route but with procedural blanks ignored."""
    return _speedisk_route(data, "all", with_phi=True, blank_correct=False)
