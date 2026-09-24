#!/usr/bin/env python3
"""
Reference solution: PRC-calibrated silicone sheets plus co-deployed Speedisks, under a
water concentration that was not steady.

  1. Blank-correct every measured mass with the procedural blanks of its own sampler
     type (sheet blanks for sheets and reference sheets, Speedisk blanks for Speedisks).
  2. Retained fraction f = (N_exposed - blank) / (N_reference - blank) per PRC.
  3. Fit one F per site by NLS on f, over all PRCs, all sheets and all deployment lengths.
  4. Per compound, ke = Rs/(Kpw*m) with Rs = F*A/M^0.47.
  5. The four sheet sets share a retrieval date but differ in start, so each sheet is an
     exponentially-weighted integral over a different window. A Speedisk is a plain
     integral over the full 30 days -- but an adsorptive disk also collects compound
     carried on dissolved organic matter, which the partitioning sheet does not see, and
     that excess scales with the compound's hydrophobicity:
         N_SD = (G/M^0.47) * TWA_free * 30 * (1 + phi*Kpw)
     with phi a site-level binding term. Fit, per site, one decline rate lambda, one
     Speedisk coefficient G, one binding term phi and 19 per-compound amplitudes to
     every sheet and every Speedisk together. The compounds the sheets still integrate
     (mirex down to HCB, log Kpw 6.4 to 4.9) span the enhancement from about 3x down to
     nearly 1, which is what pins G and phi separately; the equilibrated compounds sit
     at enhancement 1.00 and take G alone.
     For the compounds that equilibrate on the sheet in a day or two, the sheets fix
     the end-of-campaign concentration and the Speedisk fixes the 30-day integral, so
     the shared-decline model is tested by the data instead of being assumed: if those
     compounds had not followed the site profile, their Speedisk residuals would show it.
     That check is printed.
  6. Report the 30-day time-weighted average.

Residuals are on the log scale and weighted by the replicate scatter of each sampler
type, measured from the data. Nothing here reads the generating parameters.
"""
import csv, math, os
from collections import defaultdict
import numpy as np
from scipy.optimize import least_squares

DATA = "/app/data"
T = 30.0
WIN = {"A": 27.0, "B": 20.0, "C": 10.0, "D": 0.0}      # window start, day 0 = 20 Jun
AREA = 2 * 5.5 * 9.0                                    # cm2, both faces
MASS = (5.5 * 9.0 * 0.05) * 1.15 / 1000.0               # kg
EXPN = 0.47

def read(name):
    with open(os.path.join(DATA, name)) as f:
        return list(csv.DictReader(f))

def mean_by(rows, key, val="mass_ng"):
    acc = defaultdict(list)
    for r in rows:
        acc[r[key]].append(float(r[val]))
    return {k: sum(v) / len(v) for k, v in acc.items()}

props = {r["name"]: (float(r["log_Kpw"]), float(r["molar_mass_g_per_mol"]), r["role"])
         for r in read("compound_properties.csv")}
BLANK = mean_by(read("procedural_blanks.csv"), "analyte")
SD_BLANK = mean_by(read("speedisk_procedural_blanks.csv"), "analyte")
N0 = {k: v - BLANK[k] for k, v in mean_by(read("prc_reference_samplers.csv"), "prc").items()}

prc = defaultdict(list)
for r in read("prc_exposed_samplers.csv"):
    prc[(r["site"], r["prc"], float(r["deployment_days"]))].append(float(r["mass_ng"]) - BLANK[r["prc"]])
tgt = defaultdict(list)
for r in read("target_compound_masses.csv"):
    tgt[(r["site"], r["compound"], r["sampler_id"].split("-")[1][0])].append(float(r["mass_ng"]) - BLANK[r["compound"]])
sd = defaultdict(list)
for r in read("speedisk_compound_masses.csv"):
    sd[(r["site"], r["compound"])].append(float(r["mass_ng"]) - SD_BLANK[r["compound"]])

sites = sorted({k[0] for k in tgt})
compounds = [n for n, v in props.items() if v[2] == "target"]

def rs(F, M):
    return F * AREA / (M ** EXPN)

def pooled_log_sd(groups):
    """Replicate scatter on the log scale, pooled over groups of replicates."""
    ss, dof = 0.0, 0
    for v in groups:
        v = [math.log(x) for x in v if x > 0]
        if len(v) > 1:
            m = sum(v) / len(v); ss += sum((x - m) ** 2 for x in v); dof += len(v) - 1
    return math.sqrt(ss / dof)

S_SHEET = pooled_log_sd(tgt.values())
S_DISK = pooled_log_sd(sd.values())

# ---- step 3: one F per site, every PRC on every sheet
F = {}
for s in sites:
    pts = [(props[p][0], props[p][1], d, v / N0[p])
           for (site, p, d), vals in prc.items() if site == s for v in vals]
    def resid(q):
        return [math.exp(-rs(q[0], M) * d / ((10 ** lk) * MASS)) - f for lk, M, d, f in pts]
    F[s] = float(least_squares(resid, [1.0], bounds=(1e-4, 50)).x[0])

# ---- steps 4-6: shared decline + per-compound amplitude + Speedisk coefficient, per site
def sheet_pred(C0, lam, t0, ke, K):
    if abs(ke - lam) < 1e-12:
        return K * ke * C0 * (T - t0) * math.exp(-lam * T)
    return K * ke * C0 * (math.exp(-lam * T) - math.exp(-lam * t0) * math.exp(-ke * (T - t0))) / (ke - lam)

def twa_of(C0, lam):
    return C0 * (1 - math.exp(-lam * T)) / (lam * T)

rows, notes = [], []
for s in sites:
    ke, K, Mm = [], [], []
    for c in compounds:
        lk, M, _ = props[c]
        Kc = (10 ** lk) * MASS
        K.append(Kc); ke.append(rs(F[s], M) / Kc); Mm.append(M)

    Kpw = [10 ** props[c][0] for c in compounds]
    NP = 3                                                # lambda, G, phi

    def resid(p, with_disks=True):
        lam, G, phi = p[0], p[1], p[2]
        out = []
        for i, c in enumerate(compounds):
            for st, t0 in WIN.items():
                pr = sheet_pred(p[NP + i], lam, t0, ke[i], K[i])
                out += [math.log(v / pr) / S_SHEET for v in tgt[(s, c, st)] if v > 0]
            if with_disks:
                pr = G * (1 + phi * Kpw[i]) / Mm[i] ** EXPN * twa_of(p[NP + i], lam) * T
                out += [math.log(v / pr) / S_DISK for v in sd[(s, c)] if v > 0]
        return out

    lo = [1e-5, 1e-3, 0.0] + [1e-9] * len(compounds)
    hi = [0.5, 1e5, 1.0] + [1e6] * len(compounds)
    best = None
    for lam0, phi0 in ((0.005, 0.0), (0.03, 1e-7), (0.12, 1e-6)):   # the optimum must not depend on the start
        p0 = [lam0, 20.0, phi0] + [max(np.mean(tgt[(s, c, "D")]) / K[i], 1e-3) for i, c in enumerate(compounds)]
        fit = least_squares(resid, p0, bounds=(lo, hi), x_scale="jac")
        if best is None or fit.cost < best.cost:
            best = fit
    lam, G, phi = best.x[0], best.x[1], best.x[2]
    for i, c in enumerate(compounds):
        rows.append((s, c, twa_of(best.x[NP + i], lam)))

    # The check that replaces an assumption, made on data the prediction never saw: refit the
    # sheets alone, calibrate G on the compounds the sheets still integrate (30-day degree of
    # equilibrium below 0.3), then ask whether the equilibrated compounds' Speedisks land where
    # the sheet-only site profile says they should.
    so = least_squares(lambda q: resid(q, with_disks=False), best.x, bounds=(lo, hi), x_scale="jac")
    unit = lambda i: np.mean(sd[(s, compounds[i])]) * Mm[i] ** EXPN / (twa_of(so.x[NP + i], so.x[0]) * T)
    slow = sorted([i for i in range(len(compounds)) if 1 - math.exp(-ke[i] * T) < 0.3], key=lambda i: Kpw[i])
    fast = [i for i in range(len(compounds)) if 1 - math.exp(-ke[i] * 3.0) > 0.95]
    app = " ".join(f"{compounds[i]}={unit(i)/G:.2f}" for i in slow)
    dev = [math.log(unit(i) / G) for i in fast]
    notes.append(f"{s}: F={F[s]:.4f} lambda={lam:.4f}/d G={G:.2f} phi*Kpw(DDT)={phi*10**6.16:.2f} | apparent/true G "
                 f"across the integrative compounds (sheet-only profile): {app} | equilibrated compounds' disks vs G: "
                 f"mean log-deviation {np.mean(dev):+.3f}")

with open("/app/answer.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["site", "compound", "cw_ng_per_L"])
    for s, c, v in rows:
        w.writerow([s, c, f"{v:.6g}"])
with open("/app/sampling_rates.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["site", "F"])
    for s in sites:
        w.writerow([s, f"{F[s]:.6g}"])
print(f"wrote {len(rows)} concentrations")
print("\n".join(notes))
