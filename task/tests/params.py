"""
Ground-truth parameters, v3 -- sampler-memory / non-steady-concentration mechanism,
plus co-deployed Speedisk adsorption samplers as an independent integrative measurement.
NOT solver-visible. Provenance: [LIT] literature, [COR] correlation-derived, [DES] design.
"""
SHEET_L_MM=55.0; SHEET_W_MM=90.0; SHEET_T_MM=0.5      # [LIT] Abbasi & Mannaerts 2018
SILICONE_RHO=1.15                                      # [LIT] g/cm3, AlteSil silicone
N_REPLICATES=3                                         # [LIT]
RS_EXPONENT=0.47                                       # [LIT] Rs = F*A/M^0.47
KPW_SLOPE=1.06; KPW_INTERCEPT=-1.16                    # [COR] logKpw = 1.06*logKow - 1.16

# Four deployments, staggered starts, ALL retrieved together on 20 July 2016. [DES]
# End-anchored by design: start-anchored windows let a solver average two steady-state
# inversions and have the biases cancel. Windows are (start_day, end_day), day 0 = 20 Jun.
# Four windows. The extra early-period resolution matters: with only three, an
# exponential and a piecewise-constant fit diverge by ~6% on the fast-equilibrating
# compounds (measured), which penalises a solver for a defensible model choice.
DEPLOYMENTS={"A":(27.0,30.0), "B":(20.0,30.0), "C":(10.0,30.0), "D":(0.0,30.0)}
DEPLOY_DATES={"A":"17 July 2016","B":"10 July 2016","C":"30 June 2016","D":"20 June 2016"}
RETRIEVAL_DATE="20 July 2016"

# F (hydrodynamic, per site) and lambda (concentration decay, per site, /day). [DES]
# lambda chosen so each site shows a real but different decline over the 30 days.
SITES={
 "Upper Malewa":  {"F":0.92, "lam":0.026280, "code":"UM"},   # 2.2x decline
 "Middle Malewa": {"F":1.93, "lam":0.050136, "code":"MM"},   # 4.5x decline
 "Lake Naivasha": {"F":0.28, "lam":0.021402, "code":"LN"},   # 1.9x decline
}

PRCS={  # name: (logKow [LIT], molar mass [LIT])
 "BIP-D10":(4.01,164.28),"PCB001":(4.46,188.65),"PCB002":(4.69,188.65),
 "PCB003":(4.63,188.65),"PCB010":(4.84,223.10),"PCB014":(5.28,223.10),
 "PCB030":(5.44,257.54),"PCB021":(5.51,257.54),"PCB050":(5.95,291.99),
 "PCB055":(6.11,291.99),"PCB104":(6.34,326.43),"PCB078":(6.35,291.99),
 "PCB145":(6.72,360.88),"PCB204":(7.30,429.77),
}
PRC_SPIKE_NG=120.0

TARGETS={
 "endosulfan sulfate":(3.66,422.92),"gamma-HCH":(3.72,290.83),"beta-HCH":(3.78,290.83),
 "alpha-HCH":(3.80,290.83),"alpha-endosulfan":(4.74,406.93),"heptachlor epoxide":(5.00,389.32),
 "methoxychlor":(5.08,345.65),"endrin":(5.20,380.91),"dieldrin":(5.40,380.91),
 "pp-DDD":(6.02,320.04),"heptachlor":(6.10,373.32),"aldrin":(6.50,364.91),
 "pp-DDE":(6.51,318.03),"pp-DDT":(6.91,354.49),
 "HCB":(5.73,284.78),"cis-chlordane":(6.16,409.78),"trans-chlordane":(6.22,409.78),
 "trans-nonachlor":(6.35,444.22),"mirex":(7.18,545.54),
}

# C_0 = freely dissolved concentration at deployment start, ng/L. [DES]
# Per-compound site patterns, deliberately not a single global site multiplier.
TRUE_C0={
 #                      Upper Malewa  Middle Malewa  Lake Naivasha
 "endosulfan sulfate": ( 11.20,        34.80,         41.50),
 "gamma-HCH":          (  8.40,        26.10,          9.30),
 "beta-HCH":           (  4.75,         7.20,          5.60),
 "alpha-HCH":          ( 11.60,        21.40,         12.10),
 "alpha-endosulfan":   (  4.90,        18.70,          2.55),
 "heptachlor epoxide": (  0.64,         1.42,          0.91),
 "methoxychlor":       (  1.52,         6.05,          0.84),
 "endrin":             (  0.44,         0.86,          0.49),
 "dieldrin":           (  0.88,         2.41,          1.18),
 "pp-DDD":             (  0.52,         1.00,          1.07),
 "heptachlor":         (  0.31,         0.73,          0.31),
 "aldrin":             (  0.19,         0.61,          0.13),
 "pp-DDE":             (  0.74,         1.13,          1.85),
 "pp-DDT":             (  0.38,         0.95,          0.23),
 # v5: five further legacy organochlorines routinely reported in East African surveys,
 # extending the integrative end of the suite so every site carries at least five
 # compounds the sheets still integrate at 30 days
 "HCB":                (  0.42,         0.88,          0.61),
 "cis-chlordane":      (  0.27,         0.66,          0.34),
 "trans-chlordane":    (  0.33,         0.79,          0.41),
 "trans-nonachlor":    (  0.21,         0.52,          0.29),
 "mirex":              (  0.09,         0.24,          0.12),
}

SEED=20160620
RSD_PRC_REF=0.05; RSD_PRC_EXPOSED=0.07; RSD_TARGET=0.07
CARRYOVER_LO=0.015; CARRYOVER_HI=0.045
TARGET_BLANK_LO=0.05; TARGET_BLANK_HI=0.30
RSD_BLANK=0.18; N_BLANKS=3

# ---- v3: Speedisk adsorption samplers, co-deployed for the full 30 days ------------
# [LIT] Abbasi & Mannaerts 2018 mounted three Speedisks (H2O-philic DVB, 0.6 g) beside
# the three silicone sheets at every site, 20 June - 20 July 2016, without PRCs.
# Adsorptive uptake far below sorbent capacity is linear for every compound, so the
# disk mass is Rs_SD * integral(Cw dt) whatever the concentration history was.
# [DES] Rs_SD = G / M^0.47 (water-boundary-layer control, same mass dependence as the
# sheets); G is per site, tracks the hydrodynamics but is NOT a fixed multiple of F.
# Drawn from a SEPARATE random stream so every v2 file stays byte-identical.
SPEEDISK_G={"Upper Malewa":28.0, "Middle Malewa":43.0, "Lake Naivasha":11.5}   # G/(F*A): 0.31, 0.23, 0.42 -- no shared ratio

# ---- v5: DOC-associated uptake by the adsorptive disk. [LIT] Adsorption-based samplers
# (SPE disks, POCIS) take up compound carried on dissolved organic matter, which a
# partitioning sampler (silicone) does not; the excess scales with the compound's
# DOC-binding constant, itself proportional to hydrophobicity (log Kdoc ~ log Kow - 0.5).
# Disk mass = (G/M^0.47) * TWA_free * T * (1 + PHI_site * Kpw), PHI_site absorbing the
# DOC level and capture efficiency at that site. [DES] Expressed through the enhancement
# at pp-DDT (log Kpw 6.16): 1 + PHI*Kpw_DDT. No DOC data are shipped -- the site-level
# binding term is a second calibration parameter the solver has to recognise and fit.
DOC_ENH_AT_DDT={"Upper Malewa":2.3, "Middle Malewa":2.9, "Lake Naivasha":2.6}
LOGKPW_DDT=6.16
SPEEDISK_WINDOW=(0.0,30.0); SPEEDISK_DEPLOYED="20 June 2016"
SPEEDISK_SEED_OFFSET=1
RSD_SPEEDISK=0.05
SD_BLANK_LO=0.05; SD_BLANK_HI=0.40

# grading: 20% relative on the concentrations with a 95% pass bar (54 of 57), F within 10%
CW_REL_TOL=0.20; CW_MIN_PASS=54; F_REL_TOL=0.10
