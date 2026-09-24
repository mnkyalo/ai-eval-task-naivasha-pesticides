"""
Deterministic generator v3 -- non-steady concentrations, nested silicone deployments,
co-deployed Speedisks.
Run from generator/:  python3 build.py
"""
import csv, json, math, os
import numpy as np
import params as P

HERE=os.path.dirname(os.path.abspath(__file__))
OUT_ENV=os.path.join(HERE,"..","environment","data"); OUT_TRUTH=os.path.join(HERE,"_regen_truth")

AREA_CM2=2.0*(P.SHEET_L_MM/10.0)*(P.SHEET_W_MM/10.0)
VOL_CM3=(P.SHEET_L_MM/10.0)*(P.SHEET_W_MM/10.0)*(P.SHEET_T_MM/10.0)
MASS_KG=VOL_CM3*P.SILICONE_RHO/1000.0

def logkpw(lkow): return P.KPW_SLOPE*lkow+P.KPW_INTERCEPT
def rs_of(F,M):   return F*AREA_CM2/(M**P.RS_EXPONENT)

def target_mass(C0, lam, t0, t1, lk, M, F):
    """Mass on a sampler deployed [t0,t1] (starting empty) with Cw(t)=C0*exp(-lam*t)."""
    Rs=rs_of(F,M); tau=(10**lk)*MASS_KG/Rs; a=1.0/tau
    if abs(a-lam)<1e-12: return Rs*C0*(t1-t0)*math.exp(-lam*t1)
    return Rs*C0*(math.exp(-lam*t1)-math.exp(-lam*t0)*math.exp(-(t1-t0)/tau))/(a-lam)

def prc_remaining(N0, T, lk, M, F):
    """PRCs dissipate with the same time constant; unaffected by water concentration."""
    return N0*math.exp(-rs_of(F,M)*T/((10**lk)*MASS_KG))

rng=np.random.default_rng(P.SEED)
prc_blank_ng={p: rng.uniform(P.CARRYOVER_LO,P.CARRYOVER_HI)*P.PRC_SPIKE_NG for p in P.PRCS}
tgt_blank_ng={c: rng.uniform(P.TARGET_BLANK_LO,P.TARGET_BLANK_HI) for c in P.TARGETS}
def measured(true_ng, blank_ng, rsd):
    return true_ng*rng.lognormal(0.0,rsd)+blank_ng*rng.lognormal(0.0,0.20)

blank_rows=[]
for rep in range(1,P.N_BLANKS+1):
    for p in P.PRCS:
        blank_rows.append({"sampler_id":f"BLK-{rep}","analyte":p,"role":"PRC",
                           "mass_ng":round(prc_blank_ng[p]*rng.lognormal(0,P.RSD_BLANK),3)})
    for c in P.TARGETS:
        blank_rows.append({"sampler_id":f"BLK-{rep}","analyte":c,"role":"target",
                           "mass_ng":round(tgt_blank_ng[c]*rng.lognormal(0,P.RSD_BLANK),4)})

ref_rows=[]
for p in P.PRCS:
    for rep in range(1,P.N_REPLICATES+1):
        ref_rows.append({"sampler_id":f"REF-{rep}","prc":p,
                         "mass_ng":round(measured(P.PRC_SPIKE_NG,prc_blank_ng[p],P.RSD_PRC_REF),3)})

exp_rows=[]; tgt_rows=[]
for site,cfg in P.SITES.items():
    for dname,(t0,t1) in P.DEPLOYMENTS.items():
        dur=t1-t0
        for rep in range(1,P.N_REPLICATES+1):
            sid=f"{cfg['code']}-{dname}{rep}"
            for p,(lkow,M) in P.PRCS.items():
                v=measured(prc_remaining(P.PRC_SPIKE_NG,dur,logkpw(lkow),M,cfg["F"]),
                           prc_blank_ng[p],P.RSD_PRC_EXPOSED)
                exp_rows.append({"site":site,"sampler_id":sid,"deployed":P.DEPLOY_DATES[dname],
                                 "deployment_days":int(dur),"prc":p,"mass_ng":round(v,3)})
            for c,(lkow,M) in P.TARGETS.items():
                si=list(P.SITES).index(site)
                v=measured(target_mass(P.TRUE_C0[c][si],cfg["lam"],t0,t1,logkpw(lkow),M,cfg["F"]),
                           tgt_blank_ng[c],P.RSD_TARGET)
                tgt_rows.append({"site":site,"sampler_id":sid,"deployed":P.DEPLOY_DATES[dname],
                                 "deployment_days":int(dur),"compound":c,"mass_ng":round(v,4)})

prop_rows=[{"name":p,"role":"PRC","log_Kpw":round(logkpw(v[0]),2),"molar_mass_g_per_mol":v[1]}
           for p,v in P.PRCS.items()]
prop_rows+=[{"name":c,"role":"target","log_Kpw":round(logkpw(v[0]),2),"molar_mass_g_per_mol":v[1]}
            for c,v in P.TARGETS.items()]

def wcsv(path,rows,fields):
    with open(path,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
os.makedirs(OUT_ENV,exist_ok=True); os.makedirs(OUT_TRUTH,exist_ok=True)
wcsv(f"{OUT_ENV}/procedural_blanks.csv",blank_rows,["sampler_id","analyte","role","mass_ng"])
wcsv(f"{OUT_ENV}/prc_reference_samplers.csv",ref_rows,["sampler_id","prc","mass_ng"])
wcsv(f"{OUT_ENV}/prc_exposed_samplers.csv",exp_rows,["site","sampler_id","deployed","deployment_days","prc","mass_ng"])
wcsv(f"{OUT_ENV}/target_compound_masses.csv",tgt_rows,["site","sampler_id","deployed","deployment_days","compound","mass_ng"])
wcsv(f"{OUT_ENV}/compound_properties.csv",prop_rows,["name","role","log_Kpw","molar_mass_g_per_mol"])

# ---- v3: Speedisks (separate random stream; nothing above this line changes) -------
rng_sd=np.random.default_rng(P.SEED+P.SPEEDISK_SEED_OFFSET)
sd_blank_ng={c: rng_sd.uniform(P.SD_BLANK_LO,P.SD_BLANK_HI) for c in P.TARGETS}
def integral_cw(C0,lam,t0,t1): return C0*(math.exp(-lam*t0)-math.exp(-lam*t1))/lam
sd_blank_rows=[]
for rep in range(1,P.N_BLANKS+1):
    for c in P.TARGETS:
        sd_blank_rows.append({"sampler_id":f"SDB-{rep}","analyte":c,
                              "mass_ng":round(sd_blank_ng[c]*rng_sd.lognormal(0,P.RSD_BLANK),4)})
sd_rows=[]
for site,cfg in P.SITES.items():
    si=list(P.SITES).index(site); t0,t1=P.SPEEDISK_WINDOW
    for rep in range(1,P.N_REPLICATES+1):
        for c,(lkow,M) in P.TARGETS.items():
            phi=(P.DOC_ENH_AT_DDT[site]-1.0)/(10**P.LOGKPW_DDT)
            enh=1.0+phi*(10**logkpw(lkow))
            true=P.SPEEDISK_G[site]/(M**P.RS_EXPONENT)*integral_cw(P.TRUE_C0[c][si],cfg["lam"],t0,t1)*enh
            v=true*rng_sd.lognormal(0.0,P.RSD_SPEEDISK)+sd_blank_ng[c]*rng_sd.lognormal(0.0,0.20)
            sd_rows.append({"site":site,"sampler_id":f"{cfg['code']}-SD{rep}","deployed":P.SPEEDISK_DEPLOYED,
                            "deployment_days":int(t1-t0),"compound":c,"mass_ng":round(v,4)})
wcsv(f"{OUT_ENV}/speedisk_compound_masses.csv",sd_rows,["site","sampler_id","deployed","deployment_days","compound","mass_ng"])
wcsv(f"{OUT_ENV}/speedisk_procedural_blanks.csv",sd_blank_rows,["sampler_id","analyte","mass_ng"])

# time-weighted mean Cw over the 30-day period = C0*(1-exp(-lam*T))/(lam*T)
def twa(C0,lam,T): return C0*(1-math.exp(-lam*T))/(lam*T)
truth={"seed":P.SEED,
 "sampler":{"area_cm2":AREA_CM2,"volume_cm3":VOL_CM3,"mass_kg":MASS_KG},
 "deployments":P.DEPLOYMENTS,
 "true_F":{s:c["F"] for s,c in P.SITES.items()},
 "true_lambda":{s:c["lam"] for s,c in P.SITES.items()},
 "decline_factor_30d":{s:round(math.exp(c["lam"]*30),3) for s,c in P.SITES.items()},
 "prc_blank_ng":{k:round(v,4) for k,v in prc_blank_ng.items()},
 "target_blank_ng":{k:round(v,4) for k,v in tgt_blank_ng.items()},
 "true_C0_ng_per_L":{c:{s:P.TRUE_C0[c][i] for i,s in enumerate(P.SITES)} for c in P.TARGETS},
 "true_TWA_30d_ng_per_L":{c:{s:round(twa(P.TRUE_C0[c][i],P.SITES[s]["lam"],30.0),5)
                             for i,s in enumerate(P.SITES)} for c in P.TARGETS},
 "speedisk_G":P.SPEEDISK_G,
 "doc_enhancement_at_DDT":P.DOC_ENH_AT_DDT,
 "speedisk_blank_ng":{k:round(v,4) for k,v in sd_blank_ng.items()},
 "log_Kpw":{n:round(logkpw(v[0]),4) for n,v in {**P.PRCS,**P.TARGETS}.items()}}
json.dump(truth,open(f"{OUT_TRUTH}/ground_truth.json","w"),indent=2)

print(f"mass {MASS_KG*1000:.4f} g | area {AREA_CM2:.0f} cm2 | windows {P.DEPLOYMENTS}")
for s,c in P.SITES.items():
    print(f"  {s:<15} F={c['F']:.2f}  lambda={c['lam']:.5f}/d  "
          f"decline over 30d = {math.exp(c['lam']*30):.2f}x")
print(f"rows: blanks={len(blank_rows)} ref={len(ref_rows)} prc_exp={len(exp_rows)} targets={len(tgt_rows)}")

# ---- frozen answer key for the verifier: generator truth, tolerance from params
key={"cw_ng_per_L":{s:{c:round(twa(P.TRUE_C0[c][i],P.SITES[s]["lam"],30.0),5) for c in P.TARGETS} for i,s in enumerate(P.SITES)},
     "F":{s:c["F"] for s,c in P.SITES.items()},
     "tolerance":{"cw_rel":P.CW_REL_TOL,"cw_min_pass":P.CW_MIN_PASS,"cw_total":3*len(P.TARGETS),"F_rel":P.F_REL_TOL},
     "provenance":f"generated by build.py with params.py, seed {P.SEED}; values are the generating truth"}
json.dump(key,open(os.path.join(HERE,"answer_key.json"),"w"),indent=1)
print(f"answer key: {3*len(P.TARGETS)} concentrations, pass bar {P.CW_MIN_PASS}")
