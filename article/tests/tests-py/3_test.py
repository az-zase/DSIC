"""
Фит с калибраторами в качестве абсолютного якоря и явным параметром M
(абсолютная звёздная величина SN Ia).

Отличия от базового прогона:
  • калибраторы не исключаются, а используются как абсолютный якорь;
  • M фитится явно, отдельно от H0;
  • для калибраторов модель расстояния — цефеидный модуль (внешний якорь),
    для остальных объектов — космологический модуль расстояния модели;
  • H0 становится независимым выводом, а не поглотителем сдвига.

Требует предварительного запуска основного скрипта (файл 1):
из него используются C, DM_DSIC, DM_LCDM.
"""

import numpy as np
from scipy.optimize import minimize

# Чтение цефеидных расстояний калибраторов из каталога (столбец CEPH_DIST)
rows2 = []
with open("Pantheon+SH0ES.dat") as f:
    hdr = f.readline().split()
    icz   = hdr.index("zHD")
    icmu  = hdr.index("MU_SH0ES")
    iccal = hdr.index("IS_CALIBRATOR")
    iceph = hdr.index("CEPH_DIST")
    icmb  = hdr.index("m_b_corr") if "m_b_corr" in hdr else None
    for line in f:
        p = line.split()
        if p: rows2.append(p)

ceph = np.array([float(r[iceph]) for r in rows2])                  # цефеидный модуль расстояния
mb   = np.array([float(r[hdr.index("m_b_corr")]) for r in rows2])  # видимая величина SN
zall = np.array([float(r[icz]) for r in rows2])
cal  = np.array([int(r[iccal]) for r in rows2]) == 1

# Отбор: для космологии z>0.01, калибраторы включаются полностью
use = (zall > 0.01) | cal
idx2 = np.where(use)[0]
z2   = zall[idx2]
mb2  = mb[idx2]
cal2 = cal[idx2]
ceph2= ceph[idx2]
Csel2 = C[np.ix_(idx2, idx2)]
Cinv2 = np.linalg.inv(Csel2)
print(f"в фит включено {len(z2)} объектов, из них калибраторов: {cal2.sum()}")

def chi2_full(resid):
    return float(resid @ Cinv2 @ resid)

# Остатки: наблюдаемая m_b минус (модель расстояния + M).
# Для калибраторов модель расстояния — цефеидный модуль (внешний якорь);
# для остальных — космологический модуль расстояния модели.
def resid_DSIC(params):
    mu0, H0, M = params
    mod = np.empty_like(z2)
    mod[cal2] = ceph2[cal2] + M
    zc = z2[~cal2]
    mod[~cal2] = DM_DSIC(zc, mu0, H0) + M
    return mb2 - mod

def resid_LCDM(params):
    Om, H0, M = params
    mod = np.empty_like(z2)
    mod[cal2] = ceph2[cal2] + M
    zc = z2[~cal2]
    mod[~cal2] = DM_LCDM(zc, Om, H0) + M
    return mb2 - mod

# Фиты с тремя параметрами (модель + H0 + M)
def nll_d(p):
    mu0,H0,M = p
    if not(0.5<mu0<0.999 and 50<H0<90 and -20.5<M<-18.5): return 1e12
    return chi2_full(resid_DSIC(p))
rd = minimize(nll_d, [0.76,73.0,-19.3], method="Nelder-Mead",
              options=dict(xatol=1e-5,fatol=1e-4,maxiter=20000))

def nll_l(p):
    Om,H0,M = p
    if not(0.05<Om<0.6 and 50<H0<90 and -20.5<M<-18.5): return 1e12
    return chi2_full(resid_LCDM(p))
rl = minimize(nll_l, [0.33,73.0,-19.3], method="Nelder-Mead",
              options=dict(xatol=1e-5,fatol=1e-4,maxiter=20000))

dof = len(z2) - 3   # три параметра
k = 3               # для AIC/BIC
n_pts = len(z2)

print("\n"+"="*60)
print("Калибраторы-якорь, явный M, H0 как независимый вывод")
print("="*60)
print(f"DSIC : mu0={rd.x[0]:.4f}  H0={rd.x[1]:.2f}  M={rd.x[2]:.4f}")
print(f"       chi2={rd.fun:.1f}  chi2/dof={rd.fun/dof:.4f}")
print(f"ΛCDM : Om ={rl.x[0]:.4f}  H0={rl.x[1]:.2f}  M={rl.x[2]:.4f}")
print(f"       chi2={rl.fun:.1f}  chi2/dof={rl.fun/dof:.4f}")
print(f"\nΔχ² (DSIC−ΛCDM) = {rd.fun-rl.fun:+.2f}")

# Информационные критерии AIC / BIC
# (при одинаковом k разница равна Δχ², приводится явно)
def aic(chi2,k): return chi2 + 2*k
def bic(chi2,k,n): return chi2 + k*np.log(n)
print(f"\nΔAIC = {aic(rd.fun,k)-aic(rl.fun,k):+.2f}")
print(f"ΔBIC = {bic(rd.fun,k,n_pts)-bic(rl.fun,k,n_pts):+.2f}")
print("(|Δ|<2 — модели неразличимы; 2–6 — слабое предпочтение; >10 — сильное)")
print(f"\nH0 как независимый вывод:  DSIC={rd.x[1]:.2f}   ΛCDM={rl.x[1]:.2f}")
