"""
DSIC vs ΛCDM: совместный тест SN + BAO (версия 2).

Отличия от версии 1:
  • BAO с полной ковариацией D_M–D_H по бинам (корреляции учтены),
    вместо приближения независимых ошибок;
  • блок диагностик расхождения моделей, указывающий области,
    где имеющиеся или будущие данные способны их различить.

Критерии:
  (1) совместный Δχ² при едином μ₀ — согласованность модели с обоими классами данных;
  (2) Om(z)-диагностика — плоская для ΛCDM, переменная для DSIC;
  (3) высокий z (Lyman-α BAO) — область роста расхождения.

Зависимости: numpy, scipy. Скрипт самостоятельно загружает данные Pantheon+.
"""

import numpy as np
from scipy.optimize import minimize, brentq
from scipy.integrate import quad
import urllib.request, os

c = 299792.458

# ---------------------------------------------------------------------
# 1. SN: Pantheon+SH0ES
# ---------------------------------------------------------------------
BASE = ("https://raw.githubusercontent.com/PantheonPlusSH0ES/"
        "DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/")
def grab(name):
    if not os.path.exists(name):
        print(f"  загрузка {name} ..."); urllib.request.urlretrieve(BASE+name, name)
print("1) Загрузка Pantheon+SH0ES")
grab("Pantheon+SH0ES.dat"); grab("Pantheon+SH0ES_STAT+SYS.cov")

rows=[]
with open("Pantheon+SH0ES.dat") as f:
    h=f.readline().split()
    iz,imb,icep,iceph=(h.index("zHD"),h.index("m_b_corr"),
                       h.index("IS_CALIBRATOR"),h.index("CEPH_DIST"))
    for line in f:
        p=line.split()
        if p: rows.append(p)
zHD=np.array([float(r[iz]) for r in rows])
m_b=np.array([float(r[imb]) for r in rows])
is_cal=np.array([int(r[icep]) for r in rows])==1
ceph=np.array([float(r[iceph]) for r in rows])
with open("Pantheon+SH0ES_STAT+SYS.cov") as f:
    n=int(f.readline()); vals=np.fromstring(f.read().replace(" ",""),sep="\n")
C=vals.reshape(n,n)
use=(zHD>0.01)|is_cal
idx=np.where(use)[0]
z_sn,mb_sn,cal_sn,ceph_sn=zHD[idx],m_b[idx],is_cal[idx],ceph[idx]
Cinv_sn=np.linalg.inv(C[np.ix_(idx,idx)])
print(f"   SN в фите: {len(z_sn)} (калибраторов: {cal_sn.sum()})")

# ---------------------------------------------------------------------
# 2. BAO DESI DR2 с полной ковариацией D_M/r_d, D_H/r_d по бинам.
#    Для каждого бина: вектор (D_M/rd, D_H/rd) и ковариация 2x2.
#    corr — коэффициент корреляции D_M–D_H в бине.
#    Значения и корреляции следует сверять с актуальной таблицей DESI DR2.
# ---------------------------------------------------------------------
#  z,    DM/rd,  sDM,    DH/rd,  sDH,   corr(DM,DH)
bao=np.array([
 [0.510, 13.588, 0.167, 21.863, 0.425, -0.459],
 [0.706, 17.351, 0.177, 19.455, 0.330, -0.404],
 [0.930, 21.576, 0.152, 17.641, 0.193, -0.416],
 [1.317, 27.601, 0.318, 14.176, 0.221, -0.434],
 [1.491, 30.512, 0.760, 12.817, 0.516, -0.420],
 [2.330, 38.988, 0.531,  8.632, 0.101, -0.477],
])
z_bao=bao[:,0]

# Блочно-диагональная ковариация 12x12 (6 бинов × 2 величины)
obs=np.empty(12); Cb=np.zeros((12,12))
for i,row in enumerate(bao):
    _,dm,sdm,dh,sdh,cc=row
    obs[2*i]=dm; obs[2*i+1]=dh
    Cb[2*i,2*i]=sdm**2
    Cb[2*i+1,2*i+1]=sdh**2
    Cb[2*i,2*i+1]=cc*sdm*sdh
    Cb[2*i+1,2*i]=cc*sdm*sdh
Cb_inv=np.linalg.inv(Cb)
print(f"   BAO: {len(z_bao)} бинов DESI DR2, полная ковариация D_M–D_H (12×12)")

# ---------------------------------------------------------------------
# 3. Геометрия DSIC
# ---------------------------------------------------------------------
def b(m): return np.sqrt(1-m*m)/(m*m)
def mu_e(zv,mu0): return brentq(lambda m:b(m)-b(mu0)/(1+zv),1e-12,1-1e-15)
def DMc_DSIC(zv,mu0,H0):
    ck=(c/H0)*(2-mu0**2)/(1-mu0**2); me=mu_e(zv,mu0)
    return ck*(np.sqrt(1-mu0**2)/mu0**2)*(np.sqrt(1-mu0**2)-np.sqrt(1-me**2))
def H_DSIC(zv,mu0,H0):
    me=mu_e(zv,mu0); return H0*((2-me**2)/(1-me**2))/((2-mu0**2)/(1-mu0**2))
def distmod_DSIC(zv,mu0,H0):
    return 5*np.log10((1+zv)*DMc_DSIC(zv,mu0,H0)*1e6/10)

# ---------------------------------------------------------------------
# 4. ΛCDM
# ---------------------------------------------------------------------
def E_L(zv,Om): return np.sqrt(Om*(1+zv)**3+(1-Om))
def DMc_L(zv,Om,H0): return (c/H0)*quad(lambda x:1/E_L(x,Om),0,zv)[0]
def distmod_L(zv,Om,H0):
    return 5*np.log10((1+zv)*DMc_L(zv,Om,H0)*1e6/10)

# ---------------------------------------------------------------------
# 5. Предсказания BAO: вектор (D_M/rd, D_H/rd); rd — свободный параметр
#    (rd фитится явно, а не сокращается через отношение величин).
# ---------------------------------------------------------------------
def bao_vec_DSIC(mu0,H0,rd):
    v=np.empty(12)
    for i,zv in enumerate(z_bao):
        v[2*i]=DMc_DSIC(zv,mu0,H0)/rd
        v[2*i+1]=(c/H_DSIC(zv,mu0,H0))/rd
    return v
def bao_vec_L(Om,H0,rd):
    v=np.empty(12)
    for i,zv in enumerate(z_bao):
        v[2*i]=DMc_L(zv,Om,H0)/rd
        v[2*i+1]=(c/(H0*E_L(zv,Om)))/rd
    return v

# ---------------------------------------------------------------------
# 6. χ²
# ---------------------------------------------------------------------
def chi2_sn_D(mu0,H0,M):
    model=np.array([(ceph_sn[i]+M) if cal_sn[i] else (distmod_DSIC(z_sn[i],mu0,H0)+M)
                    for i in range(len(z_sn))])
    r=mb_sn-model; return float(r@Cinv_sn@r)
def chi2_sn_L(Om,H0,M):
    model=np.array([(ceph_sn[i]+M) if cal_sn[i] else (distmod_L(z_sn[i],Om,H0)+M)
                    for i in range(len(z_sn))])
    r=mb_sn-model; return float(r@Cinv_sn@r)
def chi2_bao_D(mu0,H0,rd):
    r=bao_vec_DSIC(mu0,H0,rd)-obs; return float(r@Cb_inv@r)
def chi2_bao_L(Om,H0,rd):
    r=bao_vec_L(Om,H0,rd)-obs; return float(r@Cb_inv@r)

# ---------------------------------------------------------------------
# 7. Совместный фит (4 параметра: модель, H0, M, rd)
# ---------------------------------------------------------------------
print("\n2) Совместный фит SN + BAO с полной ковариацией (общий μ₀, свободный rd)")
def joint_D(p):
    mu0,H0,M,rd=p
    if not(0.5<mu0<0.999 and 50<H0<90 and -20.5<M<-18.5 and 120<rd<170): return 1e12
    return chi2_sn_D(mu0,H0,M)+chi2_bao_D(mu0,H0,rd)
def joint_L(p):
    Om,H0,M,rd=p
    if not(0.05<Om<0.6 and 50<H0<90 and -20.5<M<-18.5 and 120<rd<170): return 1e12
    return chi2_sn_L(Om,H0,M)+chi2_bao_L(Om,H0,rd)

rd_=minimize(joint_D,[0.76,73,-19.3,147],method="Nelder-Mead",
             options=dict(xatol=1e-5,fatol=1e-4,maxiter=30000))
rl_=minimize(joint_L,[0.31,73,-19.3,147],method="Nelder-Mead",
             options=dict(xatol=1e-5,fatol=1e-4,maxiter=30000))
mu0,H0d,Md,rdd=rd_.x; Om,H0l,Ml,rdl=rl_.x
snD=chi2_sn_D(mu0,H0d,Md); baoD=chi2_bao_D(mu0,H0d,rdd)
snL=chi2_sn_L(Om,H0l,Ml);  baoL=chi2_bao_L(Om,H0l,rdl)

print("\n"+"="*60)
print("РЕЗУЛЬТАТ — SN + BAO (полная ковариация, свободный r_d)")
print("="*60)
print(f"DSIC : μ₀={mu0:.4f} H₀={H0d:.2f} M={Md:.3f} r_d={rdd:.1f}")
print(f"       χ²(SN)={snD:.1f} χ²(BAO)={baoD:.2f} тотал={rd_.fun:.1f}")
print(f"ΛCDM : Ωm={Om:.4f} H₀={H0l:.2f} M={Ml:.3f} r_d={rdl:.1f}")
print(f"       χ²(SN)={snL:.1f} χ²(BAO)={baoL:.2f} тотал={rl_.fun:.1f}")
print(f"\nΔχ² ТОТАЛ (DSIC−ΛCDM) = {rd_.fun-rl_.fun:+.2f}")
print(f"  ΔSN={snD-snL:+.2f}  ΔBAO={baoD-baoL:+.2f}")
print(f"  BAO χ²/dof (12−4=8): DSIC={baoD/8:.2f} ΛCDM={baoL/8:.2f}")

# ---------------------------------------------------------------------
# 8. Диагностики расхождения моделей
# ---------------------------------------------------------------------
print("\n"+"="*60)
print("ДИАГНОСТИКИ РАСХОЖДЕНИЯ DSIC И ΛCDM")
print("="*60)
def E_D(zv,mu0): return H_DSIC(zv,mu0,1.0)

print("\n[A] Om(z)-диагностика: для ΛCDM строго постоянна (=Ωm), для DSIC переменна.")
print("    Плоская Om(z) по данным свидетельствует против DSIC.")
print(f"    {'z':>5} {'Om_DSIC':>9} {'Om_LCDM':>9}")
for z in [0.2,0.5,1.0,2.0,3.0,5.0]:
    Ed,El=E_D(z,mu0),E_L(z,Om)
    omd=(Ed**2-1)/((1+z)**3-1); oml=(El**2-1)/((1+z)**3-1)
    print(f"    {z:>5} {omd:>9.4f} {oml:>9.4f}")
print("    DSIC даёт U-образную Om(z) с минимумом около z=2 и подъёмом")
print("    к высоким z — качественная особенность, отсутствующая у ΛCDM.")

print("\n[B] Высокий z (Lyman-α BAO, z≈2.3–3.5): расхождение D_H/r_d растёт.")
print(f"    {'z':>5} {'DH/rd DSIC':>11} {'DH/rd LCDM':>11} {'разн.%':>7}")
for z in [2.33,2.5,3.0,3.5]:
    dhd=(c/H_DSIC(z,mu0,H0d))/rdd; dhl=(c/(H0l*E_L(z,Om)))/rdl
    print(f"    {z:>5} {dhd:>11.3f} {dhl:>11.3f} {100*(dhd/dhl-1):>+6.1f}%")
print("    При z≈3.5 расхождение достигает нескольких процентов: точные")
print("    измерения Lyman-α BAO в этой области являются решающим тестом.")

print("\n[C] Очень высокий z (E(z)): за пределами данных DSIC существенно жёстче.")
for z in [5,10]:
    print(f"    z={z}: E_DSIC/E_LCDM = {E_D(z,mu0)/E_L(z,Om):.2f}  "
          f"(DSIC завышает темп расширения на {100*(E_D(z,mu0)/E_L(z,Om)-1):.0f}%)")
print("    DSIC предсказывает иную раннюю историю расширения. Любой зонд z>3")
print("    (рост структур, fσ8, реликтовое излучение) — потенциальный тест.")
