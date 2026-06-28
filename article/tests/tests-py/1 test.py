"""
DSIC vs ΛCDM: проверка на каталоге Pantheon+ (1701 SN Ia)
с ковариационной матрицей STAT+SYS.

Скрипт загружает каталог Pantheon+SH0ES и ковариационную матрицу,
выполняет космологический фит двух моделей (DSIC и плоская ΛCDM)
и сравнивает их по χ². Зависимости: numpy, scipy, matplotlib.
"""

import numpy as np
from scipy.optimize import brentq, minimize
import urllib.request, os

# ---------------------------------------------------------------------
# 1. Загрузка данных
# ---------------------------------------------------------------------
BASE = ("https://raw.githubusercontent.com/PantheonPlusSH0ES/"
        "DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/")
DAT = "Pantheon+SH0ES.dat"
COV = "Pantheon+SH0ES_STAT+SYS.cov"

def grab(name):
    if os.path.exists(name):
        print(f"  файл уже присутствует локально: {name}")
        return
    print(f"  загрузка {name} ...")
    urllib.request.urlretrieve(BASE + name, name)
    print(f"  готово: {os.path.getsize(name)/1e6:.1f} МБ")

print("1) Загрузка данных Pantheon+")
try:
    grab(DAT); grab(COV)
except Exception as e:
    print("\nАвтоматическая загрузка не выполнена:", e)
    print("Загрузите файлы вручную из каталога репозитория")
    print("Pantheon+_Data/4_DISTANCES_AND_COVAR и поместите их")
    print("в рабочую директорию.")
    raise SystemExit

# ---------------------------------------------------------------------
# 2. Чтение каталога
# ---------------------------------------------------------------------
print("\n2) Чтение каталога")
rows = []
with open(DAT) as f:
    header = f.readline().split()
    iz   = header.index("zHD")
    imu  = header.index("MU_SH0ES") if "MU_SH0ES" in header else header.index("MU")
    icep = header.index("IS_CALIBRATOR") if "IS_CALIBRATOR" in header else None
    for line in f:
        p = line.split()
        if not p: continue
        rows.append(p)
zHD = np.array([float(r[iz])  for r in rows])
MU  = np.array([float(r[imu]) for r in rows])
is_cal = (np.array([int(r[icep]) for r in rows]) == 1) if icep is not None \
         else np.zeros(len(rows), bool)
N = len(zHD)
print(f"  прочитано {N} объектов")

# ---------------------------------------------------------------------
# 3. Чтение ковариационной матрицы
#    (формат: первая строка — размер, далее значения столбцом)
# ---------------------------------------------------------------------
print("\n3) Чтение ковариационной матрицы")
with open(COV) as f:
    n = int(f.readline())
    vals = np.fromstring(f.read().replace(" ", ""), sep="\n")
assert n == N, f"размер матрицы {n} не совпадает с числом SN {N}"
C = vals.reshape(n, n)
print(f"  матрица {n}x{n} загружена")

# ---------------------------------------------------------------------
# 4. Отбор объектов для космологического фита
#    Исключаются калибраторы (MU заменён цефеидным расстоянием)
#    и объекты с очень низким z (вклад пекулярных скоростей).
# ---------------------------------------------------------------------
mask = (~is_cal) & (zHD > 0.01)
idx = np.where(mask)[0]
z = zHD[idx]; mu = MU[idx]
Csel = C[np.ix_(idx, idx)]
Cinv = np.linalg.inv(Csel)
print(f"\n4) В фит включено {len(z)} объектов (z = {z.min():.3f}…{z.max():.3f})")

# ---------------------------------------------------------------------
# 5. Модели
# ---------------------------------------------------------------------
c = 299792.458  # км/с

# DSIC: b(mu) = sqrt(1-mu^2)/mu^2 ;  1+z = b(mu0)/b(mu_e)
def b(m): return np.sqrt(1 - m*m) / (m*m)
def mu_e_of_z(zz, mu0, b0):
    target = b0 / (1 + zz)
    return brentq(lambda m: b(m) - target, 1e-9, 1 - 1e-15)
def DM_DSIC(zz, mu0, H0):
    b0 = b(mu0)
    ck = (c / H0) * (2 - mu0**2) / (1 - mu0**2)
    pref = np.sqrt(1 - mu0**2) / mu0**2
    out = np.empty_like(zz)
    for i, zi in enumerate(zz):
        me = mu_e_of_z(zi, mu0, b0)
        chi = pref * (np.sqrt(1 - mu0**2) - np.sqrt(1 - me**2)) * ck
        out[i] = 5 * np.log10((1 + zi) * chi * 1e6 / 10)
    return out

# Плоская ΛCDM
from scipy.integrate import quad
def DM_LCDM(zz, Om, H0):
    out = np.empty_like(zz)
    for i, zi in enumerate(zz):
        I = quad(lambda x: 1/np.sqrt(Om*(1+x)**3 + (1-Om)), 0, zi)[0]
        out[i] = 5*np.log10((c/H0)*(1+zi)*I*1e6/10)
    return out

# χ² с полной ковариацией
def chi2(resid):
    return float(resid @ Cinv @ resid)

# ---------------------------------------------------------------------
# 6. Фиты
# ---------------------------------------------------------------------
print("\n5) Подгонка моделей (минимизация χ² с полной ковариацией)")

# DSIC: свободные параметры mu0 и H0
def nll_dsic(p):
    mu0, H0 = p
    if not (0.5 < mu0 < 0.999) or not (50 < H0 < 90): return 1e12
    return chi2(mu - DM_DSIC(z, mu0, H0))
r1 = minimize(nll_dsic, [0.76, 67.4], method="Nelder-Mead",
              options=dict(xatol=1e-4, fatol=1e-3))
mu0_fit, H0_d = r1.x; chi2_d = r1.fun

# ΛCDM: свободные параметры Om и H0
def nll_lcdm(p):
    Om, H0 = p
    if not (0.05 < Om < 0.6) or not (50 < H0 < 90): return 1e12
    return chi2(mu - DM_LCDM(z, Om, H0))
r2 = minimize(nll_lcdm, [0.3, 67.4], method="Nelder-Mead",
              options=dict(xatol=1e-4, fatol=1e-3))
Om_fit, H0_l = r2.x; chi2_l = r2.fun

dof_d = len(z) - 2
dof_l = len(z) - 2

# ---------------------------------------------------------------------
# 7. Результаты
# ---------------------------------------------------------------------
print("\n" + "="*55)
print("РЕЗУЛЬТАТЫ (полный Pantheon+, ковариация STAT+SYS)")
print("="*55)
print(f"DSIC : mu0 = {mu0_fit:.4f},  H0 = {H0_d:.2f}")
print(f"       chi2 = {chi2_d:.1f},  chi2/dof = {chi2_d/dof_d:.4f}")
print(f"ΛCDM : Om  = {Om_fit:.4f},  H0 = {H0_l:.2f}")
print(f"       chi2 = {chi2_l:.1f},  chi2/dof = {chi2_l/dof_l:.4f}")
print(f"\nΔχ² (DSIC − ΛCDM) = {chi2_d - chi2_l:+.2f}")
print("(положительное значение — предпочтение ΛCDM; |Δχ²|<~2 — модели неразличимы)")

# ---------------------------------------------------------------------
# 8. График остатков Хаббла
# ---------------------------------------------------------------------
print("\n6) Построение графика остатков Хаббла")
import matplotlib.pyplot as plt
order = np.argsort(z)
res_d = mu - DM_DSIC(z, mu0_fit, H0_d)
res_l = mu - DM_LCDM(z, Om_fit, H0_l)
plt.figure(figsize=(9,5))
plt.axhline(0, color="k", lw=0.8)

# Биннирование для наглядности
def binned(zz, rr, nb=18):
    bins = np.logspace(np.log10(zz.min()), np.log10(zz.max()), nb+1)
    bc, bm, be = [], [], []
    for i in range(nb):
        m = (zz>=bins[i])&(zz<bins[i+1])
        if m.sum()>2:
            bc.append(np.median(zz[m])); bm.append(np.mean(rr[m]))
            be.append(rr[m].std()/np.sqrt(m.sum()))
    return np.array(bc),np.array(bm),np.array(be)
for rr, lab, col in [(res_d,"DSIC","crimson"),(res_l,"ΛCDM","steelblue")]:
    bc,bm,be = binned(z,rr)
    plt.errorbar(bc,bm,yerr=be,fmt="o-",color=col,label=lab,capsize=3)
plt.xscale("log"); plt.xlabel("z"); plt.ylabel("μ_data − μ_model (binned)")
plt.title("Pantheon+: остатки относительно DSIC и ΛCDM")
plt.legend(); plt.tight_layout(); plt.savefig("residuals.png", dpi=130)
plt.show()
print("\nГотово. Совпадение обеих кривых с нулём означает неразличимость")
print("моделей на данных SN; отклонение кривой DSIC при больших z указывает")
print("на расхождение моделей в этой области.")
