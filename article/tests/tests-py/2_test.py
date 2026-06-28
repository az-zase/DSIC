"""
Дополнительный прогон с матрицей STAT-ONLY (только статистические ошибки,
без систематик).

Цель: оценить устойчивость результата. Если χ²/dof и оценки параметров
существенно не изменяются относительно прогона STAT+SYS, результат не
обусловлен моделью систематик.

Требует предварительного запуска основного скрипта (файл 1):
из него используются idx, z, mu, C, DM_DSIC, DM_LCDM, mu0_fit, Om_fit,
H0_d, H0_l, chi2_d, chi2_l.
"""

import numpy as np, urllib.request, os
from scipy.optimize import minimize

COV2 = "Pantheon+SH0ES_STATONLY.cov"
BASE = ("https://raw.githubusercontent.com/PantheonPlusSH0ES/"
        "DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/")

# Загрузка матрицы STAT-ONLY
if not os.path.exists(COV2):
    print(f"загрузка {COV2} ...")
    urllib.request.urlretrieve(BASE + COV2, COV2)
    print(f"готово: {os.path.getsize(COV2)/1e6:.1f} МБ")
else:
    print(f"файл уже присутствует: {COV2}")

# Чтение (формат: первая строка — размер, далее значения столбцом)
with open(COV2) as f:
    n2 = int(f.readline())
    vals2 = np.fromstring(f.read().replace(" ", ""), sep="\n")
Cstat = vals2.reshape(n2, n2)
print(f"матрица STAT-ONLY {n2}x{n2} загружена")

# Тот же отбор объектов, что и в основном прогоне
Cstat_sel = Cstat[np.ix_(idx, idx)]
Cinv_stat = np.linalg.inv(Cstat_sel)

def chi2_stat(resid):
    return float(resid @ Cinv_stat @ resid)

# Фиты с ковариацией STAT-ONLY
print("\nПодгонка с матрицей STAT-ONLY ...")

def nll_dsic_s(p):
    mu0, H0 = p
    if not (0.5 < mu0 < 0.999) or not (50 < H0 < 90): return 1e12
    return chi2_stat(mu - DM_DSIC(z, mu0, H0))
rs1 = minimize(nll_dsic_s, [0.76, 73.0], method="Nelder-Mead",
               options=dict(xatol=1e-4, fatol=1e-3))

def nll_lcdm_s(p):
    Om, H0 = p
    if not (0.05 < Om < 0.6) or not (50 < H0 < 90): return 1e12
    return chi2_stat(mu - DM_LCDM(z, Om, H0))
rs2 = minimize(nll_lcdm_s, [0.33, 73.0], method="Nelder-Mead",
               options=dict(xatol=1e-4, fatol=1e-3))

dof = len(z) - 2
print("\n" + "="*55)
print("СРАВНЕНИЕ: STAT+SYS vs STAT-ONLY")
print("="*55)
print(f"{'':12}{'mu0/Om':>10}{'H0':>9}{'chi2/dof':>11}")
print("-- STAT+SYS (основной прогон) --")
print(f"{'DSIC':12}{mu0_fit:>10.4f}{H0_d:>9.2f}{chi2_d/dof:>11.4f}")
print(f"{'LCDM':12}{Om_fit:>10.4f}{H0_l:>9.2f}{chi2_l/dof:>11.4f}")
print("-- STAT-ONLY (данный прогон) --")
print(f"{'DSIC':12}{rs1.x[0]:>10.4f}{rs1.x[1]:>9.2f}{rs1.fun/dof:>11.4f}")
print(f"{'LCDM':12}{rs2.x[0]:>10.4f}{rs2.x[1]:>9.2f}{rs2.fun/dof:>11.4f}")
print(f"\nΔχ² (DSIC−LCDM), STAT-ONLY = {rs1.fun - rs2.fun:+.2f}")
print("\nИнтерпретация:")
print("• близость mu0 и формы к значениям STAT+SYS указывает на устойчивость;")
print("• χ²/dof для STAT-ONLY обычно выше (систематики увеличивают заявленные")
print("  ошибки; при их исключении разброс относительно модели растёт);")
print("• ключевой критерий — сопоставимость DSIC и ΛCDM (малое |Δχ²|).")
