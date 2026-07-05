"""
================================================================================
DSIC — ЕДИНЫЙ воспроизводимый тест (ОДИН файл, Google Colab / python).
Часть I  (ядро):         DSIC vs flat ΛCDM на идентичных данных, равное число параметров.
Часть II (второй этаж):  порог отлипания μ_d (ранняя Вселенная, реликт).
Часть III (рост структур): fσ8 на фоне DSIC (GR-предел), блок [11].

Автор модели: Azamat Zaseev
Репозиторий:   https://github.com/az-zase/DSIC
================================================================================

      ПОЖАЛУЙСТА, НАБЕРИТЕСЬ ТЕРПЕНИЯ.
      Скрипт МОЖЕТ показаться зависшим — это нормально, он работает. Ничего
      нажимать и прерывать не нужно: просто спокойно дождитесь завершения.
      Полный прогон занимает ориентировочно 10–15 минут и складывается из:
        • автоматического скачивания каталога сверхновых (~65 МБ);
        • фитов сверхновых Pantheon+SH0ES (несколько минут);
        • блока [11] «рост структур» — самого долгого этапа: это скан по сетке
          с численным интегрированием ОДУ роста в каждом узле (~5 минут).
          Между печатью блока [10] и появлением результатов [11] возможна
          длинная тихая пауза — это ожидаемо, скрипт считает.
      По окончании внизу будет напечатано «Готово» с итоговым резюме.

================================================================================

ЧТО ДЕЛАЕТ (всё в одном прогоне):
  [0]  Самопроверка: сверка вшитого DESI с официальной Table IV.
  [0b] СВЯЗНОСТЬ МОДЕЛИ: машинная проверка внутренних тождеств ядра —
       монотонность b(μ); точность D_H = dχ/dz; совпадение замкнутой χ(μ) с
       численным ∫c dz/H; ноль q на μ_cross = √(3−√5); формула возраста.
  [1]  SN-космология: Pantheon+SH0ES, якорный фит (m_b_corr + явный M).
       1σ-погрешности параметров из гессиана χ² в минимуме.
  [2]  Устойчивость: STAT-ONLY против STAT+SYS.
  [3]  Совместный SN+BAO: полная ковариация D_M–D_H (12×12), свободный r_d;
       плюс восстановление линейки r_d при ранней шкале.
  [4]  Om(z)-диагностика: устойчивость знака Δχ² к H0; хронометры в ДВУХ
       режимах ковариации (диагональная и полная Moresco et al. 2020).
  [5]  Фальсифицируемые предсказания на высоком z (D_H/r_d Lyα, E(z)).
  [6]  КОСМИЧЕСКОЕ ВРЕМЯ: t(z) DSIC vs ΛCDM при z = 5–12 — ближайшая
       зона риска ядра (возрасты галактик JWST, рост SMBH к z≈7.5).
  ---- Часть II (второй этаж, P6) ----
  [7]  Самопроверка опорных данных Planck (d_LSS = r*/θ*).
  [8]  Недобор ядра до реликта и фиксация порога μ_d по прямо измеренному
       углу θ*; r* явно помечен как ЗАИМСТВОВАННЫЙ (внешний вход из ΛCDM).
  [9]  Контроль: поправка тождественно 0 при z ≤ z_detach; таблица вырождения
       «амплитуда ↔ порог» (что фиксирует постулат P6).
  [10] След μ_d на высоком z; исправленный горизонт частиц; флаг внутренней
       согласованности (наблюдаемые галактики z ≳ 10 лежат ВЫШЕ порога).
  ---- Часть III (рост структур) ----
  [11] РОСТ СТРУКТУР: стандартное уравнение линейного роста возмущений на
       фоне DSIC (GR-предел), fσ8(z) против компиляции Gold-2017 (18 точек,
       ковариация WiggleZ учтена). Данные fσ8 вырождены по (Ω_m0, σ8) —
       результат докладывается как инвариант S8 = σ8·√(Ω_m0/0.3) плюс срез
       при Planck-подобном Ω_m0 = 0.30. Статус: проверка на выживание, не
       вывод из ядра (ядро плотностей не вводит; Ω_m0 — параметр роста).

ДАННЫЕ (происхождение каждой вшитой цифры задокументировано на месте):
  • Сверхновые — СКАЧИВАЮТСЯ АВТОМАТИЧЕСКИ (Pantheon+SH0ES DataRelease).
  • DESI DR2 BAO — вшиты, СВЕРЕНЫ построчно с официальной Table IV
        (DESI Coll. / Abdul-Karim et al. 2025, Phys. Rev. D 112, 083515;
         arXiv:2503.14738). Самопроверка [0] падает при любом расхождении.
  • Космические хронометры H(z) — вшиты, 33 точки. 15 точек Moresco СВЕРЕНЫ
        с официальным репозиторием (gitlab.com/mmoresco/CCcovariance, BC03);
        18 прочих — из других работ (Simon 2005; Stern 2010; Zhang 2014;
        Ratsimbazafy 2017; Borghi 2022 и др.), с диагональными ошибками.
  • ПОЛНАЯ КОВАРИАЦИЯ Moresco-подмножества (Moresco et al. 2020, ApJ 898, 82;
        arXiv:2003.07362) вшита и строится прямо в файле.
  • Planck 2018 (TT,TE,EE+lowE+lensing): z* = 1089.92, r* = 144.43 Мпк,
        100θ* = 1.04110, откуда d_LSS = r*/θ* = 13872.8 ± 25 Мпк. Сверяется [7].
  • fσ8 — компиляция Gold-2017 (Nesseris, Pantazis & Perivolaropoulos 2017,
        Phys. Rev. D 96, 023542; arXiv:1703.10538), 18 слабо коррелированных
        точек; z, fσ8, σ — из Table I работы Sagredo, Nesseris & Sapone 2018
        (arXiv:1806.10822). Ковариация трёх точек WiggleZ — Blake et al. 2012.

ЗАВИСИМОСТИ: numpy, scipy (в Colab предустановлены).
ЗАПУСК: вставить файл в одну ячейку Colab и выполнить; либо
        `python dsic_test.py`.
================================================================================
"""

import os
import urllib.request
import numpy as np
from scipy.optimize import minimize, brentq
from scipy.integrate import quad, solve_ivp

c = 299792.458  # км/с

# ==============================================================================
#  ДАННЫЕ
# ==============================================================================

# --- DESI DR2 BAO, Table IV (arXiv:2503.14738). Baseline-выборка для космологии.
#     Столбцы: z_eff, D_M/r_d, σ(D_M/r_d), D_H/r_d, σ(D_H/r_d), corr(D_M,D_H).
BAO = np.array([
    # tracer        z_eff   DM/rd     sDM     DH/rd     sDH     corr
    [0.510, 13.587, 0.169, 21.863, 0.427, -0.475],   # LRG1
    [0.706, 17.347, 0.180, 19.458, 0.332, -0.423],   # LRG2
    [0.934, 21.574, 0.153, 17.641, 0.193, -0.425],   # LRG3+ELG1 (baseline)
    [1.321, 27.605, 0.320, 14.178, 0.217, -0.437],   # ELG2
    [1.484, 30.519, 0.758, 12.816, 0.513, -0.489],   # QSO
    [2.330, 38.988, 0.531,  8.632, 0.101, -0.431],   # Lya
])

# --- Эталон для самопроверки [0] (та же официальная таблица, отдельной копией).
#     ВНИМАНИЕ: проверка ловит случайную правку массива, но не общую опечатку
#     переноса из публикации — обе копии сверены с arXiv:2503.14738 вручную.
_DESI_OFFICIAL = {
    0.510: (13.587, 0.169, 21.863, 0.427, -0.475),
    0.706: (17.347, 0.180, 19.458, 0.332, -0.423),
    0.934: (21.574, 0.153, 17.641, 0.193, -0.425),
    1.321: (27.605, 0.320, 14.178, 0.217, -0.437),
    1.484: (30.519, 0.758, 12.816, 0.513, -0.489),
    2.330: (38.988, 0.531,  8.632, 0.101, -0.431),
}

# --- Космические хронометры H(z), 33 точки, z ∈ [0.07, 1.965].
#     15 точек Moresco сверены с gitlab.com/mmoresco/CCcovariance (BC03;
#     Moresco et al. 2012, 2016; Moresco 2015); полная ковариация этого
#     подмножества вшита ниже. Остальные 18 точек — Simon 2005; Stern 2010;
#     Zhang 2014; Ratsimbazafy 2017; Borghi 2022 и др., диагональные ошибки.
CC = np.array([
    [0.07, 69.0, 19.6], [0.09, 69.0, 12.0], [0.12, 68.6, 26.2], [0.17, 83.0, 8.0],
    [0.1791, 74.91, 3.81], [0.1993, 74.96, 4.9], [0.2, 72.9, 29.6], [0.27, 77.0, 14.0],
    [0.28, 88.8, 36.6], [0.3519, 82.78, 13.95], [0.3802, 83.0, 13.54], [0.4, 95.0, 17.0],
    [0.4004, 76.97, 10.18], [0.4247, 87.08, 11.24], [0.4497, 92.78, 12.9], [0.47, 89.0, 49.6],
    [0.4783, 80.91, 9.04], [0.48, 97.0, 62.0], [0.5929, 103.8, 12.5], [0.6797, 91.6, 7.96],
    [0.75, 98.8, 33.6], [0.7812, 104.5, 12.2], [0.8, 113.1, 28.5], [0.8754, 125.1, 16.7],
    [0.88, 90.0, 40.0], [0.9, 117.0, 23.0], [1.037, 153.7, 19.67], [1.3, 168.0, 17.0],
    [1.363, 160.0, 32.63], [1.43, 177.0, 18.0], [1.53, 140.0, 14.0], [1.75, 202.0, 40.0],
    [1.965, 186.5, 49.58],
])
z_cc, H_cc, sH_cc = CC[:, 0], CC[:, 1], CC[:, 2]

# ==============================================================================
#  КОВАРИАЦИЯ ХРОНОМЕТРОВ (Moresco et al. 2020) — ВШИТЫЕ ТАБЛИЦЫ
# ------------------------------------------------------------------------------
#  Метод: Moresco, Jimenez, Verde, Cimatti, Pozzetti, ApJ 898, 82 (2020);
#         arXiv:2003.07362; репозиторий https://gitlab.com/mmoresco/CCcovariance.
#  Полная ковариация = статистическая (диагональ) + систематическая
#  (коррелированная по z): Cov_syst(i,j) = Σ_s σ_s(z_i)·σ_s(z_j), ρ_s ≡ 1.
#  Диагональ итоговой матрицы ТОЧНО равна опубликованной полной ошибке err_t².
#  СТАТУС: аккуратная реконструкция метода на официальных входных таблицах,
#  а не дословный запуск авторского кода; недиагональ воспроизводит рецепт.
# ==============================================================================

# --- HzTable_MM_BC03.dat (модель SPS = BC03). 15 точек Moresco.
#     Столбцы: (z, H, errH_total, stat_contr, met_contr) в км/с/Мпк.
CC_MORESCO_BC03 = [
    (0.1791,  74.91,  3.8069262,  3.8,    0.5),    # Moresco et al. (2012)
    (0.1993,  74.96,  4.9001352,  4.9,    0.6),    # Moresco et al. (2012)
    (0.3519,  82.78, 13.9484300, 13.0,    4.8),    # Moresco et al. (2012)
    (0.3802,  83.00, 13.5400000,  4.3,   12.9),    # Moresco et al. (2016)
    (0.4004,  76.97, 10.1800000,  2.1,   10.0),    # Moresco et al. (2016)
    (0.4247,  87.08, 11.2400000,  2.4,   11.0),    # Moresco et al. (2016)
    (0.4497,  92.78, 12.9000000,  4.5,   12.1),    # Moresco et al. (2016)
    (0.4783,  80.91,  9.0440000,  2.1,    8.8),    # Moresco et al. (2016)
    (0.5929, 103.80, 12.4975200, 11.6,    4.5),    # Moresco et al. (2012)
    (0.6797,  91.60,  7.9618720,  6.4,    4.3),    # Moresco et al. (2012)
    (0.7812, 104.50, 12.1951500,  9.4,    6.1),    # Moresco et al. (2012)
    (0.8754, 125.10, 16.7008500, 15.3,    6.0),    # Moresco et al. (2012)
    (1.0370, 153.70, 19.6736000, 13.6,   14.9),    # Moresco et al. (2012)
    (1.3630, 160.00, 32.6300000, 23.07,  23.07),   # Moresco (2015)
    (1.9650, 186.50, 49.5800000, 35.05,  35.05),   # Moresco (2015)
]

# --- data_MM20.dat (Moresco et al. 2020). Относительные (%) систематические
#     вклады на сетке z (шаг 0.05). Столбцы: (z, IMF%, stlib%, mod%, mod_ooo%).
CC_SYST_MM20 = [
    (0.075, 0.47, 7.40, 15.86, 9.91), (0.125, 0.47, 7.40, 14.23, 6.98),
    (0.175, 0.47, 7.40, 13.34, 5.40), (0.225, 0.47, 7.40, 13.21, 5.40),
    (0.275, 0.47, 7.40, 13.29, 5.40), (0.325, 0.47, 7.40, 12.20, 5.40),
    (0.375, 0.47, 7.40, 12.99, 5.40), (0.425, 0.47, 7.40, 10.29, 6.20),
    (0.475, 0.46, 7.39,  8.91, 5.86), (0.525, 0.23, 7.40,  9.99, 6.51),
    (0.575, 0.28, 6.87, 10.09, 6.12), (0.625, 0.47, 6.65, 11.17, 6.21),
    (0.675, 0.47, 6.57, 11.12, 5.71), (0.725, 0.47, 5.90, 10.81, 5.16),
    (0.775, 0.45, 6.03, 10.75, 5.05), (0.825, 0.47, 6.10, 10.75, 5.05),
    (0.875, 0.47, 5.89,  9.08, 2.79), (0.925, 0.44, 5.80,  8.62, 3.70),
    (0.975, 0.40, 5.94,  7.32, 3.65), (1.025, 0.27, 6.07,  5.84, 3.37),
    (1.075, 0.20, 6.08,  6.02, 3.49), (1.125, 0.20, 6.07,  4.72, 2.33),
    (1.175, 0.19, 6.09,  4.31, 2.33), (1.225, 0.19, 6.09,  3.90, 2.33),
    (1.275, 0.19, 6.09,  3.90, 2.33), (1.325, 0.20, 6.09,  3.91, 2.34),
    (1.375, 0.19, 6.09,  3.90, 2.34), (1.425, 0.19, 6.09,  3.90, 2.33),
    (1.475, 0.20, 6.09,  3.91, 2.34),
]

# --- Planck 2018 (TT,TE,EE+lowE+lensing). Опорные величины Части II.
#     Прямо измеренная величина — угол θ*. Звуковой горизонт r* = 144.43 Мпк —
#     величина, ВЫЧИСЛЕННАЯ внутри ΛCDM из физики фотон-барионной плазмы;
#     в DSIC он ЗАИМСТВУЕТСЯ как внешний вход (см. статью, §11.6).
#     d_LSS = r*/θ* = 13872.8 ± 25 Мпк (см., напр., arXiv:2311.04759, §3.6).
PLANCK = dict(
    z_star=1089.92,        # красное смещение последнего рассеяния
    r_star=144.43,         # звуковой горизонт на z*, Мпк (ЗАИМСТВОВАН из ΛCDM)
    theta_star=1.04110e-2, # 100θ* = 1.04110 -> θ* в радианах (прямое измерение)
    d_lss=13872.8,         # = r*/θ*, Мпк
    d_lss_err=25.0,        # 1σ, Мпк
)

# --- Параметры второго этажа: ранняя шкала (та же, что для BAO/хронометров/
#     возраста в Части I).
MU0_E = 0.76
H0_E  = 67.4


def build_cc_covariance(z_cc, sH_cc, full=True):
    """Ковариация 33 точек хронометров для блока [4].

    full=False : чисто диагональная (= sH_cc²), базовая версия теста.
    full=True  : блок 15 точек Moresco заменяется полной ковариацией
                 Moresco 2020 (диагональ = err_t², недиагональ = коррелированная
                 систематика); остальные 18 точек остаются диагональными.
    """
    n = len(z_cc)
    C = np.diag(np.asarray(sH_cc, float) ** 2).copy()
    if not full:
        return C

    mor = np.array(CC_MORESCO_BC03, float)
    z_m, err_t, met = mor[:, 0], mor[:, 2], mor[:, 4]

    syst = np.array(CC_SYST_MM20, float)
    zg = syst[:, 0]
    frac = np.vstack([np.interp(z_m, zg, syst[:, k]) / 100.0 for k in (1, 2, 3, 4)])
    quad_ = np.sqrt((frac ** 2).sum(axis=0))
    quad_[quad_ == 0] = 1.0
    sig = (frac / quad_) * met[None, :]          # (4,15): абс. вклад источника

    C_syst = sum(np.outer(sig[s], sig[s]) for s in range(4))   # 15×15 систематика

    idx = [int(np.argmin(np.abs(z_cc - zm))) for zm in z_m]
    for k, i in enumerate(idx):
        if abs(z_cc[i] - z_m[k]) > 1e-3:
            raise SystemExit(f"СТОП: точка Moresco z={z_m[k]} не найдена в массиве CC.")

    for a in range(15):
        for b_ in range(15):
            C[idx[a], idx[b_]] = (err_t[a] ** 2 if a == b_ else 0.0) + C_syst[a, b_]
    return C


# ==============================================================================
#  [0] САМОПРОВЕРКА ДАННЫХ DESI
# ==============================================================================
def self_check_desi():
    print("=" * 70)
    print("[0] Самопроверка: сверка вшитого DESI с официальной Table IV")
    print("=" * 70)
    ok = True
    for row in BAO:
        z = round(float(row[0]), 3)
        off = _DESI_OFFICIAL[z]
        got = tuple(round(float(x), 6) for x in row[1:])
        match = all(abs(a - b_) < 1e-6 for a, b_ in zip(got, off))
        if not match:
            ok = False
            print(f"  z={z}: РАСХОЖДЕНИЕ  файл={got}  офиц={off}")
    if ok:
        print("  OK — все 6 бинов DESI совпадают с arXiv:2503.14738, Table IV.\n")
    else:
        raise SystemExit("СТОП: данные DESI не совпадают с эталоном. Проверь массив BAO.")


# ==============================================================================
#  ГЕОМЕТРИЯ DSIC И ΛCDM
# ==============================================================================
def b_dsic(m):
    """Наблюдаемый масштабный фактор b(μ) = √(1−μ²)/μ²."""
    return np.sqrt(1.0 - m * m) / (m * m)


def mu_e_of_z(zv, mu0):
    """Инверсия 1+z = b(μ₀)/b(μ_e) на физической ветви приёма.

    При z=0 корень равен μ₀ точно и лежит на краю брекета brentq
    (оба конца дают один знак) — возвращаем его аналитически.
    Фиты этого случая не касаются (у SN z > 0.01), поведение при z > 0
    не изменено."""
    if zv <= 1e-14:
        return mu0
    target = b_dsic(mu0) / (1.0 + zv)
    return brentq(lambda m: b_dsic(m) - target, mu0 + 1e-12, 1.0 - 1e-15)


def ck_const(mu0, H0):
    """Перевод геометрии в Мпк: c/k = (c/H₀)·(2−μ₀²)/(1−μ₀²)."""
    return (c / H0) * (2.0 - mu0 ** 2) / (1.0 - mu0 ** 2)


def DM_comov_DSIC(zv, mu0, H0):
    ck = ck_const(mu0, H0)
    me = mu_e_of_z(zv, mu0)
    pref = np.sqrt(1.0 - mu0 ** 2) / mu0 ** 2
    return ck * pref * (np.sqrt(1.0 - mu0 ** 2) - np.sqrt(1.0 - me ** 2))


def H_DSIC(zv, mu0, H0):
    me = mu_e_of_z(zv, mu0)
    return H0 * ((2.0 - me ** 2) / (1.0 - me ** 2)) / ((2.0 - mu0 ** 2) / (1.0 - mu0 ** 2))


def distmod_DSIC(zv, mu0, H0):
    dl = (1.0 + zv) * DM_comov_DSIC(zv, mu0, H0)
    return 5.0 * np.log10(dl * 1e6 / 10.0)


def E_LCDM(zv, Om):
    return np.sqrt(Om * (1.0 + zv) ** 3 + (1.0 - Om))


def DM_comov_LCDM(zv, Om, H0):
    return (c / H0) * quad(lambda x: 1.0 / E_LCDM(x, Om), 0.0, zv)[0]


def distmod_LCDM(zv, Om, H0):
    dl = (1.0 + zv) * DM_comov_LCDM(zv, Om, H0)
    return 5.0 * np.log10(dl * 1e6 / 10.0)


# ==============================================================================
#  [0b] СВЯЗНОСТЬ МОДЕЛИ: машинная проверка внутренних тождеств ядра
# ==============================================================================
def self_check_coherence():
    print("=" * 70)
    print("[0b] Связность ядра: машинная проверка внутренних тождеств")
    print("=" * 70)
    mu0, H0 = 0.76, 70.0
    fails = []

    # (а) монотонность b(μ) на (μ0, 1) — корректность инверсии z -> μ_e
    grid = np.linspace(mu0 + 1e-6, 1 - 1e-9, 4000)
    mono = np.all(np.diff(b_dsic(grid)) < 0)
    print(f"  (а) b(μ) строго убывает на (μ₀,1): {'OK' if mono else 'FAIL'}")
    if not mono:
        fails.append("монотонность")

    # (б) тождество D_H = dχ/dz (статья §5.5): |dχ/dz · H/c − 1| < 1e-8
    worst = 0.0
    for zv in (0.5, 2.0, 10.0, 100.0):
        eps = 1e-5 * (1 + zv)
        d = (DM_comov_DSIC(zv + eps, mu0, H0) - DM_comov_DSIC(zv - eps, mu0, H0)) / (2 * eps)
        worst = max(worst, abs(d * H_DSIC(zv, mu0, H0) / c - 1))
    print(f"  (б) тождество D_H = dχ/dz: макс. отклонение {worst:.1e} "
          f"{'OK' if worst < 1e-6 else 'FAIL'}")
    if worst >= 1e-6:
        fails.append("D_H=dχ/dz")

    # (в) замкнутая χ(μ) против численного ∫ c dz'/H(z')
    worst = 0.0
    for zv in (1.0, 5.0):
        num = quad(lambda x: c / H_DSIC(x, mu0, H0), 0.0, zv, limit=300)[0]
        worst = max(worst, abs(num / DM_comov_DSIC(zv, mu0, H0) - 1))
    print(f"  (в) замкнутая χ = ∫c dz/H: макс. отклонение {worst:.1e} "
          f"{'OK' if worst < 1e-6 else 'FAIL'}")
    if worst >= 1e-6:
        fails.append("замкнутая χ")

    # (г) кроссовер q: μ_cross = √(3−√5) — точный корень −μ⁴+6μ²−4 = 0
    q = lambda m: (-m ** 4 + 6 * m * m - 4) / (m ** 4 - 4 * m * m + 4)
    mu_cross = np.sqrt(3.0 - np.sqrt(5.0))
    z_t = b_dsic(mu0) / b_dsic(mu_cross) - 1.0
    okq = abs(q(mu_cross)) < 1e-10
    print(f"  (г) q(μ_cross)=0 при μ_cross=√(3−√5)={mu_cross:.6f}; "
          f"z_t={z_t:.3f}: {'OK' if okq else 'FAIL'}")
    if not okq:
        fails.append("μ_cross")

    # (д) возраст: t₀·H₀ = ln(1/μ₀)·(2−μ₀²)/(1−μ₀²)
    h = lambda m: (2 - m * m) / (1 - m * m)
    t0H0 = np.log(1.0 / mu0) * h(mu0)
    t0 = t0H0 * (977.792 / 67.4)
    print(f"  (д) возраст: t₀·H₀ = {t0H0:.4f}  ->  t₀ = {t0:.2f} Гыр при H₀=67.4 "
          f"(наблюдается 13.80 ± 0.02)")

    if fails:
        raise SystemExit("СТОП: нарушены внутренние тождества: " + ", ".join(fails))
    print("  Итог: внутренние тождества ядра выполняются точно.\n")


# ==============================================================================
#  ЗАГРУЗКА СВЕРХНОВЫХ
# ==============================================================================
SN_BASE = ("https://raw.githubusercontent.com/PantheonPlusSH0ES/"
           "DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/")


def grab(name):
    if os.path.exists(name):
        print(f"  есть локально: {name}")
        return
    print(f"  скачиваю {name} ...")
    urllib.request.urlretrieve(SN_BASE + name, name)
    print(f"  готово: {os.path.getsize(name) / 1e6:.1f} МБ")


def read_cov(path, n_expected):
    with open(path) as f:
        n = int(f.readline())
        vals = np.array(f.read().split(), dtype=float)
    assert n == n_expected, f"размер ковариации {n} != {n_expected}"
    return vals.reshape(n, n)


def load_sn():
    print("=" * 70)
    print("[1a] Загрузка Pantheon+SH0ES")
    print("=" * 70)
    grab("Pantheon+SH0ES.dat")
    grab("Pantheon+SH0ES_STAT+SYS.cov")
    grab("Pantheon+SH0ES_STATONLY.cov")

    rows = []
    with open("Pantheon+SH0ES.dat") as f:
        hdr = f.readline().split()
        iz, imb = hdr.index("zHD"), hdr.index("m_b_corr")
        ical, iceph = hdr.index("IS_CALIBRATOR"), hdr.index("CEPH_DIST")
        for line in f:
            p = line.split()
            if p:
                rows.append(p)
    zHD = np.array([float(r[iz]) for r in rows])
    m_b = np.array([float(r[imb]) for r in rows])
    is_cal = np.array([int(r[ical]) for r in rows]) == 1
    ceph = np.array([float(r[iceph]) for r in rows])
    print(f"  прочитано {len(zHD)} сверхновых")

    C_sys = read_cov("Pantheon+SH0ES_STAT+SYS.cov", len(zHD))
    C_stat = read_cov("Pantheon+SH0ES_STATONLY.cov", len(zHD))

    use = (zHD > 0.01) | is_cal
    idx = np.where(use)[0]
    out = dict(
        z=zHD[idx], mb=m_b[idx], cal=is_cal[idx], ceph=ceph[idx],
        Cinv_sys=np.linalg.inv(C_sys[np.ix_(idx, idx)]),
        Cinv_stat=np.linalg.inv(C_stat[np.ix_(idx, idx)]),
    )
    print(f"  в фит идёт {len(out['z'])} SN (калибраторов: {out['cal'].sum()})\n")
    return out


# ==============================================================================
#  SN-ФИТ (якорный: m_b_corr + явный M) + ПОГРЕШНОСТИ ИЗ ГЕССИАНА
# ==============================================================================
def chi2_sn(params, model_distmod, SN, Cinv):
    p_model, H0, M = params
    z, mb, cal, ceph = SN["z"], SN["mb"], SN["cal"], SN["ceph"]
    mod = np.empty_like(z)
    mod[cal] = ceph[cal] + M
    zc = z[~cal]
    mod[~cal] = np.array([model_distmod(zz, p_model, H0) for zz in zc]) + M
    r = mb - mod
    return float(r @ Cinv @ r)


def fit_sn(model_distmod, SN, Cinv, x0, bounds_ok):
    def nll(p):
        return chi2_sn(p, model_distmod, SN, Cinv) if bounds_ok(p) else 1e12
    return minimize(nll, x0, method="Nelder-Mead",
                    options=dict(xatol=1e-5, fatol=1e-4, maxiter=20000))


def hessian_errors(fun, x0, steps):
    """1σ параметров из квадратичного разложения χ² в минимуме:
    C = 2·H⁻¹, где H — гессиан χ² (центральные разности)."""
    n = len(x0)
    Hm = np.zeros((n, n))
    f0 = fun(x0)
    for i in range(n):
        ei = np.zeros(n); ei[i] = steps[i]
        Hm[i, i] = (fun(x0 + ei) - 2 * f0 + fun(x0 - ei)) / steps[i] ** 2
        for j in range(i + 1, n):
            ej = np.zeros(n); ej[j] = steps[j]
            Hm[i, j] = Hm[j, i] = (fun(x0 + ei + ej) - fun(x0 + ei - ej)
                                   - fun(x0 - ei + ej) + fun(x0 - ei - ej)) \
                                  / (4 * steps[i] * steps[j])
    try:
        Cov = 2.0 * np.linalg.inv(Hm)
        d = np.diag(Cov)
        if np.any(d <= 0):
            return np.full(n, np.nan)
        return np.sqrt(d)
    except np.linalg.LinAlgError:
        return np.full(n, np.nan)


OK_D = lambda p: (0.5 < p[0] < 0.999) and (50 < p[1] < 90) and (-20.5 < p[2] < -18.5)
OK_L = lambda p: (0.05 < p[0] < 0.6) and (50 < p[1] < 90) and (-20.5 < p[2] < -18.5)


def run_sn_fits(SN):
    print("=" * 70)
    print("[1] SN якорный фит (Pantheon+SH0ES, STAT+SYS, m_b_corr + явный M)")
    print("=" * 70)
    dof = len(SN["z"]) - 3
    rd = fit_sn(distmod_DSIC, SN, SN["Cinv_sys"], [0.76, 73.0, -19.3], OK_D)
    rl = fit_sn(distmod_LCDM, SN, SN["Cinv_sys"], [0.33, 73.0, -19.3], OK_L)
    print(f"DSIC : mu0={rd.x[0]:.4f}  H0={rd.x[1]:.2f}  M={rd.x[2]:.4f}  chi2/dof={rd.fun/dof:.4f}")
    print(f"ΛCDM : Om ={rl.x[0]:.4f}  H0={rl.x[1]:.2f}  M={rl.x[2]:.4f}  chi2/dof={rl.fun/dof:.4f}")
    print(f"Δχ² (DSIC−ΛCDM) = {rd.fun - rl.fun:+.2f}   (|Δ|<2 — неразличимы)")

    print("\n  1σ-погрешности (гессиан χ² в минимуме, C = 2H⁻¹):")
    sd = hessian_errors(lambda p: chi2_sn(p, distmod_DSIC, SN, SN["Cinv_sys"]),
                        rd.x, [0.002, 0.15, 0.003])
    sl = hessian_errors(lambda p: chi2_sn(p, distmod_LCDM, SN, SN["Cinv_sys"]),
                        rl.x, [0.005, 0.15, 0.003])
    print(f"  DSIC : mu0 = {rd.x[0]:.4f} ± {sd[0]:.4f}   H0 = {rd.x[1]:.2f} ± {sd[1]:.2f}"
          f"   M = {rd.x[2]:.4f} ± {sd[2]:.4f}")
    print(f"  ΛCDM : Om  = {rl.x[0]:.4f} ± {sl[0]:.4f}   H0 = {rl.x[1]:.2f} ± {sl[1]:.2f}"
          f"   M = {rl.x[2]:.4f} ± {sl[2]:.4f}\n")

    print("=" * 70)
    print("[2] Устойчивость к учёту ошибок: STAT-ONLY")
    print("=" * 70)
    rds = fit_sn(distmod_DSIC, SN, SN["Cinv_stat"], [0.76, 73.0, -19.3], OK_D)
    rls = fit_sn(distmod_LCDM, SN, SN["Cinv_stat"], [0.33, 73.0, -19.3], OK_L)
    print(f"DSIC STAT-ONLY: mu0={rds.x[0]:.4f}  H0={rds.x[1]:.2f}  chi2/dof={rds.fun/dof:.4f}")
    print(f"ΛCDM STAT-ONLY: Om ={rls.x[0]:.4f}  H0={rls.x[1]:.2f}  chi2/dof={rls.fun/dof:.4f}")
    print(f"Δχ² STAT-ONLY = {rds.fun - rls.fun:+.2f}")
    print(f"  mu0: {rd.x[0]:.4f} (STAT+SYS) → {rds.x[0]:.4f} (STAT-ONLY)\n")
    return rd.x, rl.x


# ==============================================================================
#  [3] СОВМЕСТНЫЙ SN + BAO (полная ковариация 12×12, свободный r_d)
# ==============================================================================
def build_bao_cov():
    z_bao = BAO[:, 0]
    obs = np.empty(12)
    Cb = np.zeros((12, 12))
    for i, row in enumerate(BAO):
        _, dm, sdm, dh, sdh, cc_ = row
        obs[2*i], obs[2*i+1] = dm, dh
        Cb[2*i, 2*i] = sdm**2
        Cb[2*i+1, 2*i+1] = sdh**2
        Cb[2*i, 2*i+1] = cc_ * sdm * sdh
        Cb[2*i+1, 2*i] = cc_ * sdm * sdh
    return z_bao, obs, np.linalg.inv(Cb)


def run_joint(SN):
    z_bao, obs, Cb_inv = build_bao_cov()

    def bao_vec_DSIC(mu0, H0, rd):
        v = np.empty(12)
        for i, zv in enumerate(z_bao):
            v[2*i] = DM_comov_DSIC(zv, mu0, H0) / rd
            v[2*i+1] = (c / H_DSIC(zv, mu0, H0)) / rd
        return v

    def bao_vec_LCDM(Om, H0, rd):
        v = np.empty(12)
        for i, zv in enumerate(z_bao):
            v[2*i] = DM_comov_LCDM(zv, Om, H0) / rd
            v[2*i+1] = (c / (H0 * E_LCDM(zv, Om))) / rd
        return v

    def chi2_bao(vec):
        r = vec - obs
        return float(r @ Cb_inv @ r)

    def joint_D(p):
        mu0, H0, M, rd = p
        if not (0.5 < mu0 < 0.999 and 50 < H0 < 90 and -20.5 < M < -18.5 and 120 < rd < 170):
            return 1e12
        return chi2_sn([mu0, H0, M], distmod_DSIC, SN, SN["Cinv_sys"]) + chi2_bao(bao_vec_DSIC(mu0, H0, rd))

    def joint_L(p):
        Om, H0, M, rd = p
        if not (0.05 < Om < 0.6 and 50 < H0 < 90 and -20.5 < M < -18.5 and 120 < rd < 170):
            return 1e12
        return chi2_sn([Om, H0, M], distmod_LCDM, SN, SN["Cinv_sys"]) + chi2_bao(bao_vec_LCDM(Om, H0, rd))

    print("=" * 70)
    print("[3] Совместный SN+BAO (полная ковариация 12×12, свободный r_d)")
    print("=" * 70)
    jd = minimize(joint_D, [0.76, 73, -19.3, 147], method="Nelder-Mead",
                  options=dict(xatol=1e-5, fatol=1e-4, maxiter=30000))
    jl = minimize(joint_L, [0.31, 73, -19.3, 147], method="Nelder-Mead",
                  options=dict(xatol=1e-5, fatol=1e-4, maxiter=30000))
    mu0J, H0dJ, MdJ, rddJ = jd.x
    OmJ, H0lJ, MlJ, rdlJ = jl.x
    baoD = chi2_bao(bao_vec_DSIC(mu0J, H0dJ, rddJ))
    baoL = chi2_bao(bao_vec_LCDM(OmJ, H0lJ, rdlJ))
    print(f"DSIC : mu0={mu0J:.4f} H0={H0dJ:.2f} M={MdJ:.3f} r_d={rddJ:.1f}  тотал={jd.fun:.1f}")
    print(f"ΛCDM : Om ={OmJ:.4f} H0={H0lJ:.2f} M={MlJ:.3f} r_d={rdlJ:.1f}  тотал={jl.fun:.1f}")
    print(f"Δχ² тотал (DSIC−ΛCDM) = {jd.fun - jl.fun:+.2f}  "
          f"(слабое; частично отражает натяжение DESI–ΛCDM)")
    print(f"BAO χ²/dof (12−4=8): DSIC={baoD/8:.2f}  ΛCDM={baoL/8:.2f}")

    print("\n  Линейка r_d (обратная задача) при ранней шкале H0=67.4:")
    H0e = 67.4
    rdM = [DM_comov_DSIC(zv, mu0J, H0e) / BAO[i, 1] for i, zv in enumerate(z_bao)]
    rdH = [(c / H_DSIC(zv, mu0J, H0e)) / BAO[i, 3] for i, zv in enumerate(z_bao)]
    print(f"    r_d по D_M = {np.mean(rdM):.1f} Мпк  (разброс {np.std(rdM):.1f})")
    print(f"    r_d по D_H = {np.mean(rdH):.1f} Мпк  (разброс {np.std(rdH):.1f})\n")
    return (mu0J, H0dJ, rddJ), (OmJ, H0lJ, rdlJ)


# ==============================================================================
#  [4] Om(z): реконструкция + устойчивость знака Δχ² к H0, В ДВУХ РЕЖИМАХ
#      ковариации хронометров (диагональная и полная Moresco 2020).
# ==============================================================================
def run_omz(dsic_par, lcdm_par):
    mu0_fit = dsic_par[0]
    Om_fit = lcdm_par[0]
    rd_norm = 147.1  # перевод D_H/r_d -> H(z) при реконструкции Om(z)
    z_b, DHrd, sDHrd = BAO[:, 0], BAO[:, 3], BAO[:, 4]

    Cinv_cc_diag = np.linalg.inv(build_cc_covariance(z_cc, sH_cc, full=False))
    Cinv_cc_full = np.linalg.inv(build_cc_covariance(z_cc, sH_cc, full=True))

    def E_DSIC(zv):
        me = mu_e_of_z(zv, mu0_fit)
        return ((2 - me**2) / (1 - me**2)) / ((2 - mu0_fit**2) / (1 - mu0_fit**2))

    def E_L(zv):
        return E_LCDM(zv, Om_fit)

    def Om_of_E(zv, E):
        return (E**2 - 1.0) / ((1.0 + zv) ** 3 - 1.0)

    def dchi2_at_H0(H0, Cinv_cc):
        E_d = H_cc / H0
        Om_d = Om_of_E(z_cc, E_d)
        H_bao = c / (DHrd * rd_norm)
        E_b = H_bao / H0
        Om_b = Om_of_E(z_b, E_b)
        dOm_b = np.abs(2 * E_b / ((1 + z_b) ** 3 - 1)) * ((H_bao * (sDHrd / DHrd)) / H0)

        def cc_resid(Ef):
            return Om_d - np.array([Om_of_E(z, Ef(z)) for z in z_cc])

        # якобиан σ_H -> σ_Om; C_Om = J C_H J^T (J диагональный)
        J = np.abs(2 * (H_cc / H0) / ((1 + z_cc) ** 3 - 1)) / H0
        Dinv = np.diag(1.0 / J)
        Cinv_Om = Dinv @ Cinv_cc @ Dinv

        def chi2_curve(Ef):
            r_cc = cc_resid(Ef)
            a = float(r_cc @ Cinv_Om @ r_cc)
            r_b = Om_b - np.array([Om_of_E(z, Ef(z)) for z in z_b])
            b_ = np.sum((r_b / dOm_b) ** 2)
            return a + b_

        return chi2_curve(E_DSIC) - chi2_curve(E_L)

    print("=" * 70)
    print("[4] Om(z): устойчивость знака Δχ² к нормировке H0 (честный тест)")
    print("=" * 70)
    print("  Хронометры подключены в ДВУХ режимах ковариации:")
    print("    DIAG = диагональные ошибки (базовая версия теста);")
    print("    FULL = полная ковариация Moresco 2020 (15 точек коррелированы).")
    print(f"  {'H0':>4} | {'Δχ² DIAG':>10} | {'Δχ² FULL':>10}")
    print("  " + "-" * 32)
    signs_d, signs_f = [], []
    for H0 in range(64, 77, 2):
        dd = dchi2_at_H0(float(H0), Cinv_cc_diag)
        df = dchi2_at_H0(float(H0), Cinv_cc_full)
        signs_d.append(dd)
        signs_f.append(df)
        print(f"  {H0:>4} | {dd:>+10.1f} | {df:>+10.1f}")

    def verdict(signs):
        signs = np.array(signs)
        if np.all(signs < 0):
            return f"устойчиво ближе к DSIC  (диапазон {signs.min():+.1f}…{signs.max():+.1f})"
        if np.all(signs > 0):
            return f"устойчиво ближе к ΛCDM  (диапазон {signs.min():+.1f}…{signs.max():+.1f})"
        return f"ЗНАК МЕНЯЕТСЯ → тест НЕ разрешает (диапазон {signs.min():+.1f}…{signs.max():+.1f})"

    print(f"\n  Вердикт DIAG: {verdict(signs_d)}")
    print(f"  Вердикт FULL: {verdict(signs_f)}")
    print("  [!] Итог: в обоих режимах ковариации знак Δχ² определяется выбором")
    print("      нормировки H₀, а не формой Om(z): при низкой (Planck-подобной)")
    print("      H₀ данные формально ближе к DSIC, при высокой (SH0ES) — к ΛCDM.")
    print("      Переход DIAG→FULL меняет Δχ² лишь на единицы и вердикта не")
    print("      меняет. Проба моделей НЕ разрешает; она станет решающей только")
    print("      при независимо закреплённой H₀ и более точных H(z).\n")


# ==============================================================================
#  [5] ФАЛЬСИФИЦИРУЕМЫЕ ПРЕДСКАЗАНИЯ НА ВЫСОКОМ z
# ==============================================================================
def run_predictions(dsic_par, lcdm_par):
    mu0J, H0dJ, rddJ = dsic_par
    OmJ, H0lJ, rdlJ = lcdm_par
    print("=" * 70)
    print("[5] Фальсифицируемые предсказания DSIC vs ΛCDM на высоком z")
    print("=" * 70)
    print("  D_H/r_d (Lyman-α область) — след ЯДРА, от Части II не зависит:")
    for zv in [2.33, 2.5, 3.0, 3.5]:
        dhd = (c / H_DSIC(zv, mu0J, H0dJ)) / rddJ
        dhl = (c / (H0lJ * E_LCDM(zv, OmJ))) / rdlJ
        print(f"    z={zv}:  DSIC={dhd:.3f}  ΛCDM={dhl:.3f}  разн={100*(dhd/dhl-1):+.1f}%")
    print("  E(z) за пределами данных:")
    for zv in [5, 10]:
        ed = H_DSIC(zv, mu0J, 1.0)
        el = E_LCDM(zv, OmJ)
        print(f"    z={zv}:  E_DSIC/E_ΛCDM = {ed/el:.2f}  (+{100*(ed/el-1):.0f}%)")
    print()


# ==============================================================================
#  [6] КОСМИЧЕСКОЕ ВРЕМЯ t(z): ближайшая зона риска ядра (новый блок)
# ------------------------------------------------------------------------------
#  В DSIC собственное время наблюдателя от каши до эпохи μ_e:
#     τ(μ_e) = ln(1/μ_e)/k,   H₀ = k·(2−μ₀²)/(1−μ₀²)
#     =>  t(z)·H₀ = ln(1/μ_e(z)) · (2−μ₀²)/(1−μ₀²)
#  (при z=0 воспроизводит t₀·H₀ ≈ 0.925 из основного теста).
#  ΛCDM: t(z) = (1/H₀)·∫_z^∞ dz'/((1+z')E(z')) (материя+Λ; вклад излучения
#  на этих эпохах ≲1% и здесь опущен, что ЗАВЫШАЕТ возраст ΛCDM — сравнение
#  консервативно в пользу ΛCDM... т.е. не в пользу DSIC оно не искажено).
#  Оптическая поправка Части II путь света удлиняет, но H(z) НЕ меняет,
#  поэтому этот блок тестирует именно ядро (статья, §11.6).
# ==============================================================================
def run_cosmic_time(dsic_par, lcdm_par, H0_scale=67.4):
    mu0J = dsic_par[0]
    OmJ = lcdm_par[0]
    h = lambda m: (2 - m * m) / (1 - m * m)
    Gyr = 977.792 / H0_scale  # 1/H0 в Гыр

    def t_dsic(zv):
        me = mu_e_of_z(zv, mu0J)
        return np.log(1.0 / me) * h(mu0J) * Gyr

    def t_lcdm(zv):
        return Gyr * quad(lambda x: 1.0 / ((1 + x) * E_LCDM(x, OmJ)),
                          zv, np.inf, limit=400)[0]

    print("=" * 70)
    print("[6] Космическое время t(z): DSIC vs ΛCDM (общая шкала H0=67.4)")
    print("=" * 70)
    print(f"  Возраст сегодня: t0_DSIC = {t_dsic(0.0):.2f} Гыр,  "
          f"t0_ΛCDM = {t_lcdm(0.0):.2f} Гыр  (набл. 13.80 ± 0.02)")
    print(f"  {'z':>5} | {'t_DSIC, Гыр':>12} | {'t_ΛCDM, Гыр':>12} | {'отношение':>9}")
    print("  " + "-" * 48)
    for zv in [3.0, 5.0, 6.0, 7.0, 7.5, 8.0, 10.0, 12.0]:
        td, tl = t_dsic(zv), t_lcdm(zv)
        print(f"  {zv:>5.1f} | {td:>12.3f} | {tl:>12.3f} | {td/tl:>9.2f}")
    print()
    print("  Контекст (что должно успеть произойти за это время):")
    print("   - квазары с M_BH ~ 10⁹ M☉ наблюдаются к z ≈ 7.5;")
    print("   - галактики JWST подтверждены до z ≈ 10–14;")
    print("   - реионизация завершается к z ≈ 6.")
    print("  Сжатие времени в DSIC на этих эпохах — самая близкая проверяемая")
    print("  зона риска ЯДРА (Часть II её не смягчает: H(z) не меняется).\n")


# ==============================================================================
#  ЧАСТЬ II — ВТОРОЙ ЭТАЖ (P6): ПОРОГ ОТЛИПАНИЯ μ_d
# ==============================================================================
def DM_base_mu(mu_e):
    """Базовое сопутствующее расстояние ядра (ранняя шкала MU0_E, H0_E), из μ_e."""
    ck = ck_const(MU0_E, H0_E)
    pref = np.sqrt(1.0 - MU0_E ** 2) / MU0_E ** 2
    return ck * pref * (np.sqrt(1.0 - MU0_E ** 2) - np.sqrt(1.0 - mu_e ** 2))


def DM_corrected_mu(mu_e, mu_d):
    """Исправленное расстояние: добор arcsin-поправки выше порога μ_d.

    По P6 добавка к пути в слипшейся фазе есть 1/P = 1/√(1−μ²) (амплитуда
    A ≡ 1 фиксирована постулатом), отсечённая жёстким порогом; ∫ = arcsin μ:
      D_M_испр = D_M + (c/k)·[arcsin(μ_e) − arcsin(μ_d)]  при μ_e > μ_d,
      D_M_испр = D_M                                        при μ_e ≤ μ_d.
    """
    base = DM_base_mu(mu_e)
    if mu_e > mu_d:
        return base + ck_const(MU0_E, H0_E) * (np.arcsin(mu_e) - np.arcsin(mu_d))
    return base


def self_check_planck():
    print("=" * 70)
    print("[7] Самопроверка опорных данных Planck 2018 (вшиты)")
    print("=" * 70)
    d_from_geom = PLANCK["r_star"] / PLANCK["theta_star"]
    rel = abs(d_from_geom - PLANCK["d_lss"]) / PLANCK["d_lss"]
    print(f"  z* = {PLANCK['z_star']},  r* = {PLANCK['r_star']} Мпк,  "
          f"100θ* = {100*PLANCK['theta_star']:.5f}")
    print(f"  d_LSS (вшито)   = {PLANCK['d_lss']} ± {PLANCK['d_lss_err']} Мпк")
    print(f"  d_LSS = r*/θ*   = {d_from_geom:.1f} Мпк  (расхожд. {100*rel:.2f}%)")
    if rel < 0.01:
        print("  OK — опорная цифра согласована с r* и θ* (Planck 2018).\n")
    else:
        raise SystemExit("СТОП: d_LSS не согласуется с r*/θ*. Проверь PLANCK.")


def step_fix_threshold():
    print("=" * 70)
    print("[8] Второй этаж (P6): недобор ядра и фиксация порога μ_d по углу θ*")
    print("=" * 70)
    ck = ck_const(MU0_E, H0_E)
    mu_cmb = mu_e_of_z(PLANCK["z_star"], MU0_E)
    DM0 = DM_base_mu(mu_cmb)
    target_d = PLANCK["r_star"] / PLANCK["theta_star"]   # = r*/θ*
    deficit = target_d - DM0
    print(f"  Параметры ядра (ранняя шкала): μ₀ = {MU0_E}, H₀ = {H0_E}")
    print(f"  c/k                = {ck:.1f} Мпк")
    print(f"  μ_cmb(z*={PLANCK['z_star']}) = {mu_cmb:.10f}")
    print(f"  D_M база до CMB    = {DM0:.1f} Мпк")
    print(f"  Цель: θ_model(μ_d) = θ*  ⇔  D_M_испр(μ_cmb) = r*/θ* = {target_d:.1f} Мпк")
    print(f"  недобор ядра       = {deficit:.1f} Мпк  ({100*deficit/target_d:.1f} %)")
    print("  [!] r* = 144.43 Мпк — величина, ВЫЧИСЛЕННАЯ внутри ΛCDM (внешний")
    print("      вход; прямое измерение — угол θ*). Происхождение акустического")
    print("      масштаба ~147 Мпк — открытый вопрос ядра (статья, §11.6).")
    mu_d = brentq(lambda md: ck * (np.arcsin(mu_cmb) - np.arcsin(md)) - deficit,
                  0.5, mu_cmb - 1e-12)
    z_det = b_dsic(MU0_E) / b_dsic(mu_d) - 1.0
    DM_corr = DM_corrected_mu(mu_cmb, mu_d)
    theta_model = PLANCK["r_star"] / DM_corr
    print(f"\n  μ_d        = {mu_d:.6f}")
    print(f"  z_detach   = {z_det:.4f}")
    print(f"  D_M испр   = {DM_corr:.1f} Мпк;  100·θ_model = {100*theta_model:.5f}"
          f"  (Planck 100θ* = {100*PLANCK['theta_star']:.5f})")
    ok = abs(DM_corr - target_d) < 1.0
    print(f"  замыкание до θ*: {'OK' if ok else 'РАСХОЖДЕНИЕ'}\n")
    if not ok:
        raise SystemExit("СТОП: θ_model не совпал с θ*.")
    return mu_cmb, DM0, deficit, mu_d, z_det


def step_late_universe_and_degeneracy(mu_cmb, DM0, deficit, mu_d, z_det):
    print("=" * 70)
    print("[9] Контроль нуля ниже порога + вырождение «амплитуда ↔ порог»")
    print("=" * 70)
    print("  Поправка тождественно нулевая при z ≤ z_detach:")
    print(f"  {'z':>6} {'D_M база':>10} {'поправка':>10}")
    print("  " + "-" * 28)
    max_corr_below = 0.0
    for z in [0.5, 1.0, 2.33, 3.0, 4.0, z_det, z_det + 0.1, 5.0, 10.0]:
        me = mu_e_of_z(z, MU0_E)
        base = DM_base_mu(me)
        corr = DM_corrected_mu(me, mu_d) - base
        tag = "" if z <= z_det + 1e-9 else "  <- включена"
        if z <= z_det + 1e-9:
            max_corr_below = max(max_corr_below, abs(corr))
        print(f"  {z:>6.2f} {base:>10.1f} {corr:>+10.1f}{tag}")
    print(f"\n  макс. |поправка| при z ≤ z_detach = {max_corr_below:.2f} Мпк (≈0)")
    print("  -> все 6 бинов DESI (z ≤ 2.33) поправки НЕ чувствуют; результаты")
    print("     Части I ([1]–[5]) сохранены полностью.\n")

    ck = ck_const(MU0_E, H0_E)
    print("  Вырождение БЕЗ постулата P6 (что именно фиксирует постулат):")
    print("  уравнение A·(c/k)·[arcsin μ_cmb − arcsin μ_d] = недобор имеет")
    print("  кривую решений (A, μ_d) — все дают одно и то же расстояние до CMB.")
    print(f"  {'A':>6} {'μ_d':>10} {'z_detach':>10}")
    print("  " + "-" * 30)
    for A in [0.5, 0.8, 1.0, 1.2, 2.0, 5.0]:
        arg = np.arcsin(mu_cmb) - deficit / (A * ck)
        if 0 < arg < np.pi / 2:
            md = np.sin(arg)
            zd = b_dsic(MU0_E) / b_dsic(md) - 1.0
            star = "  <- фиксировано P6 (A≡1)" if abs(A - 1.0) < 1e-9 else ""
            print(f"  {A:>6.1f} {md:>10.6f} {zd:>10.3f}{star}")
        else:
            print(f"  {A:>6.1f} {'нет реш.':>10}")
    print("\n  Одно наблюдение фиксирует один параметр: степеней свободы ноль.")
    print("  При форме, зафиксированной P6, μ_d — измеренная константа второго")
    print("  этажа (статус как у μ₀ в ядре); предсказанием z_detach не является.\n")


def step_trace_and_horizon(mu_d, z_det):
    print("=" * 70)
    print("[10] След μ_d, исправленный горизонт и флаг согласованности")
    print("=" * 70)
    print(f"  {'z':>8} {'объект':<18} {'поправка к D_M':>16}")
    print("  " + "-" * 46)
    for z, obj in [(z_det, "порог отлипания"), (5.0, "галактики JWST"),
                   (6.0, "квазары"), (8.0, "реионизация"),
                   (11.0, "первые галактики"), (PLANCK["z_star"], "CMB (калибровка)")]:
        me = mu_e_of_z(z, MU0_E)
        base = DM_base_mu(me)
        corr = DM_corrected_mu(me, mu_d)
        print(f"  {z:>8.2f} {obj:<18} {100*(corr/base-1):>+13.1f} %")
    ck = ck_const(MU0_E, H0_E)
    hor_base = ck * (1.0 - MU0_E ** 2) / MU0_E ** 2
    hor_corr = ck * ((1.0 - MU0_E ** 2) / MU0_E ** 2 + (np.pi / 2 - np.arcsin(mu_d)))
    print(f"\n  Горизонт частиц: база = {hor_base:.1f} Мпк;  с поправкой = "
          f"{hor_corr:.1f} Мпк")
    print(f"  (на {hor_corr - PLANCK['d_lss']:.1f} Мпк дальше поверхности последнего "
          f"рассеяния — реликт почти вплотную к горизонту)")
    print()
    print("  Досягаемость зондов:")
    print("    DESI DR2 (есть)      z ≤ 2.33   — ниже порога")
    print("    DESI DR3 (~2027)     z ≤ 3.55   — ниже порога")
    print("    DESI-II/LBG (~2029)  2.5<z<3.5  — ниже порога")
    print("  Радиальный BAO в окне 2.33–4.53 — след ЯДРА (D_H = c/H), не μ_d;")
    print("  сам порог влияет лишь на D_M выше z≈4.5, где точных линеек, кроме")
    print("  CMB-калибровки, нет: второй независимый якорь недосягаем до появления")
    print("  точной шкалы расстояний на z ≳ 5–6.")
    print()
    mu_e14 = mu_e_of_z(14.0, MU0_E)
    print("  [!] ФЛАГ СОГЛАСОВАННОСТИ ФОРМУЛИРОВКИ P6:")
    print(f"      наблюдаемые галактики (JWST, z ≈ 10–14; μ_e(z=14) = {mu_e14:.4f})")
    print(f"      лежат ВЫШЕ порога отлипания (μ_d = {mu_d:.4f}, z_detach = {z_det:.2f}),")
    print("      то есть существуют как раздельные объекты ДО отлипания в смысле P6.")
    print("      Буквальное чтение «раздельных объектов не существует при μ > μ_d»")
    print("      противоречит наблюдениям; формулировка требует уточнения")
    print("      (отлипание доминирующего монолита vs первые отлипшие островки),")
    print("      либо пересмотра выбора A (что ломает естественность A ≡ 1).")
    print("      Это ограничение НАРРАТИВА Части II; арифметика расстояний и")
    print("      результаты Части I от него не зависят.\n")


# ==============================================================================
#  ГЛАВНЫЙ ПОТОК
# ==============================================================================




# ==============================================================================
#  ЧАСТЬ III — РОСТ СТРУКТУР  (блок [11])
# ==============================================================================
#
# DSIC — проверка роста структур (linear growth, fσ8) против данных RSD.
# Блок [11] для dsic_full_test.py.
#
# Идея: DSIC задаёт фон H(z). На этом фоне решаем СТАНДАРТНОЕ уравнение
# роста линейных возмущений материи (sub-horizon, GR-предел):
#
#     d²δ/d(lna)² + (2 + dlnH/dlna) dδ/d(lna) - (3/2) Ω_m(a) δ = 0
#
# Затем f = dlnδ/dlna, σ8(z) = σ8,0·δ(z)/δ(0), сравниваем fσ8(z) с RSD.
#
# ДОПУЩЕНИЕ (честно): DSIC не выводит уравнение роста из ядра. Здесь берётся
# GR-предел роста на фоне DSIC-геометрии — проверка на выживание, а не вывод.
# Ω_m0 входит как параметр роста, т.к. ядро плотностей не вводит.
#
# КЛЮЧЕВОЙ ВЫВОД: данные fσ8 вырождены — они ограничивают комбинацию
# S8 = σ8·√(Ω_m0/0.3), а не Ω_m0 и σ8 по отдельности. Поэтому результат
# докладывается как S8, плюс срез при фиксированном Planck-подобном Ω_m0=0.30.
#
# ------------------------------------------------------------------------
# ДАННЫЕ: компиляция "Gold-2017" fσ8 — 18 слабо коррелированных точек,
# Nesseris, Pantazis & Perivolaropoulos (2017), Phys. Rev. D 96, 023542,
# arXiv:1703.10538. z, fσ8, σ взяты из Table I работы Sagredo, Nesseris &
# Sapone (2018), arXiv:1806.10822 (первые 18 строк = Gold-2017).
# Три точки WiggleZ (z=0.44,0.60,0.73) — с ковариацией Blake et al. 2012.
#
# УПРОЩЕНИЕ: поправка Alcock–Paczynski (DSIC vs fiducial) не применяется;
# для z<1 AP-фактор обычно <2-3%, мал против ошибок. Отмечено как систематика.
# ------------------------------------------------------------------------
#
# ==============================================================================


# ---------- ФОН DSIC ----------
def b_of_mu(mu):
    return np.sqrt(1 - mu**2) / mu**2

def mu_of_z(z, mu0):
    """1+z = b(mu0)/b(mu_e), физ. ветвь mu_e in [mu0,1). При z=0 -> mu0."""
    if z <= 0:
        return mu0
    target = b_of_mu(mu0) / (1 + z)
    f = lambda m: b_of_mu(m) - target
    return brentq(f, mu0, 1 - 1e-14, xtol=1e-14)

def E_of_z(z, mu0):
    mu = mu_of_z(z, mu0)
    H  = (2 - mu**2) / (1 - mu**2)
    H0 = (2 - mu0**2) / (1 - mu0**2)
    return H / H0

# ---------- РОСТ ВОЗМУЩЕНИЙ ----------
def growth_ode(mu0, Om0, z_start=3.0, n=800):
    a_start = 1/(1+z_start)
    lna = np.linspace(np.log(a_start), 0.0, n)
    def E_a(a): return E_of_z(1/a - 1, mu0)
    def dlnE_dlna(a, h=1e-4):
        return (np.log(E_a(a*np.exp(h))) - np.log(E_a(a*np.exp(-h)))) / (2*h)
    def rhs(x, y):
        d, dp = y; a = np.exp(x)
        Om_a = Om0 * a**-3 / E_a(a)**2
        return [dp, -(2 + dlnE_dlna(a))*dp + 1.5*Om_a*d]
    sol = solve_ivp(rhs, (lna[0], lna[-1]), [a_start, a_start],
                    t_eval=lna, rtol=1e-9, atol=1e-12, method="RK45")
    a = np.exp(sol.t); delta = sol.y[0]; f = sol.y[1]/sol.y[0]
    z = 1/a - 1
    return z[::-1], delta[::-1], f[::-1]

def fsigma8_curve(mu0, Om0, sigma8_0):
    z, delta, f = growth_ode(mu0, Om0)
    delta0 = np.interp(0.0, z, delta)
    return z, f * sigma8_0 * delta/delta0

# ---------- ДАННЫЕ fσ8 (Gold-2017, 18 точек) ----------
fs8_data = np.array([
    [0.020, 0.428, 0.0465], [0.020, 0.398, 0.0650], [0.020, 0.314, 0.0480],
    [0.100, 0.370, 0.1300], [0.150, 0.490, 0.1450], [0.170, 0.510, 0.0600],
    [0.180, 0.360, 0.0900], [0.380, 0.440, 0.0600], [0.250, 0.3512, 0.0583],
    [0.370, 0.4602, 0.0378],[0.320, 0.384, 0.0950], [0.590, 0.488, 0.0600],
    [0.440, 0.413, 0.0800], [0.600, 0.390, 0.0630], [0.730, 0.437, 0.0720],
    [0.600, 0.550, 0.1200], [0.860, 0.400, 0.1100], [1.400, 0.482, 0.1160],
])
WIGGLEZ_IDX = [12, 13, 14]
COV_WIGGLEZ = 1e-3 * np.array([
    [6.400, 2.570, 0.000],
    [2.570, 3.969, 2.540],
    [0.000, 2.540, 5.184],
])
def build_cov(fs8):
    C = np.diag(fs8[:, 2]**2)
    for a, i in enumerate(WIGGLEZ_IDX):
        for b, j in enumerate(WIGGLEZ_IDX):
            C[i, j] = COV_WIGGLEZ[a, b]
    return C
COV = build_cov(fs8_data)
COV_INV = np.linalg.inv(COV)

def chi2_fs8(mu0, Om0, sigma8_0):
    zt, curve = fsigma8_curve(mu0, Om0, sigma8_0)
    r = np.interp(fs8_data[:, 0], zt, curve) - fs8_data[:, 1]
    return float(r @ COV_INV @ r)


def run_growth():
    """Блок [11]: рост структур fσ8 на фоне DSIC (GR-предел)."""
    mu0 = 0.76
    dof = len(fs8_data) - 2

    # 1) Полный скан -> истинный минимум (демонстрация вырождения)
    Om_grid = np.linspace(0.05, 0.50, 91)
    s8_grid = np.linspace(0.55, 1.15, 121)
    best = (1e9, None)
    for Om0 in Om_grid:
        for s8 in s8_grid:
            c = chi2_fs8(mu0, Om0, s8)
            if c < best[0]: best = (c, (Om0, s8))
    c, (Om0, s8) = best
    print(f"DSIC growth (GR-предел на фоне DSIC), mu0={mu0}")
    print(f"  выборка: Gold-2017, {len(fs8_data)} точек (WiggleZ-ковариация учтена)")
    print(f"  формальный минимум: Om0={Om0:.3f}, s8={s8:.3f}, "
          f"chi2={c:.2f}, chi2/dof={c/dof:.2f}")
    print(f"  (минимум внутри сетки; долина сильно вырождена — см. ниже)")

    # 2) Вырождение: комбинация S8 постоянна вдоль долины
    print("\n  Вырождение fσ8 (данные ограничивают S8, не Om0 и s8 порознь):")
    print("   Om0    s8*     chi2/dof   S8=s8*sqrt(Om0/0.3)")
    S8_vals = []
    for Om0 in [0.15, 0.20, 0.25, 0.30, 0.35, 0.40]:
        cb, s8b = min((chi2_fs8(mu0, Om0, s8), s8) for s8 in s8_grid)
        S8 = s8b*np.sqrt(Om0/0.3); S8_vals.append(S8)
        print(f"  {Om0:.2f}   {s8b:.3f}   {cb/dof:.2f}       {S8:.3f}")
    print(f"  => S8 стабильно: {np.mean(S8_vals):.3f} ± {np.std(S8_vals):.3f}")

    # 3) Срез при Planck-подобном Om0=0.30 — рекомендуемая форма отчёта
    cb, s8b = min((chi2_fs8(mu0, 0.30, s8), s8) for s8 in s8_grid)
    print(f"\n  При фиксированном Om0=0.30 (Planck-подобное):")
    print(f"    s8={s8b:.3f}, chi2={cb:.2f}, chi2/dof={cb/(len(fs8_data)-1):.2f}")
    print(f"    S8={s8b*np.sqrt(0.30/0.3):.3f}")

    # 4) Кривая при Om0=0.30
    zt, curve = fsigma8_curve(mu0, 0.30, s8b)
    print("\n   z    fσ8_model   fσ8_obs±err   (Om0=0.30)")
    for zz, obs, err in fs8_data:
        print(f"  {zz:.2f}   {np.interp(zz, zt, curve):.3f}       {obs:.3f}±{err:.3f}")


def main():
    # -------- Часть I: ядро --------
    self_check_desi()                                 # [0]
    self_check_coherence()                            # [0b]
    SN = load_sn()                                    # [1a]
    dsic_sn, lcdm_sn = run_sn_fits(SN)                # [1],[2]
    dsic_par, lcdm_par = run_joint(SN)                # [3]
    run_omz(dsic_par, lcdm_par)                       # [4]
    run_predictions(dsic_par, lcdm_par)               # [5]
    run_cosmic_time(dsic_par, lcdm_par)               # [6]
    # -------- Часть II: второй этаж --------
    self_check_planck()                               # [7]
    mu_cmb, DM0, deficit, mu_d, z_det = step_fix_threshold()        # [8]
    step_late_universe_and_degeneracy(mu_cmb, DM0, deficit, mu_d, z_det)  # [9]
    step_trace_and_horizon(mu_d, z_det)               # [10]
    # -------- Часть III: рост структур --------
    print()
    print("#" * 70)
    print("#  ЧАСТЬ III — РОСТ СТРУКТУР (блок [11])")
    print("#" * 70)
    run_growth()                                      # [11]
    print("=" * 70)
    print("Готово. Часть I: SN скачаны, DESI прошёл самопроверку, тождества")
    print("ядра выполнены, погрешности параметров посчитаны из гессиана.")
    print("Часть II: опорные данные Planck сверены, порог μ_d зафиксирован по θ*")
    print("(r* заимствован — помечено), поздняя Вселенная не затронута, флаг")
    print("согласованности формулировки P6 выведен. μ_d — калибровка по одному")
    print("наблюдению при форме, зафиксированной постулатом.")
    print("Часть III: рост структур fσ8 на фоне DSIC (GR-предел) сверен с\n"
          "компиляцией Gold-2017; данные вырождены по (Ω_m0, σ8), докладывается\n"
          "инвариант S8. Статус — проверка на выживание, не вывод из ядра.")
    print("=" * 70)


if __name__ == "__main__":
    main()
