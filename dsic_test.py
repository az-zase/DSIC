#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dsic_test.py — независимая проверка всех чисел статьи
"Dual Scale Inversion Cosmology" (А. Засеев).

Запуск в Google Colab:
    !python dsic_test.py                # всё, кроме полного спектра реликта (около минуты)
    !python dsic_test.py --cmb          # плюс CAMB + plik-lite: четыре подбора таблицы 7 и проба f (около двух часов)
    !python dsic_test.py --cmb-f        # только проба f с полным спектром реликта (от 20 минут до часа)

Нужны numpy и scipy. Для --cmb скрипт сам ставит camb и скачивает planck-lite-py.
Все наблюдательные данные скачиваются из официальных репозиториев и сверяются
с копиями, вшитыми в этот файл:
  * сверхновые Pantheon+SH0ES — github.com/PantheonPlusSH0ES/DataRelease
  * BAO DESI DR2 (13 точек, полная ковариация) — github.com/CobayaSampler/bao_data
  * хронометры H(z), fσ8, чёрные дыры — таблицы ниже, источники в комментариях.

Соглашения (одни и те же для DSIC и ΛCDM):
  * плоская геометрия, безмассовые нейтрино, N_eff = 3.046, T_CMB = 2.7255 K;
  * фотометрическое расстояние d_L = (1 + z_hel)·D_M(z_HD) — стандарт Pantheon+;
  * звуковой горизонт: r_s(z) = ∫ c_s dz / H, z* = 1089.92, z_d = 1059.94 (Planck 2018);
    поправка приближения откалибрована на ΛCDM Planck и одинакова для обеих моделей.
"""
import os, sys, json, time, urllib.request
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import minimize, minimize_scalar, brentq

T_START = time.time()
RUN_CMB = "--cmb" in sys.argv
RUN_CMB_F = "--cmb-f" in sys.argv                # только проба f с полным спектром реликта
DATA_DIR = os.environ.get("DSIC_DATA", "dsic_data")
os.makedirs(DATA_DIR, exist_ok=True)
R = {}                                   # все итоговые числа — сюда, в конце печатаются как JSON


def head(title):
    print("\n" + "=" * 78 + "\n" + title + "\n" + "=" * 78)


# ==============================================================================
#  КОНСТАНТЫ
# ==============================================================================
C_KMS = 299792.458                       # км/с
GYR = 977.792                            # 1/(км/с/Мпк) в млрд лет
W_GAMMA = 2.4728e-5                      # ω_γ при T = 2.7255 K
N_EFF = 3.046
W_R = W_GAMMA * (1 + 0.2271 * N_EFF)     # фотоны + безмассовые нейтрино: 4.183e-5
Z_STAR, Z_DRAG = 1089.92, 1059.94        # Planck 2018, TT,TE,EE+lowE+lensing
THETA100, THETA100_ERR = 1.04110, 0.00031
WB = 0.02237                             # ω_b, Planck 2018 (то же значение для обеих моделей)
PLANCK_LCDM = dict(Om=0.3153, h=0.6736)  # для калибровки приближения r_s

# ==============================================================================
#  ДАННЫЕ, ВШИТЫЕ В ФАЙЛ
# ==============================================================================
# --- DESI DR2 BAO: официальный вектор данных и ковариация (13 точек).
#     Источник: CobayaSampler/bao_data, desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_{mean,cov}.txt
#     (то, чем сама коллаборация DESI пользуется в космологических подгонках;
#     arXiv:2503.14738, таблица IV — те же измерения, округлённые).
DESI_Z = np.array([0.295, 0.510, 0.510, 0.706, 0.706, 0.934, 0.934,
                   1.321, 1.321, 1.484, 1.484, 2.330, 2.330])
DESI_KIND = ["DV", "DM", "DH", "DM", "DH", "DM", "DH", "DM", "DH", "DM", "DH", "DH", "DM"]
DESI_VAL = np.array([7.94167639, 13.58758434, 21.86294686, 17.35069094, 19.45534918,
                     21.57563956, 17.64149464, 27.60085612, 14.17602155, 30.51190063,
                     12.81699964, 8.631545674846294, 38.988973961958784])
DESI_COV = np.zeros((13, 13))
DESI_COV[0, 0] = 5.78998687e-03
for (i, j, a, b_, c_) in [(1, 2, 2.83473742e-02, 1.83928040e-01, -3.26062007e-02),
                          (3, 4, 3.23752442e-02, 1.11469198e-01, -2.37445646e-02),
                          (5, 6, 2.61732816e-02, 4.04183878e-02, -1.12938006e-02),
                          (7, 8, 1.05336516e-01, 5.04233092e-02, -2.90308418e-02),
                          (9, 10, 5.83020277e-01, 2.68336193e-01, -1.95215562e-01),
                          (11, 12, 1.02136194e-02, 2.82685779e-01, -2.31395216e-02)]:
    DESI_COV[i, i], DESI_COV[j, j], DESI_COV[i, j], DESI_COV[j, i] = a, b_, c_, c_
DESI_TRACER = {0.295: "BGS", 0.510: "LRG1", 0.706: "LRG2", 0.934: "LRG3+ELG1",
               1.321: "ELG2", 1.484: "QSO", 2.330: "Lya"}
DESI_URL = ("https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2/"
            "desi_gaussian_bao_ALL_GCcomb_")

# --- Космические хронометры H(z): 33 точки, z от 0.07 до 1.965, км/с/Мпк.
#     15 точек Moresco (BC03): Moresco et al. 2012, Moresco 2015, Moresco et al. 2016.
#     Остальные 18: Zhang et al. 2014 (0.07, 0.12, 0.20, 0.28); Jimenez et al. 2003 (0.09);
#     Simon et al. 2005 (0.17, 0.27, 0.40, 0.90, 1.30, 1.43, 1.53, 1.75);
#     Stern et al. 2010 (0.48, 0.88); Ratsimbazafy et al. 2017 (0.47);
#     Borghi et al. 2022 (0.75); Jiao et al. 2023 (0.80).
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
# Ковариация 15 точек Moresco по рецепту Moresco et al. 2020 (ApJ 898, 82):
# (z, H, полная ошибка, статистическая часть, вклад металличности), км/с/Мпк
CC_MORESCO_BC03 = [
    (0.1791, 74.91, 3.8069262, 3.8, 0.5), (0.1993, 74.96, 4.9001352, 4.9, 0.6),
    (0.3519, 82.78, 13.94843, 13.0, 4.8), (0.3802, 83.00, 13.54, 4.3, 12.9),
    (0.4004, 76.97, 10.18, 2.1, 10.0), (0.4247, 87.08, 11.24, 2.4, 11.0),
    (0.4497, 92.78, 12.9, 4.5, 12.1), (0.4783, 80.91, 9.044, 2.1, 8.8),
    (0.5929, 103.80, 12.49752, 11.6, 4.5), (0.6797, 91.60, 7.961872, 6.4, 4.3),
    (0.7812, 104.50, 12.19515, 9.4, 6.1), (0.8754, 125.10, 16.70085, 15.3, 6.0),
    (1.0370, 153.70, 19.6736, 13.6, 14.9), (1.3630, 160.00, 32.63, 23.07, 23.07),
    (1.9650, 186.50, 49.58, 35.05, 35.05),
]
# Относительные (%) систематические вклады на сетке по z: (z, IMF, stlib, mod, mod_ooo)
CC_SYST_MM20 = [
    (0.075, 0.47, 7.40, 15.86, 9.91), (0.125, 0.47, 7.40, 14.23, 6.98),
    (0.175, 0.47, 7.40, 13.34, 5.40), (0.225, 0.47, 7.40, 13.21, 5.40),
    (0.275, 0.47, 7.40, 13.29, 5.40), (0.325, 0.47, 7.40, 12.20, 5.40),
    (0.375, 0.47, 7.40, 12.99, 5.40), (0.425, 0.47, 7.40, 10.29, 6.20),
    (0.475, 0.46, 7.39, 8.91, 5.86), (0.525, 0.23, 7.40, 9.99, 6.51),
    (0.575, 0.28, 6.87, 10.09, 6.12), (0.625, 0.47, 6.65, 11.17, 6.21),
    (0.675, 0.47, 6.57, 11.12, 5.71), (0.725, 0.47, 5.90, 10.81, 5.16),
    (0.775, 0.45, 6.03, 10.75, 5.05), (0.825, 0.47, 6.10, 10.75, 5.05),
    (0.875, 0.47, 5.89, 9.08, 2.79), (0.925, 0.44, 5.80, 8.62, 3.70),
    (0.975, 0.40, 5.94, 7.32, 3.65), (1.025, 0.27, 6.07, 5.84, 3.37),
    (1.075, 0.20, 6.08, 6.02, 3.49), (1.125, 0.20, 6.07, 4.72, 2.33),
    (1.175, 0.19, 6.09, 4.31, 2.33), (1.225, 0.19, 6.09, 3.90, 2.33),
    (1.275, 0.19, 6.09, 3.90, 2.33), (1.325, 0.20, 6.09, 3.91, 2.34),
    (1.375, 0.19, 6.09, 3.90, 2.34), (1.425, 0.19, 6.09, 3.90, 2.33),
    (1.475, 0.20, 6.09, 3.91, 2.34),
]

# --- Рост структур fσ8: набор Gold-2017, 18 точек (Nesseris, Pantazis & Perivolaropoulos 2017,
#     arXiv:1703.10538; значения по таблице I из Sagredo, Nesseris & Sapone 2018, arXiv:1806.10822).
FS8 = np.array([
    [0.020, 0.428, 0.0465], [0.020, 0.398, 0.0650], [0.020, 0.314, 0.0480],
    [0.100, 0.370, 0.1300], [0.150, 0.490, 0.1450], [0.170, 0.510, 0.0600],
    [0.180, 0.360, 0.0900], [0.380, 0.440, 0.0600], [0.250, 0.3512, 0.0583],
    [0.370, 0.4602, 0.0378], [0.320, 0.384, 0.0950], [0.590, 0.488, 0.0600],
    [0.440, 0.413, 0.0800], [0.600, 0.390, 0.0630], [0.730, 0.437, 0.0720],
    [0.600, 0.550, 0.1200], [0.860, 0.400, 0.1100], [1.400, 0.482, 0.1160],
])
FS8_WIGGLEZ_IDX = [12, 13, 14]           # ковариация WiggleZ: Blake et al. 2012
FS8_WIGGLEZ_COV = 1e-3 * np.array([[6.400, 2.570, 0.000], [2.570, 3.969, 2.540], [0.000, 2.540, 5.184]])

# --- Ранние сверхмассивные чёрные дыры: (имя, z, масса в массах Солнца).
BH_OBJECTS = [
    ("UHZ1", 10.073, 4.0e7),             # Bogdan et al. 2024; Natarajan et al. 2024; Goulding et al. 2023
    ("GN-z11", 10.603, 1.6e6),           # Maiolino et al. 2024
    ("CEERS 1019", 8.679, 1.0e7),        # Larson et al. 2023
    ("J0313-1806", 7.642, 1.6e9),        # Wang et al. 2021
    ("J1342+0928", 7.541, 7.8e8),        # Banados et al. 2018
    ("J1120+0641", 7.085, 2.0e9),        # Mortlock et al. 2011
]
T_SALPETER = 0.045                        # млрд лет на e-складывание (предел Эддингтона, ε = 0.1)

# ==============================================================================
#  МОДЕЛИ
# ==============================================================================
MU_X2 = 3 - np.sqrt(5)                    # μ² в точке q = 0
MU_H2 = (8 - 2 * np.sqrt(7)) / 3          # μ² в точке q = 1/2
KAPPA_H2 = 4 + 2 * np.sqrt(7)             # κ² в той же точке


def b_of_mu(mu):
    """Наблюдаемый масштаб, формула (8): b = sqrt(1-μ²)/μ²."""
    return np.sqrt(1 - mu * mu) / (mu * mu)


def q_of_y(y):
    """Параметр замедления ядра как функция y = μ², формула (12)."""
    return (-y * y + 6 * y - 4) / (2 - y) ** 2


def G_of_k2(k2):
    """Масса связи ядра H²b³/k² как функция κ², формула (25)."""
    return (2 + k2) ** 2 * (1 + k2) ** 1.5 / k2 ** 3


class DSIC:
    """Полная модель DSIC: одна формула темпа (27) для всей истории.
    mu0 — сегодняшняя доля объектности, h = H0/100, wc — потолок уравнения состояния
    запертой объектности (принцип P6: wc = 0)."""

    def __init__(self, mu0, h, wc=0.0):
        self.mu0, self.h, self.wc = mu0, h, wc
        self.b0 = b_of_mu(mu0)
        self.K = (2 - mu0 ** 2) / (1 - mu0 ** 2)          # H0/k
        self.Or = W_R / h ** 2
        qc = (1 + 3 * wc) / 2                               # w_core = (2q - 1)/3 = wc
        yc = MU_H2 if wc == 0.0 else brentq(lambda y: q_of_y(y) - qc, 0.5, 1 - 1e-12)
        self.zc = self.b0 / b_of_mu(np.sqrt(yc)) - 1        # точка включения ядра
        self.Fc = float(self.Ecore(self.zc) ** 2 / (1 + self.zc) ** 3)
        self.Am = (1 - self.Or) * self.Fc                   # амплитуда пылевой фазы (при wc = 0)
        self.brk = self.zc

    def mu_of_z(self, z):
        B = self.b0 / (1 + np.asarray(z, float))
        s = np.sqrt(1 + 4 * B * B)
        return np.sqrt(2 / (1 + s))

    def Ecore(self, z):
        """Темп одного ядра, H_core/H0: формула (13)."""
        B = self.b0 / (1 + np.asarray(z, float))
        s = np.sqrt(1 + 4 * B * B)
        return s * (s + 1) / (2 * B * B) / self.K

    def kappa(self, z):
        """Степень сжатия κ = Q/P, формула (23); P/R = 2B/(s+1) — без потери точности при μ → 1."""
        B = self.b0 / (1 + np.asarray(z, float))
        s = np.sqrt(1 + 4 * B * B)
        return np.sqrt(2 / (1 + s)) / (2 * B / (s + 1))

    def F(self, z):
        z = np.asarray(z, float)
        zl = np.minimum(z, self.zc)
        Fl = self.Ecore(zl) ** 2 / (1 + zl) ** 3
        if self.wc == 0.0:
            return Fl
        return np.where(z <= self.zc, Fl, self.Fc * ((1 + z) / (1 + self.zc)) ** (3 * self.wc))

    def E(self, z):
        z = np.asarray(z, float)
        return np.sqrt((1 - self.Or) * (1 + z) ** 3 * self.F(z) + self.Or * (1 + z) ** 4)

    def Sigma(self, z):
        """Ускорение внутреннего времени: H_core/H."""
        return self.Ecore(z) / self.E(z)

    def Sigma_rule(self, z):
        """То же в замкнутом виде, формула (29)."""
        k2 = self.kappa(z) ** 2
        return np.sqrt(G_of_k2(k2) / G_of_k2(KAPPA_H2) / (1 - self.Or) / (1 + self.Or * (1 + z) / self.Am))

    def X(self, z):
        """Пространство, которое прибавляет ядро, формула (31)."""
        z = np.asarray(z, float)
        return self.E(z) ** 2 - self.Am * (1 + z) ** 3 - self.Or * (1 + z) ** 4


class LCDM:
    def __init__(self, Om, h):
        self.Om, self.h = Om, h
        self.Or = W_R / h ** 2
        self.brk = None

    def E(self, z):
        z = np.asarray(z, float)
        return np.sqrt(self.Om * (1 + z) ** 3 + self.Or * (1 + z) ** 4 + 1 - self.Om - self.Or)


class Variant:
    """Вариант позднего закона без ранней части (для проверок принципов)."""

    def __init__(self, Efun, h=0.7):
        self.E, self.h, self.brk = Efun, h, None


def E_theta_uniform(mu0):
    """Если бы равномерно тёк угол θ, а не μ: H ∝ (2 - μ²)/sqrt(1 - μ²), см. формулу (17)."""
    b0 = b_of_mu(mu0)
    Hf = lambda y: (2 - y) / np.sqrt(1 - y)

    def E(z):
        B = b0 / (1 + np.asarray(z, float))
        s = np.sqrt(1 + 4 * B * B)
        return Hf(2 / (1 + s)) / Hf(mu0 * mu0)
    return E


def E_no_clock(mu0):
    """Правило измерения без множителя часов: b = P/Q, тогда H = k/(1 - μ²)."""
    b0 = np.sqrt(1 - mu0 * mu0) / mu0

    def E(z):
        B = b0 / (1 + np.asarray(z, float))
        return ((1 + B * B) / (B * B)) / ((1 + b0 * b0) / (b0 * b0))
    return E


# ---------- расстояния, звуковой горизонт, возраст ----------
def DM(M, z):
    """Сопутствующее расстояние, Мпк."""
    f = lambda x: np.exp(x) / M.E(np.expm1(x))
    brk = getattr(M, "brk", None)
    pts = [np.log1p(brk)] if (brk is not None and 0 < brk < z) else None
    return C_KMS / (100 * M.h) * quad(f, 0, np.log1p(z), points=pts, limit=400, epsabs=0, epsrel=1e-10)[0]


def DH(M, z):
    return C_KMS / (100 * M.h * float(M.E(z)))


def r_sound(M, z, wb=WB):
    """Звуковой горизонт к эпохе z, формула (32), Мпк."""
    def f(x):
        zz = np.expm1(x)
        cs = C_KMS / np.sqrt(3 * (1 + 3 * wb / (4 * W_GAMMA) / (1 + zz)))
        return cs / (100 * M.h * M.E(zz)) * np.exp(x)
    return quad(f, np.log1p(z), np.log1p(1e10), limit=400, epsabs=0, epsrel=1e-10)[0]


def age(M, z):
    """Собственное время от нуля пространства до эпохи z, млрд лет."""
    f = lambda x: 1 / M.E(np.expm1(x))
    brk = getattr(M, "brk", None)
    pts = [np.log1p(brk)] if (brk is not None and brk > z) else None
    return GYR / (100 * M.h) * quad(f, np.log1p(z), np.log1p(1e13), points=pts, limit=400, epsabs=0, epsrel=1e-10)[0]


_P = LCDM(**PLANCK_LCDM)
RS_P, DM_P = r_sound(_P, Z_STAR), DM(_P, Z_STAR)
CAL = THETA100 / (100 * RS_P / DM_P)       # поправка приближения, одна для обеих моделей


def theta100(M, wb=WB):
    return 100 * r_sound(M, Z_STAR, wb) / DM(M, Z_STAR) * CAL


def solve_h_theta(make, wb=WB):
    """h, при котором угол реликта модели равен измеренному."""
    return brentq(lambda h: theta100(make(h), wb) - THETA100, 0.5, 0.9, xtol=1e-8)


def growth_factor(M, Am, z_start):
    """Линейный рост возмущений вещества в пределе ОТО от z_start до 0: D(0)/D(z_start)
    и функция, возвращающая (D, f) на сетке."""
    def rhs(x, y):
        z = np.exp(-x) - 1
        d = 1e-4
        dlnH = (np.log(M.E(np.exp(-(x + d)) - 1)) - np.log(M.E(np.exp(-(x - d)) - 1))) / (2 * d)
        return [y[1], -(2 + dlnH) * y[1] + 1.5 * Am * (1 + z) ** 3 / M.E(z) ** 2 * y[0]]
    x0 = -np.log(1 + z_start)
    sol = solve_ivp(rhs, [x0, 0], [1.0, 1.0], rtol=1e-9, atol=1e-12, method="LSODA", dense_output=True)
    return sol


# ==============================================================================
#  [1] ЯДРО: ТОЖДЕСТВА И ОСОБЫЕ ТОЧКИ
# ==============================================================================
def block_core():
    head("[1] Ядро: формулы (6)–(19), таблица 1")
    mu = np.linspace(0.05, 0.995, 2000)
    k, d = 1.0, 1e-6
    Hf = lambda m: k * (2 - m ** 2) / (1 - m ** 2)
    # dμ/dτ = -kμ; численно проверяем H = d ln b/dτ и q = -1 - (dH/dτ)/H²
    dlnb = (np.log(b_of_mu(mu + d)) - np.log(b_of_mu(mu - d))) / (2 * d) * (-k * mu)
    dH = (Hf(mu + d) - Hf(mu - d)) / (2 * d) * (-k * mu)
    q_num = -1 - dH / Hf(mu) ** 2
    e1 = np.max(np.abs(dlnb / Hf(mu) - 1))
    e2 = np.max(np.abs(q_num - q_of_y(mu ** 2)))
    print(f"  H = k(2-μ²)/(1-μ²): макс. отклонение численной производной {e1:.1e}")
    print(f"  q = (-μ⁴+6μ²-4)/(2-μ²)²: макс. отклонение {e2:.1e}")
    tau = np.linspace(0.01, 3, 500)
    e3 = np.max(np.abs(b_of_mu(np.exp(-tau)) ** 2 / (np.exp(4 * tau) - np.exp(2 * tau)) - 1))
    print(f"  b²(τ) = e^(4kτ) - e^(2kτ): макс. отклонение {e3:.1e}")
    mux, muh = np.sqrt(MU_X2), np.sqrt(MU_H2)
    print(f"  q = 0   при μ² = 3-√5      → μ = {mux:.6f}, q = {q_of_y(MU_X2):+.1e}")
    print(f"  q = 1/2 при μ² = (8-2√7)/3 → μ = {muh:.6f}, q-1/2 = {q_of_y(MU_H2) - 0.5:+.1e}")
    k2 = np.linspace(2, 30, 200001)
    k2min = k2[np.argmin(G_of_k2(k2))]
    print(f"  минимум G(κ): κ² = {k2min:.4f} (формула: 4+2√7 = {KAPPA_H2:.4f}); "
          f"μ² = κ²/(1+κ²) = {KAPPA_H2 / (1 + KAPPA_H2):.6f} = (8-2√7)/3 = {MU_H2:.6f}")
    # расстояние в замкнутом виде против численного интеграла
    M = DSIC(0.76, 0.674)
    core = Variant(M.Ecore, 0.674)
    errs = []
    for z in (0.1, 0.5, 1.0, 2.0, 5.0, 100.0):
        mue = M.mu_of_z(z)
        chi_closed = C_KMS / 67.4 * M.K * M.b0 * (np.sqrt(1 - 0.76 ** 2) - np.sqrt(1 - mue ** 2))
        errs.append(abs(DM(core, z) / chi_closed - 1))
    print(f"  χ(μ) в замкнутом виде (18) против интеграла c/H: макс. отклонение {max(errs):.1e}")
    hor = C_KMS / 67.4 * M.K * (1 - 0.76 ** 2) / 0.76 ** 2
    print(f"  горизонт ядра при μ₀ = 0.76, H₀ = 67.4: {hor:.0f} Мпк")
    print("\n  Зависимость особых точек от μ₀:")
    print("     μ₀      b₀      K=H₀/k   z_t(q=0)  z_h(q=½)  κ₀")
    tab = []
    for m0 in (0.740, 0.744, 0.750, 0.752, 0.756, 0.760, 0.761, 0.765, 0.770):
        b0 = b_of_mu(m0)
        zt, zh = b0 / b_of_mu(mux) - 1, b0 / b_of_mu(muh) - 1
        tab.append([m0, zt, zh])
        print(f"    {m0:.3f}  {b0:.4f}  {(2 - m0 ** 2) / (1 - m0 ** 2):.4f}   {zt:.3f}     {zh:.3f}   {m0 / np.sqrt(1 - m0 ** 2):.3f}")
    R["core"] = dict(mu_x=mux, mu_h=muh, horizon_mu076_H674=hor, zt_zh_table=tab,
                     err_H=e1, err_q=e2, err_b2=e3, err_chi=max(errs))
    # врезка про Лоренца
    m115 = float(M.mu_of_z(1.15))
    print(f"\n  Врезка §6: P₀/R = {np.sqrt(1 - 0.76 ** 2):.3f}, 1/μ₀ = {1 / 0.76:.3f}; μ(z=1.15) = {m115:.4f}; "
          f"c̃: {100 * ((0.76 / m115) ** 2 - 1):+.0f} %, α: {100 * ((m115 / 0.76) ** 2 - 1):+.0f} %")
    R["lorentz"] = dict(mu_z115=m115, c_change=(0.76 / m115) ** 2 - 1, alpha_change=(m115 / 0.76) ** 2 - 1)
    # время наблюдателя
    for H0 in (67.9,):
        invk = GYR / H0 * M.K
        print(f"  Часы наблюдателя: 1/k = {invk:.1f} млрд лет; μ уменьшается вдвое каждые {np.log(2) * invk:.1f} млрд лет; "
              f"до равенства полюсов (μ = 1/√2): {np.log(0.76 * np.sqrt(2)) * invk:.2f} млрд лет")
        R["clock"] = dict(inv_k=invk, halving=np.log(2) * invk, to_equality=np.log(0.76 * np.sqrt(2)) * invk)


# ==============================================================================
#  [2] РАННЯЯ ВСЕЛЕННАЯ: ФОРМУЛЫ (20)–(32)
# ==============================================================================
def block_early():
    head("[2] Ранняя Вселенная: формулы (20)–(32), таблицы 3–5")
    print(f"  Самопроверка на ΛCDM Planck (h = 0.6736, Ωₘ = 0.3153): r* = {RS_P:.2f} Мпк, "
          f"100θ* = {100 * RS_P / DM_P:.5f} (измерено 1.04110) → поправка {CAL:.6f}")
    rdP, ageP = r_sound(_P, Z_DRAG), age(_P, 0)
    print(f"  ΛCDM Planck: r_d = {rdP:.2f} Мпк, возраст {ageP:.3f} млрд лет, t(z=7.5) = {age(_P, 7.5):.3f}")
    R["lcdm_planck"] = dict(r_star=RS_P, r_d=rdP, age=ageP, t75=age(_P, 7.5), cal=CAL)
    print("\n  Производные величины DSIC при H₀ из угла реликта:")
    print("     μ₀     H₀     z_t    z_h    A_m     A_m·h²  z_eq   r*      r_d     возраст  t(7.5)  Σ∞")
    rows = []
    for m0 in (0.760, 0.752, 0.744):
        h = solve_h_theta(lambda hh: DSIC(m0, hh))
        M = DSIC(m0, h)
        zt = M.b0 / b_of_mu(np.sqrt(MU_X2)) - 1
        row = dict(mu0=m0, H0=100 * h, z_t=zt, z_h=M.zc, Am=M.Am, Amh2=M.Am * h * h, z_eq=M.Am / M.Or - 1,
                   r_star=r_sound(M, Z_STAR), r_d=r_sound(M, Z_DRAG), age=age(M, 0), t75=age(M, 7.5),
                   Sigma_inf=1 / (M.K * M.b0 ** 2 * np.sqrt(M.Or)))
        rows.append(row)
        print(f"    {m0:.3f}  {row['H0']:.2f}  {zt:.3f}  {M.zc:.3f}  {M.Am:.4f}  {row['Amh2']:.4f}  {row['z_eq']:.0f}  "
              f"{row['r_star']:.2f}  {row['r_d']:.2f}  {row['age']:.3f}   {row['t75']:.3f}   {row['Sigma_inf']:.2f}")
    R["early_table"] = rows
    h = solve_h_theta(lambda hh: DSIC(0.76, hh))
    M = DSIC(0.76, h)
    core = Variant(M.Ecore, 0.674)
    print(f"\n  μ₀ = 0.76: доля атомов в веществе ω_b/(A_m h²) = {WB / (M.Am * h * h):.3f}")
    # требование (20) и предел ядра (21)
    need_lo, need_hi = 0.980 / (THETA100 / 100 * DESI_VAL[12]), 0.983 / (THETA100 / 100 * DESI_VAL[12])
    hor = M.K * M.b0 * np.sqrt(1 - 0.76 ** 2)
    chi233 = M.K * M.b0 * (np.sqrt(1 - 0.76 ** 2) - np.sqrt(1 - M.mu_of_z(2.33) ** 2))
    print(f"  (20) нужное D_M(z*)/D_M(2.33) = {need_lo:.3f}–{need_hi:.3f}; (21) предел ядра = {hor / chi233:.3f}")
    print(f"  D_M(z*)/r_d, которое нужно: {0.9826 / (THETA100 / 100):.1f}; горизонт ядра / 147 Мпк: {C_KMS / 67.4 * hor / 147:.1f}")
    einf = float(M.Ecore(1e7) / (1e7 + 1) ** 2)
    Or674 = W_R / 0.674 ** 2
    print(f"  ядро у нуля пространства: E → {einf:.4f}(1+z)²; плотность в {einf ** 2 / Or674:.0f} раз больше реального света, "
          f"темп в {einf / np.sqrt(Or674):.1f} раза (при h = 0.674)")
    print(f"  темп при нуклеосинтезе (z = 4·10⁸) относительно стандартного: полная модель "
          f"{float(M.E(4e8)) / np.sqrt(M.Or * (4e8 + 1) ** 4):.4f}, одно ядро {float(M.Ecore(4e8)) / np.sqrt(M.Or * (4e8 + 1) ** 4):.1f}")
    print("\n  Сжатие и ускорение внутреннего времени (μ₀ = 0.76):")
    print("       z        κ        Σ по (27)    Σ по (29)")
    sig = []
    for z in (2.0, M.zc, 5, 10, 100, 1090, round(M.Am / M.Or - 1), 1e6, 4e8):
        s1, s2 = float(M.Sigma(z)), float(M.Sigma_rule(z)) if z >= M.zc else 1.0
        sig.append([z, float(M.kappa(z)), s1, s2])
        print(f"    {z:9.4g}  {float(M.kappa(z)):9.4g}  {s1:9.4f}   {s2:9.4f}")
    zz = np.logspace(np.log10(M.zc + 0.01), 8.6, 3000)
    dev = np.max(np.abs(M.Sigma_rule(zz) / M.Sigma(zz) - 1))
    print(f"  замкнутый вид (29) против определения через формулу (27): макс. отклонение {dev:.1e}")
    zl = np.linspace(0, M.zc, 500)
    print(f"  ниже точки включения |E/E_core - 1| ≤ {np.max(np.abs(M.E(zl) / M.Ecore(zl) - 1)):.1e}")
    R["sigma_table"] = sig
    R["kappa_eq_z"] = M.b0 / b_of_mu(1 / np.sqrt(2)) - 1
    print(f"  κ = 1 (равенство полюсов) при z = {R['kappa_eq_z']:.3f}")
    Rb = 3 * WB / (4 * W_GAMMA * (1 + Z_STAR))
    print(f"  скорость звука на рекомбинации: {1 / np.sqrt(3 * (1 + Rb)):.3f} c")
    # одно ядро
    age_core = -np.log(0.76) * M.K * GYR / 67.4
    t75_core = -np.log(float(M.mu_of_z(7.5))) * M.K * GYR / 67.4
    print(f"  одно ядро (H₀ = 67.4): возраст {age_core:.2f}, t(z=7.5) = {t75_core:.2f} млрд лет")
    g = growth_factor(M, M.Am, Z_STAR).y[0, -1]
    gP = growth_factor(_P, 0.3153, Z_STAR).y[0, -1]
    print(f"  рост возмущений D(0)/D(z*): DSIC ×{g:.0f}, ΛCDM Planck ×{gP:.0f}")
    R["early_misc"] = dict(need=[need_lo, need_hi], core_limit=hor / chi233, e_inf=einf, sigma_dev=dev,
                           age_core=age_core, t75_core=t75_core, growth_dsic=g, growth_lcdm=gP,
                           cs_rec=1 / np.sqrt(3 * (1 + Rb)), baryon_frac=WB / (M.Am * h * h))
    # сгусток, «застывший» на сжатии рекомбинации
    f = lambda x: float(M.Sigma(Z_STAR) / M.Sigma(np.expm1(x)) / M.E(np.expm1(x)))
    t_clump = GYR / (100 * h) * quad(f, np.log1p(7.3), np.log1p(Z_STAR), limit=400)[0]
    print(f"  сгусток с замороженным Σ рекомбинации: к z = 7.3 прожил {t_clump:.2f} млрд лет против {age(M, 7.3):.2f} у фона")
    R["clump"] = dict(t_clump=t_clump, t_bg=age(M, 7.3))
    return M


# ==============================================================================
#  [3] ЯЗЫК ТЁМНОЙ ЭНЕРГИИ И Om(z)
# ==============================================================================
def block_de(M):
    head("[3] Язык тёмной энергии: формулы (31) и (34), μ₀ = 0.76")

    def wX(z, d=1e-4):
        return -1 + (np.log(M.X(z + d)) - np.log(M.X(z - d))) / (np.log1p(z + d) - np.log1p(z - d)) / 3
    w0 = float(wX(1e-3))
    zx = brentq(lambda z: wX(z) + 1, 0.05, 1.0)
    print(f"  w_X(0) = {w0:.3f}; пересечение w = -1 при z = {zx:.3f}")
    print(f"  ρ_X(1)/ρ_X(0) = {float(M.X(1) / M.X(0)):.3f}; ρ_X(2)/ρ_X(0) = {float(M.X(2) / M.X(0)):.3f}")
    cpl = {}
    for zmax in (1.0, 1.5):
        zg = np.linspace(0.001, zmax, 400)
        w = np.array([wX(z) for z in zg])
        A = np.vstack([np.ones_like(zg), zg / (1 + zg)]).T
        a0, a1 = np.linalg.lstsq(A, w, rcond=None)[0]
        cpl[zmax] = (a0, a1)
        print(f"  аппроксимация w₀ + w_a(1-a) на z ≤ {zmax}: w₀ = {a0:.2f}, w_a = {a1:.2f}")
    Om = lambda Mod, z: float((Mod.E(z) ** 2 - 1) / ((1 + z) ** 3 - 1))
    print(f"  Om(0.2) = {Om(M, 0.2):.3f}, Om(1) = {Om(M, 1.0):.3f}, Om(2) = {Om(M, 2.0):.3f}")
    R["de"] = dict(w0=w0, z_cross=zx, X1=float(M.X(1) / M.X(0)), X2=float(M.X(2) / M.X(0)),
                   cpl_z1=cpl[1.0], cpl_z15=cpl[1.5], Om02=Om(M, 0.2), Om1=Om(M, 1.0), Om2=Om(M, 2.0))


# ==============================================================================
#  ЗАГРУЗКА ДАННЫХ
# ==============================================================================
SN_BASE = ("https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
           "Pantheon%2B_Data/4_DISTANCES_AND_COVAR/")


def fetch(url, name):
    path = os.path.join(DATA_DIR, name)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"  скачиваю {name} ...")
        urllib.request.urlretrieve(url, path)
    return path


def load_desi_check():
    head("[4] Данные DESI DR2: сверка вшитой копии с официальным файлом")
    try:
        pm = fetch(DESI_URL + "mean.txt", "desi_dr2_mean.txt")
        pc = fetch(DESI_URL + "cov.txt", "desi_dr2_cov.txt")
        with open(pm, encoding="utf-8") as fh:
            rows = [l.split() for l in fh if not l.startswith("#")]
        z = np.array([float(r[0]) for r in rows])
        v = np.array([float(r[1]) for r in rows])
        kind = [r[2][:2] for r in rows]
        cov = np.loadtxt(pc)
        ok = (np.allclose(z, DESI_Z) and np.allclose(v, DESI_VAL, rtol=1e-9) and kind == DESI_KIND
              and np.allclose(cov, DESI_COV, rtol=1e-8, atol=1e-14))
        print("  официальный файл совпадает с вшитой копией:", "ДА" if ok else "НЕТ — ПРОВЕРЬТЕ ДАННЫЕ")
        R["desi_official_match"] = bool(ok)
    except Exception as ex:                                  # нет сети — работаем по вшитой копии
        print("  не удалось скачать официальный файл, использую вшитую копию:", ex)
        R["desi_official_match"] = None
    sd = np.sqrt(np.diag(DESI_COV))
    print("  трассер      z      величина   значение   ошибка   корреляция D_M,D_H")
    for i in range(13):
        corr = ""
        if DESI_KIND[i] == "DM":
            j = i + 1 if i + 1 < 13 and DESI_Z[i + 1] == DESI_Z[i] else i - 1
            corr = f"{DESI_COV[i, j] / (sd[i] * sd[j]):+.3f}"
        print(f"  {DESI_TRACER[round(float(DESI_Z[i]), 3)]:10s} {DESI_Z[i]:.3f}   {DESI_KIND[i]}/r_d    {DESI_VAL[i]:8.3f}  {sd[i]:7.3f}   {corr}")


def load_sn():
    head("[5] Сверхновые Pantheon+SH0ES: загрузка")
    pd = fetch(SN_BASE + "Pantheon%2BSH0ES.dat", "PantheonPlusSH0ES.dat")
    pc = fetch(SN_BASE + "Pantheon%2BSH0ES_STAT%2BSYS.cov", "PantheonPlusSH0ES_STAT+SYS.cov")
    ps = fetch(SN_BASE + "Pantheon%2BSH0ES_STATONLY.cov", "PantheonPlusSH0ES_STATONLY.cov")
    with open(pd, encoding="utf-8") as fh:
        names = fh.readline().split()
    cols = [names.index(n) for n in ("zHD", "zHEL", "m_b_corr", "CEPH_DIST", "IS_CALIBRATOR")]
    zHD, zHEL, mb, ceph, cal = np.genfromtxt(pd, skip_header=1, usecols=cols).T

    def read_cov(p):
        raw = np.loadtxt(p)
        n = int(raw[0])
        return raw[1:].reshape(n, n)
    Ctot, Cstat = read_cov(pc), read_cov(ps)

    def subset(mask, C):
        return dict(z=zHD[mask], zh=zHEL[mask], mb=mb[mask], ceph=ceph[mask], cal=cal[mask].astype(bool),
                    Ci=np.linalg.inv(C[np.ix_(mask, mask)]))
    m_all = (zHD > 0.01) | (cal == 1)
    m_hf = (zHD > 0.01) & (cal == 0)
    SN = dict(all=subset(m_all, Ctot), hf=subset(m_hf, Ctot), stat=subset(m_all, Cstat))
    print(f"  всего строк {len(zHD)}; с калибраторами SH0ES: {m_all.sum()}; только хаббловский поток: {m_hf.sum()}; "
          f"калибраторов {int(cal[m_all].sum())}")
    R["sn_counts"] = dict(total=len(zHD), with_cal=int(m_all.sum()), hubble_flow=int(m_hf.sum()))
    return SN


# ==============================================================================
#  ФУНКЦИИ χ²
# ==============================================================================
_ZG = np.concatenate([[0], np.logspace(-4, np.log10(2.6), 6000)])


def mu_theory_h1(Efun, S):
    """Модуль расстояния при h = 1: 5 lg[(1+z_hel)·D_M(z_HD)/Мпк] + 25."""
    f = 1 / Efun(_ZG)
    I = np.concatenate([[0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(_ZG))])
    dm = C_KMS / 100 * np.interp(S["z"], _ZG, I)
    return 5 * np.log10(np.maximum((1 + S["zh"]) * dm, 1e-12)) + 25


def chi2_sn_cal(Efun, S, return_lin=False):
    """SN с калибраторами: M и η = 5 lg h входят линейно и находятся аналитически.
    Остаток: калибратор m - M - μ_цефеид; поток m - M - μ_теор(h=1) + η."""
    hf = (~S["cal"]).astype(float)
    d = np.where(S["cal"], S["mb"] - S["ceph"], S["mb"] - mu_theory_h1(Efun, S))
    A = np.vstack([np.ones_like(d), -hf]).T            # d = M·1 - η·hf + остаток
    N = A.T @ S["Ci"] @ A
    p = np.linalg.solve(N, A.T @ S["Ci"] @ d)
    r = d - A @ p
    c2 = float(r @ S["Ci"] @ r)
    if return_lin:
        return c2, p[0], 10 ** (p[1] / 5), np.linalg.inv(N)     # χ², M, h, ков. (M, η)
    return c2


def chi2_sn_marg(Efun, S):
    """SN без калибраторов, аналитическая маргинализация по M + 5 lg h."""
    r = S["mb"] - mu_theory_h1(Efun, S)
    one = np.ones_like(r)
    a, b = one @ S["Ci"] @ r, one @ S["Ci"] @ one
    return float(r @ S["Ci"] @ r - a * a / b)


DESI_CI = np.linalg.inv(DESI_COV)
DESI_CI12 = np.linalg.inv(DESI_COV[1:, 1:])


def bao_vector(M):
    out = []
    for z, kind in zip(DESI_Z, DESI_KIND):
        dm, dh = DM(M, z), DH(M, z)
        out.append(dm if kind == "DM" else dh if kind == "DH" else (z * dm * dm * dh) ** (1 / 3))
    return np.array(out)


def chi2_bao(M, rd=None, with_bgs=True, pulls=False):
    """BAO DESI DR2. rd=None — r_d (точнее, r_d·h) свободен и находится аналитически."""
    t = bao_vector(M)
    Ci, d = (DESI_CI, DESI_VAL) if with_bgs else (DESI_CI12, DESI_VAL[1:])
    if not with_bgs:
        t = t[1:]
    r = (t @ Ci @ t) / (t @ Ci @ d) if rd is None else rd
    res = t / r - d
    c2 = float(res @ Ci @ res)
    if pulls:
        return c2, r, res / np.sqrt(np.diag(DESI_COV if with_bgs else DESI_COV[1:, 1:]))
    return c2, r


def chi2_theta(M, wb=WB):
    return ((theta100(M, wb) - THETA100) / THETA100_ERR) ** 2


def hessian_err(fun, x, steps):
    n = len(x)
    Hm = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            ei, ej = np.eye(n)[i] * steps[i], np.eye(n)[j] * steps[j]
            Hm[i, j] = (fun(x + ei + ej) - fun(x + ei - ej) - fun(x - ei + ej) + fun(x - ei - ej)) / (4 * steps[i] * steps[j])
    return np.sqrt(np.diag(np.linalg.inv(0.5 * Hm)))


def profile_err(fun, x0, c0, step):
    """Ошибка одного параметра по профилю: Δχ² = 1 (среднее по двум сторонам)."""
    out = []
    for sgn in (+1, -1):
        g = lambda d: fun(x0 + sgn * d) - c0 - 1.0
        hi = step
        while g(hi) < 0 and hi < 50 * step:
            hi *= 2
        out.append(brentq(g, 1e-6 * step, hi))
    return 0.5 * (out[0] + out[1])


# ==============================================================================
#  [6] ПОЗДНЯЯ ВСЕЛЕННАЯ: ТАБЛИЦА 1
# ==============================================================================
H_LATE = 0.70          # h в члене Ω_r для поздних проб; влияние на χ² < 0.01


def block_late(SN):
    head("[6] Поздняя Вселенная (таблица 2)")
    mk_d = lambda m: DSIC(m, H_LATE)
    mk_l = lambda o: LCDM(o, H_LATE)
    out = {}
    # --- сверхновые с калибраторами
    for tag, key, S in (("STAT+SYS", "sn", SN["all"]), ("STAT-ONLY", "sn_stat", SN["stat"])):
        fd = lambda m: chi2_sn_cal(mk_d(m).E, S)
        fl = lambda o: chi2_sn_cal(mk_l(o).E, S)
        rd_ = minimize_scalar(fd, bounds=(0.60, 0.90), method="bounded", options=dict(xatol=1e-7))
        rl_ = minimize_scalar(fl, bounds=(0.10, 0.60), method="bounded", options=dict(xatol=1e-7))
        cd, Md, hd, covd = chi2_sn_cal(mk_d(rd_.x).E, S, True)
        cl, Ml, hl, covl = chi2_sn_cal(mk_l(rl_.x).E, S, True)
        ed, el = profile_err(fd, rd_.x, cd, 0.01), profile_err(fl, rl_.x, cl, 0.02)

        # ошибка H0: полный гессиан по (форма, η, M)
        def full(make, S=S):
            def f(p):
                hf = (~S["cal"]).astype(float)
                d = np.where(S["cal"], S["mb"] - S["ceph"], S["mb"] - mu_theory_h1(make(p[0]).E, S))
                r = d - p[2] + p[1] * hf
                return float(r @ S["Ci"] @ r)
            return f
        e3d = hessian_err(full(mk_d), np.array([rd_.x, 5 * np.log10(hd), Md]), [2e-3, 2e-3, 2e-3])
        e3l = hessian_err(full(mk_l), np.array([rl_.x, 5 * np.log10(hl), Ml]), [4e-3, 2e-3, 2e-3])
        eHd, eHl = 100 * hd * np.log(10) / 5 * e3d[1], 100 * hl * np.log(10) / 5 * e3l[1]
        print(f"  Сверхновые {tag} ({len(S['z'])} SN):")
        print(f"    DSIC : μ₀ = {rd_.x:.4f} ± {ed:.4f}, H₀ = {100 * hd:.2f} ± {eHd:.2f}, M = {Md:.4f}, χ² = {cd:.2f}")
        print(f"    ΛCDM : Ωₘ = {rl_.x:.4f} ± {el:.4f}, H₀ = {100 * hl:.2f} ± {eHl:.2f}, M = {Ml:.4f}, χ² = {cl:.2f}")
        print(f"    Δχ² (DSIC - ΛCDM) = {cd - cl:+.2f}")
        zt = b_of_mu(rd_.x) / b_of_mu(np.sqrt(MU_X2)) - 1
        zt_err = abs(b_of_mu(rd_.x + ed) - b_of_mu(rd_.x - ed)) / 2 / b_of_mu(np.sqrt(MU_X2))
        out[key] = dict(mu0=rd_.x, mu0_err=ed, H0_d=100 * hd, H0_d_err=eHd, chi2_d=cd,
                        Om=rl_.x, Om_err=el, H0_l=100 * hl, H0_l_err=eHl, chi2_l=cl, dchi2=cd - cl,
                        n=len(S["z"]), zt=zt, zt_err=zt_err, zh=b_of_mu(rd_.x) / b_of_mu(np.sqrt(MU_H2)) - 1)
    S = SN["all"]
    c_l_sn = out["sn"]["chi2_l"]
    # --- проверки принципов на сверхновых
    ft = lambda m: chi2_sn_cal(E_theta_uniform(m), S)
    rt = minimize_scalar(ft, bounds=(0.80, 0.995), method="bounded", options=dict(xatol=1e-7))
    fn = lambda m: chi2_sn_cal(E_no_clock(m), S)
    rn = minimize_scalar(fn, bounds=(0.30, 0.80), method="bounded", options=dict(xatol=1e-7))
    print(f"  Равномерный угол θ вместо μ: μ₀ = {rt.x:.4f}, χ² = {rt.fun:.2f}, Δχ² к ΛCDM = {rt.fun - c_l_sn:+.2f}")
    print(f"  Без множителя часов (b = P/Q), только сверхновые: μ₀ = {rn.x:.4f}, χ² = {rn.fun:.2f}, "
          f"Δχ² к DSIC = {rn.fun - out['sn']['chi2_d']:+.2f}")
    out["theta_uniform"] = dict(mu0=rt.x, chi2=rt.fun, dchi2_lcdm=rt.fun - c_l_sn)
    out["no_clock_sn"] = dict(mu0=rn.x, chi2=rn.fun, dchi2_dsic=rn.fun - out["sn"]["chi2_d"])
    # --- сверхновые + BAO со свободным r_d
    res = {}
    for nb, with_bgs in ((13, True), (12, False)):
        tot_d = lambda m: chi2_sn_cal(mk_d(m).E, S) + chi2_bao(mk_d(m), None, with_bgs)[0]
        tot_l = lambda o: chi2_sn_cal(mk_l(o).E, S) + chi2_bao(mk_l(o), None, with_bgs)[0]
        tot_n = lambda m: chi2_sn_cal(E_no_clock(m), S) + chi2_bao(Variant(E_no_clock(m)), None, with_bgs)[0]
        a = minimize_scalar(tot_d, bounds=(0.70, 0.82), method="bounded", options=dict(xatol=1e-7))
        b = minimize_scalar(tot_l, bounds=(0.20, 0.42), method="bounded", options=dict(xatol=1e-7))
        n = minimize_scalar(tot_n, bounds=(0.35, 0.70), method="bounded", options=dict(xatol=1e-7))
        ba, bb = chi2_bao(mk_d(a.x), None, with_bgs)[0], chi2_bao(mk_l(b.x), None, with_bgs)[0]
        ea, eb = profile_err(tot_d, a.x, a.fun, 0.005), profile_err(tot_l, b.x, b.fun, 0.01)
        print(f"  Сверхновые + BAO ({nb} точек DESI), r_d свободен:")
        print(f"    DSIC : μ₀ = {a.x:.4f} ± {ea:.4f}, χ²(SN) = {a.fun - ba:.2f}, χ²(BAO) = {ba:.2f}, итого {a.fun:.2f}")
        print(f"    ΛCDM : Ωₘ = {b.x:.4f} ± {eb:.4f}, χ²(SN) = {b.fun - bb:.2f}, χ²(BAO) = {bb:.2f}, итого {b.fun:.2f}")
        print(f"    Δχ² (DSIC - ΛCDM) = {a.fun - b.fun:+.2f};  без множителя часов: μ₀ = {n.x:.4f}, Δχ² к DSIC = {n.fun - a.fun:+.2f}")
        res[f"n{nb}"] = dict(mu0=a.x, mu0_err=ea, Om=b.x, Om_err=eb, bao_d=ba, bao_l=bb, tot_d=a.fun, tot_l=b.fun,
                       dchi2=a.fun - b.fun, noclock_mu0=n.x, noclock_dchi2=n.fun - a.fun)
    out["sn_bao"] = res
    mu_ref, om_ref = res["n13"]["mu0"], res["n13"]["Om"]
    Md_, Ml_ = mk_d(mu_ref), mk_l(om_ref)
    # --- только BAO, свободный r_d
    a = minimize_scalar(lambda m: chi2_bao(mk_d(m))[0], bounds=(0.70, 0.82), method="bounded")
    b = minimize_scalar(lambda o: chi2_bao(mk_l(o))[0], bounds=(0.20, 0.42), method="bounded")
    print(f"  Только BAO (13 точек), r_d свободен: DSIC μ₀ = {a.x:.4f}, χ² = {a.fun:.2f}; ΛCDM Ωₘ = {b.x:.4f}, χ² = {b.fun:.2f}; Δχ² = {a.fun - b.fun:+.2f}")
    out["bao_only"] = dict(mu0=a.x, chi2_d=a.fun, Om=b.x, chi2_l=b.fun, dchi2=a.fun - b.fun)
    # --- отношение D_M/D_H: не зависит ни от r_d, ни от H0
    def chi_ratio(M):
        s = 0.0
        for i in range(13):
            if DESI_KIND[i] != "DM":
                continue
            j = i + 1 if (i + 1 < 13 and DESI_Z[i + 1] == DESI_Z[i]) else i - 1
            dm, dh = DESI_VAL[i], DESI_VAL[j]
            r = dm / dh
            sr = r * np.sqrt(DESI_COV[i, i] / dm ** 2 + DESI_COV[j, j] / dh ** 2 - 2 * DESI_COV[i, j] / (dm * dh))
            s += ((DM(M, DESI_Z[i]) / DH(M, DESI_Z[i]) - r) / sr) ** 2
        return s
    cr_d, cr_l = chi_ratio(Md_), chi_ratio(Ml_)
    print(f"  Отношение D_M/D_H, 6 точек (μ₀ = {mu_ref:.4f}, Ωₘ = {om_ref:.4f}): DSIC χ² = {cr_d:.2f}, ΛCDM χ² = {cr_l:.2f}, Δχ² = {cr_d - cr_l:+.2f}")
    out["ratio"] = dict(chi2_d=cr_d, chi2_l=cr_l, dchi2=cr_d - cr_l)
    # --- хронометры: форма закреплена теми же μ₀ и Ωₘ, H0 находится аналитически
    zc, Hc, sc = CC.T
    Cdiag = np.diag(sc ** 2)
    Cfull = cc_covariance(zc, sc)
    cc = {}
    for tag, key, Cm in (("диагональные ошибки", "diag", Cdiag), ("ковариация Moresco 2020", "cov", Cfull)):
        Ci = np.linalg.inv(Cm)
        row = {}
        for name, M in (("DSIC", Md_), ("ΛCDM", Ml_)):
            e = M.E(zc)
            H0 = (e @ Ci @ Hc) / (e @ Ci @ e)
            r = H0 * e - Hc
            row[name] = (float(r @ Ci @ r), float(H0), float(1 / np.sqrt(e @ Ci @ e)))
        print(f"  Хронометры, 33 точки, {tag}: DSIC χ² = {row['DSIC'][0]:.2f} (H₀ = {row['DSIC'][1]:.1f} ± {row['DSIC'][2]:.1f}); "
              f"ΛCDM χ² = {row['ΛCDM'][0]:.2f} (H₀ = {row['ΛCDM'][1]:.1f} ± {row['ΛCDM'][2]:.1f}); Δχ² = {row['DSIC'][0] - row['ΛCDM'][0]:+.2f}")
        cc[key] = dict(chi2_d=row["DSIC"][0], H0_d=row["DSIC"][1], chi2_l=row["ΛCDM"][0], H0_l=row["ΛCDM"][1],
                       H0_err_d=row["DSIC"][2], H0_err_l=row["ΛCDM"][2], dchi2=row["DSIC"][0] - row["ΛCDM"][0])
    out["cc"] = cc
    out["ref"] = dict(mu0=mu_ref, Om=om_ref)
    zt_l = (2 * (1 - om_ref) / om_ref) ** (1 / 3) - 1
    print(f"  Переход к ускорению: DSIC z_t = {b_of_mu(mu_ref) / b_of_mu(np.sqrt(MU_X2)) - 1:.3f} (μ₀ = {mu_ref:.4f}); "
          f"ΛCDM z_t = {zt_l:.3f} (Ωₘ = {om_ref:.4f}), {(2 * (1 - 0.3153) / 0.3153) ** (1 / 3) - 1:.3f} (Planck)")
    out["zt"] = dict(dsic=b_of_mu(mu_ref) / b_of_mu(np.sqrt(MU_X2)) - 1, lcdm=zt_l,
                     lcdm_planck=(2 * (1 - 0.3153) / 0.3153) ** (1 / 3) - 1)
    R["late"] = out


def cc_covariance(z_cc, sH_cc):
    """Ковариация 33 хронометров: блок 15 точек Moresco — по рецепту Moresco et al. 2020
    (диагональ = полная ошибка², вне диагонали — коррелированная систематика), остальные диагональны."""
    C = np.diag(np.asarray(sH_cc, float) ** 2).copy()
    mor = np.array(CC_MORESCO_BC03, float)
    z_m, err_t, met = mor[:, 0], mor[:, 2], mor[:, 4]
    syst = np.array(CC_SYST_MM20, float)
    frac = np.vstack([np.interp(z_m, syst[:, 0], syst[:, k]) / 100.0 for k in (1, 2, 3, 4)])
    quad_ = np.sqrt((frac ** 2).sum(axis=0))
    quad_[quad_ == 0] = 1.0
    sig = (frac / quad_) * met[None, :]
    C_syst = sum(np.outer(sig[s], sig[s]) for s in range(4))
    idx = [int(np.argmin(np.abs(z_cc - zm))) for zm in z_m]
    for a in range(15):
        for b_ in range(15):
            C[idx[a], idx[b_]] = (err_t[a] ** 2 if a == b_ else 0.0) + C_syst[a, b_]
    return C


# ==============================================================================
#  [7] СОВМЕСТНЫЙ ТЕСТ: ТАБЛИЦА 3
# ==============================================================================
def block_joint(SN):
    head("[7] Совместный тест (таблица 6): SN без калибраторов + BAO со своим r_d + угол реликта")
    S = SN["hf"]

    def parts(M, with_bgs=True):
        rd = r_sound(M, Z_DRAG)
        cs, cb, ct = chi2_sn_marg(M.E, S), chi2_bao(M, rd, with_bgs)[0], chi2_theta(M)
        return cs, cb, ct, rd

    def fit(make, x0, with_bgs=True, bounds=None):
        def f(p):
            if bounds is not None and any((v < lo or v > hi) for v, (lo, hi) in zip(p, bounds)):
                return 1e30
            return sum(parts(make(p), with_bgs)[:3])
        r = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-6, maxiter=4000))
        r = minimize(f, r.x, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-7, maxiter=4000))
        return r, f
    out = {}
    for nb, with_bgs in ((13, True), (12, False)):
        print(f"\n  --- DESI: {nb} точек ---")
        rl, fl = fit(lambda p: LCDM(p[0], p[1]), [0.296, 0.687], with_bgs, [(0.1, 0.6), (0.5, 0.9)])
        rd_, fd = fit(lambda p: DSIC(p[0], p[1]), [0.752, 0.691], with_bgs, [(0.6, 0.9), (0.5, 0.9)])
        rw, fw = fit(lambda p: DSIC(p[0], p[1], p[2]), [0.761, 0.63, -0.013], with_bgs,
                     [(0.6, 0.9), (0.5, 0.9), (-0.2, 0.2)])
        pl, pd_, pw = parts(LCDM(*rl.x), with_bgs), parts(DSIC(*rd_.x), with_bgs), parts(DSIC(*rw.x), with_bgs)
        el = hessian_err(fl, rl.x, [2e-3, 2e-3])
        ed = hessian_err(fd, rd_.x, [1e-3, 2e-3])
        ew = hessian_err(fw, rw.x, [1e-3, 2e-3, 1e-3])
        print("  модель                параметры                                  χ²SN     χ²BAO   χ²θ*   итого    r_d")
        print(f"  ΛCDM                  Ωₘ = {rl.x[0]:.4f} ± {el[0]:.4f}, H₀ = {100 * rl.x[1]:.2f} ± {100 * el[1]:.2f}          "
              f"{pl[0]:8.2f} {pl[1]:7.2f} {pl[2]:6.2f} {rl.fun:8.2f}  {pl[3]:.2f}")
        print(f"  DSIC, P6 (w_c = 0)    μ₀ = {rd_.x[0]:.4f} ± {ed[0]:.4f}, H₀ = {100 * rd_.x[1]:.2f} ± {100 * ed[1]:.2f}          "
              f"{pd_[0]:8.2f} {pd_[1]:7.2f} {pd_[2]:6.2f} {rd_.fun:8.2f}  {pd_[3]:.2f}")
        print(f"  DSIC, w_c свободен    μ₀ = {rw.x[0]:.4f}, H₀ = {100 * rw.x[1]:.2f} ± {100 * ew[1]:.2f}, w_c = {rw.x[2]:+.4f} ± {ew[2]:.4f}  "
              f"{pw[0]:8.2f} {pw[1]:7.2f} {pw[2]:6.2f} {rw.fun:8.2f}  {pw[3]:.2f}")
        print(f"  Δχ²: DSIC - ΛCDM = {rd_.fun - rl.fun:+.2f};  свободный w_c - ΛCDM = {rw.fun - rl.fun:+.2f};  "
              f"w_c = 0 против свободного: {rd_.fun - rw.fun:+.2f}")
        out[f"n{nb}"] = dict(
            lcdm=dict(Om=rl.x[0], Om_err=el[0], H0=100 * rl.x[1], H0_err=100 * el[1], sn=pl[0], bao=pl[1], th=pl[2], tot=rl.fun, rd=pl[3]),
            dsic=dict(mu0=rd_.x[0], mu0_err=ed[0], H0=100 * rd_.x[1], H0_err=100 * ed[1], sn=pd_[0], bao=pd_[1], th=pd_[2], tot=rd_.fun, rd=pd_[3]),
            wc=dict(mu0=rw.x[0], H0=100 * rw.x[1], H0_err=100 * ew[1], wc=rw.x[2], wc_err=ew[2], sn=pw[0], bao=pw[1], th=pw[2], tot=rw.fun, rd=pw[3]),
            dchi2=rd_.fun - rl.fun, dchi2_wc=rw.fun - rl.fun, dchi2_wc0=rd_.fun - rw.fun)
        Mj = DSIC(*rd_.x)
        out[f"n{nb}"]["dsic"].update(zt=Mj.b0 / b_of_mu(np.sqrt(MU_X2)) - 1, zh=Mj.zc, Am=Mj.Am, age=age(Mj, 0))
        out[f"n{nb}"]["lcdm"].update(age=age(LCDM(*rl.x), 0))
        if with_bgs:
            for name, M in (("DSIC", DSIC(*rd_.x)), ("ΛCDM", LCDM(*rl.x))):
                _, _, pulls = chi2_bao(M, r_sound(M, Z_DRAG), True, True)
                txt = ", ".join(f"{DESI_TRACER[round(float(z), 3)]} {k}: {p:+.1f}" for z, k, p in zip(DESI_Z, DESI_KIND, pulls))
                print(f"  отклонения точек BAO от {name}, в σ (модель минус измерение): {txt}")
                out["n13"]["pulls_" + ("dsic" if name == "DSIC" else "lcdm")] = [float(p) for p in pulls]
    # r_d, который требуют BAO, при μ₀ = 0.76 и у ΛCDM Planck
    h = solve_h_theta(lambda hh: DSIC(0.76, hh))
    M = DSIC(0.76, h)
    need_d = chi2_bao(M, None)[1]
    need_l = chi2_bao(_P, None)[1]
    hb = brentq(lambda hh: r_sound(DSIC(0.76, hh), Z_DRAG) - chi2_bao(DSIC(0.76, hh), None)[1], 0.6, 0.8)
    print(f"\n  μ₀ = 0.76, H₀ = {100 * h:.2f}: r_d модели {r_sound(M, Z_DRAG):.2f} Мпк, BAO требуют {need_d:.2f} "
          f"({100 * (1 - r_sound(M, Z_DRAG) / need_d):.1f} % меньше); H₀, при котором они сходятся: {100 * hb:.2f}")
    print(f"  ΛCDM Planck: r_d модели {r_sound(_P, Z_DRAG):.2f} Мпк, BAO требуют {need_l:.2f} ({100 * (1 - r_sound(_P, Z_DRAG) / need_l):.1f} % меньше)")
    out["rd_need"] = dict(dsic_model=r_sound(M, Z_DRAG), dsic_need=need_d, H0_bao=100 * hb,
                          lcdm_model=r_sound(_P, Z_DRAG), lcdm_need=need_l)
    R["joint"] = out


# ==============================================================================
#  [7b] ГДЕ НЕ СХОДИТСЯ ЗВУКОВОЙ ГОРИЗОНТ
# ==============================================================================
class DSICAmp(DSIC):
    """Проба: раньше точки включения веществом считается только доля f массы связи ядра в этой точке,
    остаток (1 - f) продолжается в прошлое как постоянная плотность. Темп непрерывен; f = 1 — модель статьи."""

    def __init__(self, mu0, h, f):
        super().__init__(mu0, h)
        self.f = f
        self.Xc = (1 - self.Or) * (1 + self.zc) ** 3 * self.Fc * (1 - f)
        self.Am = f * self.Am                                # вещество — доля f; остаток входит в X(z)

    def F(self, z):
        z = np.asarray(z, float)
        early = self.f * self.Fc + self.Xc / ((1 - self.Or) * (1 + z) ** 3)
        return np.where(z <= self.zc, super().F(z), early)


def omega_m_at(M, z):
    """Плотность вещества в эпоху z в единицах ω = Ω h²: множитель при (1+z)³ в формуле темпа."""
    return float((1 - M.Or) * M.F(z) * M.h ** 2)


def z_equality(M):
    return brentq(lambda z: (1 - M.Or) * (1 + z) ** 3 * M.F(z) - M.Or * (1 + z) ** 4, 100, 1e5)


def block_sound(SN):
    head("[7b] Где не сходится звуковой горизонт")
    S, J = SN["hf"], R["joint"]["n13"]
    out = {}
    # 1. фаза, которую выбирают поздние данные (SN + BAO со свободным r_d), и угол реликта
    mu_l = R["late"]["sn_bao"]["n13"]["mu0"]
    h = solve_h_theta(lambda hh: DSIC(mu_l, hh))
    M = DSIC(mu_l, h)
    rs, rd, need = r_sound(M, Z_STAR), r_sound(M, Z_DRAG), chi2_bao(M, None)[1]
    out["late_phase"] = dict(mu0=mu_l, H0=100 * h, rd=rd, rd_need=need, deficit=100 * (1 - rd / need))
    print(f"  фаза поздних данных μ₀ = {mu_l:.4f}, H₀ из угла реликта {100 * h:.2f}: r_d модели {rd:.2f} Мпк, "
          f"BAO нужен {need:.2f} — на {100 * (1 - rd / need):.2f} % короче")
    # 2. когда набирается звуковой горизонт и где идёт свет
    frac = {f"{Z:.0e}": r_sound(M, Z) / rs for Z in (1e4, 1e5, 1e6)}
    out["rs_frac"] = frac
    print("  доля r*, набранная раньше эпохи z: " + ", ".join(f"z > {k}: {100 * v:.2f} %" for k, v in frac.items()))
    dms = DM(M, Z_STAR)
    path = dict(low=DM(M, 2.33) / dms, high=1 - DM(M, 2.33) / dms)
    out["path"] = path
    print(f"  путь света до реликта: от нас до z = 2.33 — {100 * path['low']:.1f} %, от z = 2.33 до рекомбинации — {100 * path['high']:.1f} %")
    L = LCDM(J["lcdm"]["Om"], J["lcdm"]["H0"] / 100)
    ratio = dict(dsic=DM(M, Z_STAR) / DM(M, 2.33), lcdm=DM(L, Z_STAR) / DM(L, 2.33))
    out["dm_ratio"] = ratio
    lo, hi = R["early_misc"]["need"]
    print(f"  D_M(z*)/D_M(2.33): DSIC {ratio['dsic']:.4f}, ΛCDM (совместный) {ratio['lcdm']:.4f}; нужно {lo:.3f}–{hi:.3f}")
    # 2b. что дал бы сдвиг точки включения: масса связи ядра там минимальна, поэтому любой сдвиг добавляет вещество
    shift = {f"{zs:g}": float(M.Ecore(zs) ** 2 / (1 + zs) ** 3 / M.Fc) for zs in (1.2, 2.0, 3.0, 4.0)}
    out["shift"] = shift
    print(f"  точка включения z_h = {M.zc:.3f}; включение в другой точке дало бы вещества больше на: "
          + ", ".join(f"z = {k}: {100 * (v - 1):+.1f} %" for k, v in shift.items()))
    # 3. что делает свободный w_c
    Mp = DSIC(J["dsic"]["mu0"], J["dsic"]["H0"] / 100)
    Mw = DSIC(J["wc"]["mu0"], J["wc"]["H0"] / 100, J["wc"]["wc"])
    sh = (73.04 - J["wc"]["H0"]) / np.hypot(1.04, J["wc"]["H0_err"])
    out["wc"] = dict(om_p6=omega_m_at(Mp, Z_STAR), om_wc=omega_m_at(Mw, Z_STAR), zeq_p6=z_equality(Mp), zeq_wc=z_equality(Mw),
                     age_wc=age(Mw, 0), sh0es_sigma=sh)
    print(f"  свободный w_c: вещество на рекомбинации ω_m = {out['wc']['om_wc']:.4f} (при w_c = 0: {out['wc']['om_p6']:.4f}), "
          f"равенство вещества и излучения z = {out['wc']['zeq_wc']:.0f} (при w_c = 0: {out['wc']['zeq_p6']:.0f}); "
          f"H₀ от SH0ES (73.04 ± 1.04) на {sh:.1f}σ")
    # 4. проба: меньше вещества, чем масса связи в точке включения
    def parts(Mx):
        r = r_sound(Mx, Z_DRAG)
        return chi2_sn_marg(Mx.E, S), chi2_bao(Mx, r)[0], chi2_theta(Mx), r
    f = lambda p: 1e30 if not (0.6 < p[0] < 0.9 and 0.5 < p[1] < 0.9 and 0.5 < p[2] < 1.2) else sum(parts(DSICAmp(*p))[:3])
    r = minimize(f, [0.762, 0.675, 0.92], method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-6, maxiter=4000))
    r = minimize(f, r.x, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-7, maxiter=4000))
    e = hessian_err(f, r.x, [1e-3, 2e-3, 3e-3])
    Ma = DSICAmp(*r.x)
    p = parts(Ma)
    out["amp"] = dict(mu0=r.x[0], mu0_err=e[0], H0=100 * r.x[1], H0_err=100 * e[1], f=r.x[2], f_err=e[2], sn=p[0], bao=p[1], th=p[2],
                      tot=r.fun, rd=p[3], dchi2=r.fun - J["lcdm"]["tot"], dchi2_f1=J["dsic"]["tot"] - r.fun,
                      om=omega_m_at(Ma, Z_STAR), rest_zh=Ma.Xc / float(Ma.E(Ma.zc)) ** 2)
    a = out["amp"]
    print(f"  проба «вещества меньше на долю 1 − f»: μ₀ = {a['mu0']:.4f} ± {a['mu0_err']:.4f}, H₀ = {a['H0']:.2f} ± {a['H0_err']:.2f}, "
          f"f = {a['f']:.3f} ± {a['f_err']:.3f}; χ²: SN {a['sn']:.2f}, BAO {a['bao']:.2f}, θ* {a['th']:.2f}, итого {a['tot']:.2f}, r_d = {a['rd']:.2f}")
    print(f"    Δχ² к ΛCDM = {a['dchi2']:+.2f}; f = 1 хуже на {a['dchi2_f1']:.2f}; ω_m на рекомбинации = {a['om']:.4f}; "
          f"доля остатка в темпе при z_h: {100 * a['rest_zh']:.1f} %")
    R["sound"] = out


# ==============================================================================
#  [8] РОСТ СТРУКТУР fσ8
# ==============================================================================
def block_growth():
    head("[8] Рост структур: fσ8, набор Gold-2017 (18 точек), линейный рост в пределе ОТО")
    Cm = np.diag(FS8[:, 2] ** 2)
    for a, i in enumerate(FS8_WIGGLEZ_IDX):
        for b_, j in enumerate(FS8_WIGGLEZ_IDX):
            Cm[i, j] = FS8_WIGGLEZ_COV[a, b_]
    Ci = np.linalg.inv(Cm)
    h = solve_h_theta(lambda hh: DSIC(0.76, hh))
    out = {}
    for name, key, M, Am in (("DSIC (μ₀ = 0.76)", "dsic", DSIC(0.76, h), DSIC(0.76, h).Am), ("ΛCDM Planck", "lcdm", _P, 0.3153)):
        sol = growth_factor(M, Am, 100.0)
        x = -np.log(1 + FS8[:, 0])
        D, Dp = sol.sol(x)
        D0 = sol.sol(0.0)[0]
        t = (Dp / D) * (D / D0)                         # f·D/D0, умножается на σ8
        s8 = (t @ Ci @ FS8[:, 1]) / (t @ Ci @ t)
        r = s8 * t - FS8[:, 1]
        c2, err = float(r @ Ci @ r), float(1 / np.sqrt(t @ Ci @ t))
        print(f"  {name}: σ₈ = {s8:.3f} ± {err:.3f}, χ² = {c2:.1f}")
        out[key] = dict(sigma8=float(s8), err=err, chi2=c2)
    R["growth"] = out


# ==============================================================================
#  [9] РАННИЕ ЧЁРНЫЕ ДЫРЫ: БЮДЖЕТ ВРЕМЕНИ
# ==============================================================================
def block_bh():
    head("[9] Ранние чёрные дыры: хватает ли времени на рост в пределе Эддингтона")
    h = solve_h_theta(lambda hh: DSIC(0.76, hh))
    M = DSIC(0.76, h)
    out = {}
    for seed in (1e2, 1e4):
        print(f"  зародыш {seed:.0e} масс Солнца (время отсчитывается от нуля пространства):")
        okd = okl = 0
        for name, z, mbh in BH_OBJECTS:
            need = np.log(mbh / seed) * T_SALPETER
            td, tl = age(M, z), age(_P, z)
            okd += td >= need
            okl += tl >= need
            print(f"    {name:12s} z = {z:6.3f}  M = {mbh:.1e}  нужно {need:.3f}  есть: DSIC {td:.3f} {'да ' if td >= need else 'нет'}  "
                  f"ΛCDM {tl:.3f} {'да' if tl >= need else 'нет'}")
        print(f"    успевают: DSIC {okd} из 6, ΛCDM {okl} из 6")
        out[f"seed{int(seed)}"] = dict(dsic=int(okd), lcdm=int(okl))
    R["bh"] = out


# ==============================================================================
#  [10] ПОЛНЫЙ СПЕКТР РЕЛИКТА: CAMB + plik-lite (по флагу --cmb)
# ==============================================================================
def block_cmb(SN):
    head("[10] Полный спектр реликта (таблица 7): CAMB + сжатое правдоподобие Planck 2018 (plik-lite)")
    import subprocess
    try:
        import camb
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "camb"])
        import camb
    from camb.dark_energy import DarkEnergyPPF
    pl_dir = os.path.join(DATA_DIR, "planck-lite-py")
    if not os.path.exists(pl_dir):
        subprocess.check_call(["git", "clone", "--depth", "1", "https://github.com/heatherprince/planck-lite-py.git", pl_dir])
    sys.path.insert(0, pl_dir)
    from planck_lite_py import PlanckLitePy
    like = PlanckLitePy(data_directory=os.path.join(pl_dir, "data"), year=2018, spectra="TTTEEE", use_low_ell_bins=True)
    TAU, LMAX = 0.0544, 2600
    S = SN["hf"]

    def camb_params(h, ombh2, omch2, ns, de=None):
        p = camb.CAMBparams()
        if de is not None:
            p.DarkEnergy = de
        p.set_cosmology(H0=100 * h, ombh2=ombh2, omch2=omch2, mnu=0.0, nnu=N_EFF, num_massive_neutrinos=0,
                        omk=0.0, tau=TAU, TCMB=2.7255)
        p.InitPower.set_params(As=2.1e-9, ns=ns)
        p.set_for_lmax(LMAX, lens_potential_accuracy=1)
        return p

    def dsic_de(M):
        """Пространство X(z), которое прибавляет ядро, как гладкая жидкость с таблицей w(a) (PPF)."""
        zc = M.zc
        d = 1e-5
        dw = lambda z: -1 + (np.log(M.X(z + d)) - np.log(M.X(z - d))) / (np.log1p(z + d) - np.log1p(z - d)) / 3
        if getattr(M, "f", 1.0) < 1.0:
            # проба f: за точкой включения остаток X постоянен, w = -1. В точке включения w скачет с 0 на -1;
            # скачок заменён гладкой ступенькой шириной 0.02 по ln(1+z), иначе сплайн CAMB «звенит».
            lo, hi = np.log1p(zc) - 0.02, np.log1p(zc) + 0.02
            z1 = np.linspace(0, np.expm1(lo), 400)
            w1 = dw(z1)
            t = np.linspace(0, 1, 12)[1:-1]
            zr, wr = np.expm1(lo + (hi - lo) * t), w1[-1] + (-1 - w1[-1]) * (3 * t ** 2 - 2 * t ** 3)
            z3 = np.concatenate([[np.expm1(hi)], np.logspace(np.log10(np.expm1(hi) + 1) + 0.01, 5, 60) - 1])
            z = np.concatenate([z1, zr, z3])
            w = np.concatenate([w1, wr, np.full(len(z3), -1.0)])
            a = 1 / (1 + z[::-1])
            w = np.maximum(w[::-1], -60.0)
            a = np.concatenate([[1e-6], a])
            w = np.concatenate([[-1.0], w])
        else:
            z = np.concatenate([np.linspace(0, 0.98 * zc, 400), [0.985 * zc, 0.99 * zc]])
            w = dw(z)
            a = 1 / (1 + z[::-1])
            w = np.maximum(w[::-1], -60.0)
            a = np.concatenate([[1e-4, a[0] * 0.9], a])       # раньше включения плотность X гаснет
            w = np.concatenate([[w[0], w[0]], w])
        de = DarkEnergyPPF()
        de.set_w_a_table(a, w)
        return de

    def spectra(p):
        res = camb.get_results(p)
        cl = res.get_cmb_power_spectra(p, CMB_unit="muK", raw_cl=False)["total"]
        return res, cl[2:, 0], cl[2:, 3], cl[2:, 1]           # TT, TE, EE (D_l), начиная с l = 2

    def chi2_cmb(tt, te, ee):
        """Общая амплитуда находится аналитически: спектр ∝ α."""
        best = minimize_scalar(lambda al: -2 * like.loglike(al * tt, al * te, al * ee, ellmin=2),
                               bounds=(0.8, 1.2), method="bounded", options=dict(xatol=1e-5))
        return best.fun, best.x

    def h_from_theta(make_params, th100):
        """h, при котором 100θ* (по CAMB) равен заданному."""
        f = lambda h: camb.get_background(make_params(h)).get_derived_params()["thetastar"] - th100
        return brentq(f, 0.55, 0.85, xtol=1e-6)

    def eval_lcdm(x):
        ombh2, omch2, th100, ns = x
        mk = lambda h: camb_params(h, ombh2, omch2, ns)
        h = h_from_theta(mk, th100)
        res, tt, te, ee = spectra(mk(h))
        c2, al = chi2_cmb(tt, te, ee)
        return c2, dict(h=h, rd=res.get_derived_params()["rdrag"], amp=al, model=LCDM((ombh2 + omch2) / h ** 2, h),
                        Om=(ombh2 + omch2) / h ** 2)

    def eval_dsic(x):
        mu0, th100, ombh2, ns = x

        def mk(h):
            M = DSIC(mu0, h)
            return camb_params(h, ombh2, M.Am * h * h - ombh2, ns, dsic_de(M))
        h = h_from_theta(mk, th100)
        M = DSIC(mu0, h)
        res, tt, te, ee = spectra(mk(h))
        c2, al = chi2_cmb(tt, te, ee)
        return c2, dict(h=h, rd=res.get_derived_params()["rdrag"], amp=al, model=M, Om=M.Am)

    def eval_dsicf(x):
        mu0, fr, th100, ombh2, ns = x
        if not (0.6 < mu0 < 0.9 and 0.5 < fr < 0.999):
            raise ValueError("вне области")

        def mk(h):
            M = DSICAmp(mu0, h, fr)
            return camb_params(h, ombh2, M.Am * h * h - ombh2, ns, dsic_de(M))
        h = h_from_theta(mk, th100)
        M = DSICAmp(mu0, h, fr)
        res, tt, te, ee = spectra(mk(h))
        c2, al = chi2_cmb(tt, te, ee)
        return c2, dict(h=h, rd=res.get_derived_params()["rdrag"], amp=al, model=M, Om=M.Am)

    def late(info):
        M, rd = info["model"], info["rd"]
        return chi2_sn_marg(M.E, S), chi2_bao(M, rd, True)[0]

    out = {}
    # проверка перевода: темп CAMB против формулы (27)
    M = DSIC(0.744, 0.705)
    bg = camb.get_background(camb_params(0.705, WB, M.Am * 0.705 ** 2 - WB, 0.965, dsic_de(M)))
    zz = np.concatenate([np.linspace(0, 3, 60), np.logspace(0.5, 6, 60)])
    dev = np.max(np.abs(bg.hubble_parameter(zz) / (70.5 * M.E(zz)) - 1))
    print(f"  перевод фона в CAMB: макс. отклонение H(z) от формулы (27) на 0 ≤ z ≤ 10⁶: {dev:.1e}")
    out["camb_dev"] = float(dev)
    Mf = DSICAmp(0.762, 0.674, 0.92)
    bg = camb.get_background(camb_params(0.674, WB, Mf.Am * 0.674 ** 2 - WB, 0.965, dsic_de(Mf)))
    devf = np.max(np.abs(bg.hubble_parameter(zz) / (67.4 * Mf.E(zz)) - 1))
    print(f"  то же для пробы f = 0.92: {devf:.1e}")
    out["camb_dev_f"] = float(devf)
    names = {"ΛCDM": ("ω_b", "ω_c", "100θ*", "n_s"), "DSIC": ("μ₀", "100θ*", "ω_b", "n_s"),
             "DSIC, проба f": ("μ₀", "f", "100θ*", "ω_b", "n_s")}
    amp = R.get("sound", {}).get("amp", dict(mu0=0.7626, f=0.922))
    start = {"ΛCDM": [0.02237, 0.1200, 1.0396, 0.965], "DSIC": [0.744, 1.0396, 0.02237, 0.965],
             "DSIC, проба f": [amp["mu0"], amp["f"], 1.0411, 0.02240, 0.966]}
    step = {"ΛCDM": [2e-4, 2e-3, 5e-4, 6e-3], "DSIC": [4e-3, 5e-4, 2e-4, 6e-3], "DSIC, проба f": [3e-3, 0.02, 4e-4, 1.5e-4, 5e-3]}
    hstep = {"ΛCDM": [1.5e-4, 1.2e-3, 3e-4, 4e-3], "DSIC": [2.5e-3, 3e-4, 1.5e-4, 4e-3], "DSIC, проба f": [2e-3, 0.01, 3e-4, 1e-4, 3e-3]}
    evals = {"ΛCDM": eval_lcdm, "DSIC": eval_dsic, "DSIC, проба f": eval_dsicf}
    cases = [("только реликт", "ΛCDM"), ("только реликт", "DSIC"), ("реликт + SN + BAO", "ΛCDM"),
             ("реликт + SN + BAO", "DSIC"), ("реликт + SN + BAO", "DSIC, проба f")]
    only = os.environ.get("DSIC_CMB_CASES")                  # например "0" — только первый из четырёх подборов
    warm_path = os.environ.get("DSIC_CMB_START")             # необязательно: JSON с начальными точками подборов
    warm = json.load(open(warm_path)) if warm_path and os.path.exists(warm_path) else {}
    if RUN_CMB_F and not RUN_CMB and not only:
        only = "4"                                           # --cmb-f: только проба f
    for idx, (mode, name) in enumerate(cases):
        if True:
            if only and str(idx) not in only.split(","):
                continue
            t0, ev = time.time(), evals[name]
            cache, best = {}, [np.inf]

            def f(x):
                key = tuple(np.round(x, 10))
                if key not in cache:
                    try:
                        c2, info = ev(x)
                        cache[key] = c2 if mode == "только реликт" else c2 + sum(late(info))
                    except Exception:
                        cache[key] = 1e30
                    best[0] = min(best[0], cache[key])
                    if len(cache) % 25 == 0:
                        print(f"      … {name}: вызовов {len(cache)}, лучший χ² {best[0]:.2f}, прошло {time.time() - t0:.0f} с")
                        sys.stdout.flush()
                return cache[key]
            jkey = {"ΛCDM": "lcdm", "DSIC": "dsic", "DSIC, проба f": "dsicf"}[name] + ("_cmb" if mode == "только реликт" else "_joint")
            x0, scale, r = warm.get(jkey, start[name]), (0.35 if jkey in warm else 1.0), None
            n = len(x0)
            for _ in range(4):                               # симплекс с перезапусками, пока минимум не перестанет улучшаться
                simplex = np.array([x0] + [np.array(x0) + np.eye(n)[i] * step[name][i] * scale for i in range(n)])
                rr = minimize(f, x0, method="Nelder-Mead",
                              options=dict(initial_simplex=simplex, xatol=2e-6, fatol=0.005, maxiter=int(os.environ.get("DSIC_CMB_MAXITER", 400))))
                gain = (r.fun - rr.fun) if r is not None else np.inf
                if r is None or rr.fun < r.fun:
                    r = rr
                if gain < 0.02:
                    break
                x0, scale = list(r.x), scale * 0.4
            c2, info = ev(r.x)
            sn, bao = late(info)
            err = hessian_err(f, r.x, hstep[name])
            # ошибка H0: линейное распространение через градиент h(x)
            g = np.array([(ev(r.x + np.eye(n)[i] * hstep[name][i])[1]["h"] - ev(r.x - np.eye(n)[i] * hstep[name][i])[1]["h"])
                          / (2 * hstep[name][i]) for i in range(n)])
            Hm = np.zeros((n, n))
            for i in range(n):
                for j in range(n):
                    ei, ej = np.eye(n)[i] * hstep[name][i], np.eye(n)[j] * hstep[name][j]
                    Hm[i, j] = (f(r.x + ei + ej) - f(r.x + ei - ej) - f(r.x - ei + ej) + f(r.x - ei - ej)) / (4 * hstep[name][i] * hstep[name][j])
            cov = np.linalg.inv(0.5 * Hm)
            eH = 100 * float(np.sqrt(g @ cov @ g))
            key = f"{name} | {mode}"
            out[jkey] = dict(names=list(names[name]), x=[float(v) for v in r.x], err=[float(v) for v in err],
                            chi2_cmb=float(c2), chi2_sn=float(sn), chi2_bao=float(bao), total=float(r.fun),
                            H0=100 * info["h"], H0_err=eH, rd=float(info["rd"]), amp=float(info["amp"]),
                            Om_or_Am=float(info["Om"]), nfev=len(cache))
            pars = ", ".join(f"{nm} = {v:.5f} ± {e:.5f}" for nm, v, e in zip(names[name], r.x, err))
            print(f"  {key}: χ²(реликт) = {c2:.2f}, χ²(SN) = {sn:.2f}, χ²(BAO) = {bao:.2f}; минимум подбора {r.fun:.2f}")
            print(f"      {pars}; H₀ = {100 * info['h']:.2f} ± {eH:.2f}; Ωₘ или A_m = {info['Om']:.4f}; r_d = {info['rd']:.2f} Мпк; "
                  f"вызовов {len(cache)}, {time.time() - t0:.0f} с")
            sys.stdout.flush()
            if mode != "только реликт":
                _, _, pulls = chi2_bao(info["model"], info["rd"], True, True)
                txt = ", ".join(f"{DESI_TRACER[round(float(z), 3)]} {k}: {p:+.1f}" for z, k, p in zip(DESI_Z, DESI_KIND, pulls))
                print(f"      отклонения точек BAO, в σ: {txt}")
                out[jkey]["pulls"] = [float(p) for p in pulls]
    if "dsicf_joint" in out:
        # сравнение пробы f с ΛCDM и DSIC (f = 1); если они в этом запуске не считались, берутся из полного прогона с флагом --cmb (табл. 7): 16.11.2026
        ref = {"lcdm_joint": dict(total=1988.8608, chi2_cmb=587.7170, chi2_sn=1389.7777, chi2_bao=11.3661),
               "dsic_joint": dict(total=1997.6743, chi2_cmb=590.5316, chi2_sn=1387.0500, chi2_bao=20.0927)}
        a = out["dsicf_joint"]
        L, D = out.get("lcdm_joint", ref["lcdm_joint"]), out.get("dsic_joint", ref["dsic_joint"])
        om = a["Om_or_Am"] * (a["H0"] / 100) ** 2
        a.update(dchi2_lcdm=a["total"] - L["total"], dchi2_f1=D["total"] - a["total"], om_rec=om)
        print(f"\n  проба f против ΛCDM: Δχ² итого = {a['dchi2_lcdm']:+.2f} (реликт {a['chi2_cmb'] - L['chi2_cmb']:+.2f}, "
              f"SN {a['chi2_sn'] - L['chi2_sn']:+.2f}, BAO {a['chi2_bao'] - L['chi2_bao']:+.2f}); с поправкой на лишний параметр {a['dchi2_lcdm'] + 2:+.2f}")
        print(f"  проба f против DSIC с f = 1: лучше на {a['dchi2_f1']:.2f}; вещество на рекомбинации ω_m = f·A_m·h² = {om:.4f}")
    R["cmb"] = out


# ==============================================================================
#  ЗАПУСК
# ==============================================================================
def main():
    print("dsic3_test.py — проверка чисел статьи Dual Scale Inversion Cosmology")
    print(f"ω_r = {W_R:.4e}; поправка приближения r_s (калибровка на ΛCDM Planck): {CAL:.6f}")
    block_core()
    M = block_early()
    block_de(M)
    load_desi_check()
    SN = load_sn()
    block_late(SN)
    block_joint(SN)
    block_sound(SN)
    block_growth()
    block_bh()
    if RUN_CMB or RUN_CMB_F:
        block_cmb(SN)
    head("ИТОГ (JSON для сверки со статьёй)")

    def clean(o):
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, (np.floating, np.integer)):
            return float(o)
        if isinstance(o, np.ndarray):
            return clean(o.tolist())
        return o
    js = json.dumps(clean(R), ensure_ascii=False)
    print(js)
    with open(os.path.join(DATA_DIR, "dsic3_test_results.json"), "w") as fh:
        fh.write(js)
    print(f"\nГотово за {time.time() - T_START:.0f} с.")


if __name__ == "__main__":
    main()
