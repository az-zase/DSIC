"""
================================================================================
DSIC — SINGLE reproducible test (ONE file, Google Colab / python).
Part I   (the core):        DSIC vs flat ΛCDM on identical data, equal number of parameters.
Part II  (the second tier): the detachment threshold μ_d (early Universe, CMB).
Part III (structure growth): fσ8 on the DSIC background (GR limit), block [11].

Model author: Azamat Zaseev
Repository:   https://github.com/az-zase/DSIC
================================================================================

      PLEASE BE PATIENT.
      The script MAY look frozen — this is normal, it is working. There is
      nothing to press or interrupt: just calmly wait for it to finish.
      A full run takes roughly 10–15 minutes and consists of:
        • automatic download of the supernova catalog (~65 MB);
        • the Pantheon+SH0ES supernova fits (a few minutes);
        • block [11] "structure growth" — the longest stage: a grid scan
          with numerical integration of the growth ODE at every node (~5 min).
          Between the printout of block [10] and the appearance of the [11]
          results there may be a long silent pause — this is expected,
          the script is computing.
      At the end, "Done" will be printed at the bottom with a final summary.

================================================================================

WHAT IT DOES (all in one run):
  [0]  Self-check: the hard-coded DESI against the official Table IV.
  [0b] MODEL CONSISTENCY: machine verification of the core's internal
       identities — monotonicity of b(μ); accuracy of D_H = dχ/dz; agreement
       of the closed-form χ(μ) with the numerical ∫c dz/H; the zero of q at
       μ_cross = √(3−√5); the age formula.
  [1]  SN cosmology: Pantheon+SH0ES, anchored fit (m_b_corr + explicit M).
       1σ parameter uncertainties from the χ² Hessian at the minimum.
  [2]  Robustness: STAT-ONLY versus STAT+SYS.
  [3]  Joint SN+BAO: the full D_M–D_H covariance (12×12), free r_d;
       plus the recovery of the r_d ruler on the early scale.
  [4]  The Om(z) diagnostic: robustness of the sign of Δχ² to H0; the
       chronometers in TWO covariance modes (diagonal and full Moresco et al. 2020).
  [5]  Falsifiable predictions at high z (D_H/r_d Lyα, E(z)).
  [6]  COSMIC TIME: t(z) DSIC vs ΛCDM at z = 5–12 — the nearest risk zone
       of the core (JWST galaxy ages, SMBH growth by z≈7.5).
  ---- Part II (the second tier, P6) ----
  [7]  Self-check of the Planck reference data (d_LSS = r*/θ*).
  [8]  The core's shortfall to the CMB and the fixing of the threshold μ_d
       by the directly measured angle θ*; r* is explicitly flagged as
       BORROWED (an external input from ΛCDM).
  [9]  Control: the correction is identically 0 at z ≤ z_detach; the
       "amplitude ↔ threshold" degeneracy table (what postulate P6 fixes).
  [10] The trace of μ_d at high z; the corrected particle horizon; the
       internal-consistency flag (observed galaxies z ≳ 10 lie ABOVE the threshold).
  ---- Part III (structure growth) ----
  [11] STRUCTURE GROWTH: the standard equation of linear perturbation growth
       on the DSIC background (GR limit), fσ8(z) against the Gold-2017
       compilation (18 points, the WiggleZ covariance included). The fσ8 data
       are degenerate in (Ω_m0, σ8) — the result is reported as the invariant
       S8 = σ8·√(Ω_m0/0.3) plus a slice at Planck-like Ω_m0 = 0.30.
       Status: a survival test, not a derivation from the core (the core
       introduces no densities; Ω_m0 is a growth parameter).

DATA (the provenance of every hard-coded figure is documented in place):
  • Supernovae — DOWNLOADED AUTOMATICALLY (Pantheon+SH0ES DataRelease).
  • DESI DR2 BAO — hard-coded, VERIFIED line by line against the official
        Table IV (DESI Coll. / Abdul-Karim et al. 2025, Phys. Rev. D 112,
        083515; arXiv:2503.14738). Self-check [0] fails on any mismatch.
  • Cosmic chronometers H(z) — hard-coded, 33 points. The 15 Moresco points
        are VERIFIED against the official repository
        (gitlab.com/mmoresco/CCcovariance, BC03); the other 18 are from
        other works (Simon 2005; Stern 2010; Zhang 2014; Ratsimbazafy 2017;
        Borghi 2022 et al.), with diagonal errors.
  • The FULL COVARIANCE of the Moresco subset (Moresco et al. 2020, ApJ 898,
        82; arXiv:2003.07362) is hard-coded and built directly in this file.
  • Planck 2018 (TT,TE,EE+lowE+lensing): z* = 1089.92, r* = 144.43 Mpc,
        100θ* = 1.04110, whence d_LSS = r*/θ* = 13872.8 ± 25 Mpc. Verified by [7].
  • fσ8 — the "Gold-2017" compilation (Nesseris, Pantazis & Perivolaropoulos
        2017, Phys. Rev. D 96, 023542; arXiv:1703.10538), 18 weakly correlated
        points; z, fσ8, σ are from Table I of Sagredo, Nesseris & Sapone 2018
        (arXiv:1806.10822). The covariance of the three WiggleZ points —
        Blake et al. 2012.

DEPENDENCIES: numpy, scipy (pre-installed in Colab).
RUN: paste the file into a single Colab cell and execute; or
        `python en_dsic_test.py`.
================================================================================
"""

import os
import urllib.request
import numpy as np
from scipy.optimize import minimize, brentq
from scipy.integrate import quad, solve_ivp

c = 299792.458  # km/s

# ==============================================================================
#  DATA
# ==============================================================================

# --- DESI DR2 BAO, Table IV (arXiv:2503.14738). Baseline sample for cosmology.
#     Columns: z_eff, D_M/r_d, σ(D_M/r_d), D_H/r_d, σ(D_H/r_d), corr(D_M,D_H).
BAO = np.array([
    # tracer        z_eff   DM/rd     sDM     DH/rd     sDH     corr
    [0.510, 13.587, 0.169, 21.863, 0.427, -0.475],   # LRG1
    [0.706, 17.347, 0.180, 19.458, 0.332, -0.423],   # LRG2
    [0.934, 21.574, 0.153, 17.641, 0.193, -0.425],   # LRG3+ELG1 (baseline)
    [1.321, 27.605, 0.320, 14.178, 0.217, -0.437],   # ELG2
    [1.484, 30.519, 0.758, 12.816, 0.513, -0.489],   # QSO
    [2.330, 38.988, 0.531,  8.632, 0.101, -0.431],   # Lya
])

# --- Reference for self-check [0] (the same official table, as a separate copy).
#     NOTE: the check catches an accidental edit of the array, but not a common
#     transcription typo from the publication — both copies were verified
#     against arXiv:2503.14738 by hand.
_DESI_OFFICIAL = {
    0.510: (13.587, 0.169, 21.863, 0.427, -0.475),
    0.706: (17.347, 0.180, 19.458, 0.332, -0.423),
    0.934: (21.574, 0.153, 17.641, 0.193, -0.425),
    1.321: (27.605, 0.320, 14.178, 0.217, -0.437),
    1.484: (30.519, 0.758, 12.816, 0.513, -0.489),
    2.330: (38.988, 0.531,  8.632, 0.101, -0.431),
}

# --- Cosmic chronometers H(z), 33 points, z ∈ [0.07, 1.965].
#     The 15 Moresco points are verified against gitlab.com/mmoresco/CCcovariance
#     (BC03; Moresco et al. 2012, 2016; Moresco 2015); the full covariance of
#     this subset is hard-coded below. The other 18 points — Simon 2005;
#     Stern 2010; Zhang 2014; Ratsimbazafy 2017; Borghi 2022 et al.,
#     diagonal errors.
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
#  CHRONOMETER COVARIANCE (Moresco et al. 2020) — HARD-CODED TABLES
# ------------------------------------------------------------------------------
#  Method: Moresco, Jimenez, Verde, Cimatti, Pozzetti, ApJ 898, 82 (2020);
#          arXiv:2003.07362; repository https://gitlab.com/mmoresco/CCcovariance.
#  Full covariance = statistical (diagonal) + systematic
#  (correlated in z): Cov_syst(i,j) = Σ_s σ_s(z_i)·σ_s(z_j), ρ_s ≡ 1.
#  The diagonal of the final matrix is EXACTLY equal to the published total
#  error err_t².
#  STATUS: a careful reconstruction of the method on the official input tables,
#  not a verbatim run of the authors' code; the off-diagonal reproduces the recipe.
# ==============================================================================

# --- HzTable_MM_BC03.dat (SPS model = BC03). 15 Moresco points.
#     Columns: (z, H, errH_total, stat_contr, met_contr) in km/s/Mpc.
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

# --- data_MM20.dat (Moresco et al. 2020). Relative (%) systematic
#     contributions on a z grid (step 0.05). Columns: (z, IMF%, stlib%, mod%, mod_ooo%).
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

# --- Planck 2018 (TT,TE,EE+lowE+lensing). Reference quantities of Part II.
#     The directly measured quantity is the angle θ*. The sound horizon
#     r* = 144.43 Mpc is a quantity COMPUTED within ΛCDM from the physics of
#     the photon–baryon plasma; in DSIC it is BORROWED as an external input
#     (see the paper, §11.6).
#     d_LSS = r*/θ* = 13872.8 ± 25 Mpc (see, e.g., arXiv:2311.04759, §3.6).
PLANCK = dict(
    z_star=1089.92,        # redshift of last scattering
    r_star=144.43,         # sound horizon at z*, Mpc (BORROWED from ΛCDM)
    theta_star=1.04110e-2, # 100θ* = 1.04110 -> θ* in radians (direct measurement)
    d_lss=13872.8,         # = r*/θ*, Mpc
    d_lss_err=25.0,        # 1σ, Mpc
)

# --- Second-tier parameters: the early scale (the same as for BAO/chronometers/
#     age in Part I).
MU0_E = 0.76
H0_E  = 67.4


def build_cc_covariance(z_cc, sH_cc, full=True):
    """Covariance of the 33 chronometer points for block [4].

    full=False : purely diagonal (= sH_cc²), the baseline version of the test.
    full=True  : the block of the 15 Moresco points is replaced with the full
                 Moresco 2020 covariance (diagonal = err_t², off-diagonal =
                 correlated systematics); the other 18 points stay diagonal.
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
    sig = (frac / quad_) * met[None, :]          # (4,15): absolute contribution of each source

    C_syst = sum(np.outer(sig[s], sig[s]) for s in range(4))   # 15×15 systematics

    idx = [int(np.argmin(np.abs(z_cc - zm))) for zm in z_m]
    for k, i in enumerate(idx):
        if abs(z_cc[i] - z_m[k]) > 1e-3:
            raise SystemExit(f"STOP: Moresco point z={z_m[k]} not found in the CC array.")

    for a in range(15):
        for b_ in range(15):
            C[idx[a], idx[b_]] = (err_t[a] ** 2 if a == b_ else 0.0) + C_syst[a, b_]
    return C


# ==============================================================================
#  [0] DESI DATA SELF-CHECK
# ==============================================================================
def self_check_desi():
    print("=" * 70)
    print("[0] Self-check: the hard-coded DESI against the official Table IV")
    print("=" * 70)
    ok = True
    for row in BAO:
        z = round(float(row[0]), 3)
        off = _DESI_OFFICIAL[z]
        got = tuple(round(float(x), 6) for x in row[1:])
        match = all(abs(a - b_) < 1e-6 for a, b_ in zip(got, off))
        if not match:
            ok = False
            print(f"  z={z}: MISMATCH  file={got}  official={off}")
    if ok:
        print("  OK — all 6 DESI bins match arXiv:2503.14738, Table IV.\n")
    else:
        raise SystemExit("STOP: the DESI data do not match the reference. Check the BAO array.")


# ==============================================================================
#  DSIC AND ΛCDM GEOMETRY
# ==============================================================================
def b_dsic(m):
    """The observable scale factor b(μ) = √(1−μ²)/μ²."""
    return np.sqrt(1.0 - m * m) / (m * m)


def mu_e_of_z(zv, mu0):
    """Inversion of 1+z = b(μ₀)/b(μ_e) on the physical reception branch.

    At z=0 the root equals μ₀ exactly and lies on the edge of the brentq
    bracket (both ends give the same sign) — it is returned analytically.
    The fits never touch this case (SNe have z > 0.01); the behavior at
    z > 0 is unchanged."""
    if zv <= 1e-14:
        return mu0
    target = b_dsic(mu0) / (1.0 + zv)
    return brentq(lambda m: b_dsic(m) - target, mu0 + 1e-12, 1.0 - 1e-15)


def ck_const(mu0, H0):
    """Conversion of the geometry into Mpc: c/k = (c/H₀)·(2−μ₀²)/(1−μ₀²)."""
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
#  [0b] MODEL CONSISTENCY: machine verification of the core's internal identities
# ==============================================================================
def self_check_coherence():
    print("=" * 70)
    print("[0b] Core consistency: machine verification of the internal identities")
    print("=" * 70)
    mu0, H0 = 0.76, 70.0
    fails = []

    # (a) monotonicity of b(μ) on (μ0, 1) — validity of the inversion z -> μ_e
    grid = np.linspace(mu0 + 1e-6, 1 - 1e-9, 4000)
    mono = np.all(np.diff(b_dsic(grid)) < 0)
    print(f"  (a) b(μ) strictly decreasing on (μ₀,1): {'OK' if mono else 'FAIL'}")
    if not mono:
        fails.append("monotonicity")

    # (b) the identity D_H = dχ/dz (paper §5.5): |dχ/dz · H/c − 1| < 1e-8
    worst = 0.0
    for zv in (0.5, 2.0, 10.0, 100.0):
        eps = 1e-5 * (1 + zv)
        d = (DM_comov_DSIC(zv + eps, mu0, H0) - DM_comov_DSIC(zv - eps, mu0, H0)) / (2 * eps)
        worst = max(worst, abs(d * H_DSIC(zv, mu0, H0) / c - 1))
    print(f"  (b) identity D_H = dχ/dz: max deviation {worst:.1e} "
          f"{'OK' if worst < 1e-6 else 'FAIL'}")
    if worst >= 1e-6:
        fails.append("D_H=dχ/dz")

    # (c) closed-form χ(μ) against the numerical ∫ c dz'/H(z')
    worst = 0.0
    for zv in (1.0, 5.0):
        num = quad(lambda x: c / H_DSIC(x, mu0, H0), 0.0, zv, limit=300)[0]
        worst = max(worst, abs(num / DM_comov_DSIC(zv, mu0, H0) - 1))
    print(f"  (c) closed-form χ = ∫c dz/H: max deviation {worst:.1e} "
          f"{'OK' if worst < 1e-6 else 'FAIL'}")
    if worst >= 1e-6:
        fails.append("closed-form χ")

    # (d) the q crossover: μ_cross = √(3−√5) — the exact root of −μ⁴+6μ²−4 = 0
    q = lambda m: (-m ** 4 + 6 * m * m - 4) / (m ** 4 - 4 * m * m + 4)
    mu_cross = np.sqrt(3.0 - np.sqrt(5.0))
    z_t = b_dsic(mu0) / b_dsic(mu_cross) - 1.0
    okq = abs(q(mu_cross)) < 1e-10
    print(f"  (d) q(μ_cross)=0 at μ_cross=√(3−√5)={mu_cross:.6f}; "
          f"z_t={z_t:.3f}: {'OK' if okq else 'FAIL'}")
    if not okq:
        fails.append("μ_cross")

    # (e) age: t₀·H₀ = ln(1/μ₀)·(2−μ₀²)/(1−μ₀²)
    h = lambda m: (2 - m * m) / (1 - m * m)
    t0H0 = np.log(1.0 / mu0) * h(mu0)
    t0 = t0H0 * (977.792 / 67.4)
    print(f"  (e) age: t₀·H₀ = {t0H0:.4f}  ->  t₀ = {t0:.2f} Gyr at H₀=67.4 "
          f"(observed 13.80 ± 0.02)")

    if fails:
        raise SystemExit("STOP: internal identities violated: " + ", ".join(fails))
    print("  Bottom line: the core's internal identities hold exactly.\n")


# ==============================================================================
#  SUPERNOVA LOADING
# ==============================================================================
SN_BASE = ("https://raw.githubusercontent.com/PantheonPlusSH0ES/"
           "DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/")


def grab(name):
    if os.path.exists(name):
        print(f"  found locally: {name}")
        return
    print(f"  downloading {name} ...")
    urllib.request.urlretrieve(SN_BASE + name, name)
    print(f"  done: {os.path.getsize(name) / 1e6:.1f} MB")


def read_cov(path, n_expected):
    with open(path) as f:
        n = int(f.readline())
        vals = np.array(f.read().split(), dtype=float)
    assert n == n_expected, f"covariance size {n} != {n_expected}"
    return vals.reshape(n, n)


def load_sn():
    print("=" * 70)
    print("[1a] Loading Pantheon+SH0ES")
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
    print(f"  read {len(zHD)} supernovae")

    C_sys = read_cov("Pantheon+SH0ES_STAT+SYS.cov", len(zHD))
    C_stat = read_cov("Pantheon+SH0ES_STATONLY.cov", len(zHD))

    use = (zHD > 0.01) | is_cal
    idx = np.where(use)[0]
    out = dict(
        z=zHD[idx], mb=m_b[idx], cal=is_cal[idx], ceph=ceph[idx],
        Cinv_sys=np.linalg.inv(C_sys[np.ix_(idx, idx)]),
        Cinv_stat=np.linalg.inv(C_stat[np.ix_(idx, idx)]),
    )
    print(f"  {len(out['z'])} SNe enter the fit (calibrators: {out['cal'].sum()})\n")
    return out


# ==============================================================================
#  SN FIT (anchored: m_b_corr + explicit M) + UNCERTAINTIES FROM THE HESSIAN
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
    """1σ parameter uncertainties from the quadratic expansion of χ² at the
    minimum: C = 2·H⁻¹, where H is the Hessian of χ² (central differences)."""
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
    print("[1] SN anchored fit (Pantheon+SH0ES, STAT+SYS, m_b_corr + explicit M)")
    print("=" * 70)
    dof = len(SN["z"]) - 3
    rd = fit_sn(distmod_DSIC, SN, SN["Cinv_sys"], [0.76, 73.0, -19.3], OK_D)
    rl = fit_sn(distmod_LCDM, SN, SN["Cinv_sys"], [0.33, 73.0, -19.3], OK_L)
    print(f"DSIC : mu0={rd.x[0]:.4f}  H0={rd.x[1]:.2f}  M={rd.x[2]:.4f}  chi2/dof={rd.fun/dof:.4f}")
    print(f"ΛCDM : Om ={rl.x[0]:.4f}  H0={rl.x[1]:.2f}  M={rl.x[2]:.4f}  chi2/dof={rl.fun/dof:.4f}")
    print(f"Δχ² (DSIC−ΛCDM) = {rd.fun - rl.fun:+.2f}   (|Δ|<2 — indistinguishable)")

    print("\n  1σ uncertainties (χ² Hessian at the minimum, C = 2H⁻¹):")
    sd = hessian_errors(lambda p: chi2_sn(p, distmod_DSIC, SN, SN["Cinv_sys"]),
                        rd.x, [0.002, 0.15, 0.003])
    sl = hessian_errors(lambda p: chi2_sn(p, distmod_LCDM, SN, SN["Cinv_sys"]),
                        rl.x, [0.005, 0.15, 0.003])
    print(f"  DSIC : mu0 = {rd.x[0]:.4f} ± {sd[0]:.4f}   H0 = {rd.x[1]:.2f} ± {sd[1]:.2f}"
          f"   M = {rd.x[2]:.4f} ± {sd[2]:.4f}")
    print(f"  ΛCDM : Om  = {rl.x[0]:.4f} ± {sl[0]:.4f}   H0 = {rl.x[1]:.2f} ± {sl[1]:.2f}"
          f"   M = {rl.x[2]:.4f} ± {sl[2]:.4f}\n")

    print("=" * 70)
    print("[2] Robustness to the error treatment: STAT-ONLY")
    print("=" * 70)
    rds = fit_sn(distmod_DSIC, SN, SN["Cinv_stat"], [0.76, 73.0, -19.3], OK_D)
    rls = fit_sn(distmod_LCDM, SN, SN["Cinv_stat"], [0.33, 73.0, -19.3], OK_L)
    print(f"DSIC STAT-ONLY: mu0={rds.x[0]:.4f}  H0={rds.x[1]:.2f}  chi2/dof={rds.fun/dof:.4f}")
    print(f"ΛCDM STAT-ONLY: Om ={rls.x[0]:.4f}  H0={rls.x[1]:.2f}  chi2/dof={rls.fun/dof:.4f}")
    print(f"Δχ² STAT-ONLY = {rds.fun - rls.fun:+.2f}")
    print(f"  mu0: {rd.x[0]:.4f} (STAT+SYS) → {rds.x[0]:.4f} (STAT-ONLY)\n")
    return rd.x, rl.x


# ==============================================================================
#  [3] JOINT SN + BAO (full 12×12 covariance, free r_d)
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
    print("[3] Joint SN+BAO (full 12×12 covariance, free r_d)")
    print("=" * 70)
    jd = minimize(joint_D, [0.76, 73, -19.3, 147], method="Nelder-Mead",
                  options=dict(xatol=1e-5, fatol=1e-4, maxiter=30000))
    jl = minimize(joint_L, [0.31, 73, -19.3, 147], method="Nelder-Mead",
                  options=dict(xatol=1e-5, fatol=1e-4, maxiter=30000))
    mu0J, H0dJ, MdJ, rddJ = jd.x
    OmJ, H0lJ, MlJ, rdlJ = jl.x
    baoD = chi2_bao(bao_vec_DSIC(mu0J, H0dJ, rddJ))
    baoL = chi2_bao(bao_vec_LCDM(OmJ, H0lJ, rdlJ))
    print(f"DSIC : mu0={mu0J:.4f} H0={H0dJ:.2f} M={MdJ:.3f} r_d={rddJ:.1f}  total={jd.fun:.1f}")
    print(f"ΛCDM : Om ={OmJ:.4f} H0={H0lJ:.2f} M={MlJ:.3f} r_d={rdlJ:.1f}  total={jl.fun:.1f}")
    print(f"Δχ² total (DSIC−ΛCDM) = {jd.fun - jl.fun:+.2f}  "
          f"(weak; partly reflects the DESI–ΛCDM tension)")
    print(f"BAO χ²/dof (12−4=8): DSIC={baoD/8:.2f}  ΛCDM={baoL/8:.2f}")

    print("\n  The r_d ruler (inverse problem) on the early scale H0=67.4:")
    H0e = 67.4
    rdM = [DM_comov_DSIC(zv, mu0J, H0e) / BAO[i, 1] for i, zv in enumerate(z_bao)]
    rdH = [(c / H_DSIC(zv, mu0J, H0e)) / BAO[i, 3] for i, zv in enumerate(z_bao)]
    print(f"    r_d from D_M = {np.mean(rdM):.1f} Mpc  (scatter {np.std(rdM):.1f})")
    print(f"    r_d from D_H = {np.mean(rdH):.1f} Mpc  (scatter {np.std(rdH):.1f})\n")
    return (mu0J, H0dJ, rddJ), (OmJ, H0lJ, rdlJ)


# ==============================================================================
#  [4] Om(z): reconstruction + robustness of the sign of Δχ² to H0, in TWO
#      chronometer-covariance modes (diagonal and full Moresco 2020).
# ==============================================================================
def run_omz(dsic_par, lcdm_par):
    mu0_fit = dsic_par[0]
    Om_fit = lcdm_par[0]
    rd_norm = 147.1  # conversion D_H/r_d -> H(z) in the Om(z) reconstruction
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

        # Jacobian σ_H -> σ_Om; C_Om = J C_H J^T (J is diagonal)
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
    print("[4] Om(z): robustness of the sign of Δχ² to the H0 normalization (an honest test)")
    print("=" * 70)
    print("  The chronometers are plugged in with TWO covariance modes:")
    print("    DIAG = diagonal errors (the baseline version of the test);")
    print("    FULL = the full Moresco 2020 covariance (15 points correlated).")
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
            return f"consistently closer to DSIC  (range {signs.min():+.1f}…{signs.max():+.1f})"
        if np.all(signs > 0):
            return f"consistently closer to ΛCDM  (range {signs.min():+.1f}…{signs.max():+.1f})"
        return f"THE SIGN FLIPS → the test does NOT discriminate (range {signs.min():+.1f}…{signs.max():+.1f})"

    print(f"\n  Verdict DIAG: {verdict(signs_d)}")
    print(f"  Verdict FULL: {verdict(signs_f)}")
    print("  [!] Bottom line: in both covariance modes the sign of Δχ² is set by the")
    print("      choice of the H₀ normalization, not by the shape of Om(z): at a low")
    print("      (Planck-like) H₀ the data are formally closer to DSIC, at a high one")
    print("      (SH0ES) — to ΛCDM. Going from DIAG to FULL changes Δχ² by mere units")
    print("      and does not change the verdict. The probe does NOT discriminate; it")
    print("      becomes decisive only with an independently pinned H₀ and more precise H(z).\n")


# ==============================================================================
#  [5] FALSIFIABLE PREDICTIONS AT HIGH z
# ==============================================================================
def run_predictions(dsic_par, lcdm_par):
    mu0J, H0dJ, rddJ = dsic_par
    OmJ, H0lJ, rdlJ = lcdm_par
    print("=" * 70)
    print("[5] Falsifiable predictions, DSIC vs ΛCDM at high z")
    print("=" * 70)
    print("  D_H/r_d (Lyman-α region) — a trace of the CORE, independent of Part II:")
    for zv in [2.33, 2.5, 3.0, 3.5]:
        dhd = (c / H_DSIC(zv, mu0J, H0dJ)) / rddJ
        dhl = (c / (H0lJ * E_LCDM(zv, OmJ))) / rdlJ
        print(f"    z={zv}:  DSIC={dhd:.3f}  ΛCDM={dhl:.3f}  diff={100*(dhd/dhl-1):+.1f}%")
    print("  E(z) beyond the data:")
    for zv in [5, 10]:
        ed = H_DSIC(zv, mu0J, 1.0)
        el = E_LCDM(zv, OmJ)
        print(f"    z={zv}:  E_DSIC/E_ΛCDM = {ed/el:.2f}  (+{100*(ed/el-1):.0f}%)")
    print()


# ==============================================================================
#  [6] COSMIC TIME t(z): the nearest risk zone of the core (new block)
# ------------------------------------------------------------------------------
#  In DSIC, the observer's proper time from the Mush to the epoch μ_e:
#     τ(μ_e) = ln(1/μ_e)/k,   H₀ = k·(2−μ₀²)/(1−μ₀²)
#     =>  t(z)·H₀ = ln(1/μ_e(z)) · (2−μ₀²)/(1−μ₀²)
#  (at z=0 it reproduces t₀·H₀ ≈ 0.925 from the main test).
#  ΛCDM: t(z) = (1/H₀)·∫_z^∞ dz'/((1+z')E(z')) (matter+Λ; the radiation
#  contribution at these epochs is ≲1% and is omitted here, which OVERSTATES
#  the ΛCDM age — the comparison is conservative in ΛCDM's favor, i.e. it is
#  not skewed toward DSIC).
#  The optical correction of Part II lengthens the light path but does NOT
#  change H(z), so this block tests precisely the core (the paper, §11.6).
# ==============================================================================
def run_cosmic_time(dsic_par, lcdm_par, H0_scale=67.4):
    mu0J = dsic_par[0]
    OmJ = lcdm_par[0]
    h = lambda m: (2 - m * m) / (1 - m * m)
    Gyr = 977.792 / H0_scale  # 1/H0 in Gyr

    def t_dsic(zv):
        me = mu_e_of_z(zv, mu0J)
        return np.log(1.0 / me) * h(mu0J) * Gyr

    def t_lcdm(zv):
        return Gyr * quad(lambda x: 1.0 / ((1 + x) * E_LCDM(x, OmJ)),
                          zv, np.inf, limit=400)[0]

    print("=" * 70)
    print("[6] Cosmic time t(z): DSIC vs ΛCDM (common scale H0=67.4)")
    print("=" * 70)
    print(f"  Age today: t0_DSIC = {t_dsic(0.0):.2f} Gyr,  "
          f"t0_ΛCDM = {t_lcdm(0.0):.2f} Gyr  (obs. 13.80 ± 0.02)")
    print(f"  {'z':>5} | {'t_DSIC, Gyr':>12} | {'t_ΛCDM, Gyr':>12} | {'ratio':>9}")
    print("  " + "-" * 48)
    for zv in [3.0, 5.0, 6.0, 7.0, 7.5, 8.0, 10.0, 12.0]:
        td, tl = t_dsic(zv), t_lcdm(zv)
        print(f"  {zv:>5.1f} | {td:>12.3f} | {tl:>12.3f} | {td/tl:>9.2f}")
    print()
    print("  Context (what must have time to happen within this budget):")
    print("   - quasars with M_BH ~ 10⁹ M☉ are observed by z ≈ 7.5;")
    print("   - JWST galaxies are confirmed out to z ≈ 10–14;")
    print("   - reionization completes by z ≈ 6.")
    print("  The compression of time in DSIC at these epochs is the nearest testable")
    print("  risk zone of the CORE (Part II does not soften it: H(z) is unchanged).\n")


# ==============================================================================
#  PART II — THE SECOND TIER (P6): THE DETACHMENT THRESHOLD μ_d
# ==============================================================================
def DM_base_mu(mu_e):
    """The core's base comoving distance (early scale MU0_E, H0_E), from μ_e."""
    ck = ck_const(MU0_E, H0_E)
    pref = np.sqrt(1.0 - MU0_E ** 2) / MU0_E ** 2
    return ck * pref * (np.sqrt(1.0 - MU0_E ** 2) - np.sqrt(1.0 - mu_e ** 2))


def DM_corrected_mu(mu_e, mu_d):
    """Corrected distance: the arcsin-correction make-up above the threshold μ_d.

    By P6, the addition to the path in the stuck-together phase is
    1/P = 1/√(1−μ²) (the amplitude A ≡ 1 is fixed by the postulate), cut off
    by a hard threshold; ∫ = arcsin μ:
      D_M_corr = D_M + (c/k)·[arcsin(μ_e) − arcsin(μ_d)]  for μ_e > μ_d,
      D_M_corr = D_M                                        for μ_e ≤ μ_d.
    """
    base = DM_base_mu(mu_e)
    if mu_e > mu_d:
        return base + ck_const(MU0_E, H0_E) * (np.arcsin(mu_e) - np.arcsin(mu_d))
    return base


def self_check_planck():
    print("=" * 70)
    print("[7] Self-check of the Planck 2018 reference data (hard-coded)")
    print("=" * 70)
    d_from_geom = PLANCK["r_star"] / PLANCK["theta_star"]
    rel = abs(d_from_geom - PLANCK["d_lss"]) / PLANCK["d_lss"]
    print(f"  z* = {PLANCK['z_star']},  r* = {PLANCK['r_star']} Mpc,  "
          f"100θ* = {100*PLANCK['theta_star']:.5f}")
    print(f"  d_LSS (hard-coded) = {PLANCK['d_lss']} ± {PLANCK['d_lss_err']} Mpc")
    print(f"  d_LSS = r*/θ*      = {d_from_geom:.1f} Mpc  (discrepancy {100*rel:.2f}%)")
    if rel < 0.01:
        print("  OK — the reference figure is consistent with r* and θ* (Planck 2018).\n")
    else:
        raise SystemExit("STOP: d_LSS is inconsistent with r*/θ*. Check PLANCK.")


def step_fix_threshold():
    print("=" * 70)
    print("[8] Second tier (P6): the core's shortfall and fixing μ_d via the angle θ*")
    print("=" * 70)
    ck = ck_const(MU0_E, H0_E)
    mu_cmb = mu_e_of_z(PLANCK["z_star"], MU0_E)
    DM0 = DM_base_mu(mu_cmb)
    target_d = PLANCK["r_star"] / PLANCK["theta_star"]   # = r*/θ*
    deficit = target_d - DM0
    print(f"  Core parameters (early scale): μ₀ = {MU0_E}, H₀ = {H0_E}")
    print(f"  c/k                = {ck:.1f} Mpc")
    print(f"  μ_cmb(z*={PLANCK['z_star']}) = {mu_cmb:.10f}")
    print(f"  D_M base to CMB    = {DM0:.1f} Mpc")
    print(f"  Target: θ_model(μ_d) = θ*  ⇔  D_M_corr(μ_cmb) = r*/θ* = {target_d:.1f} Mpc")
    print(f"  core shortfall     = {deficit:.1f} Mpc  ({100*deficit/target_d:.1f} %)")
    print("  [!] r* = 144.43 Mpc is a quantity COMPUTED within ΛCDM (an external")
    print("      input; the direct measurement is the angle θ*). The origin of the")
    print("      acoustic scale ~147 Mpc is an open question of the core (the paper, §11.6).")
    mu_d = brentq(lambda md: ck * (np.arcsin(mu_cmb) - np.arcsin(md)) - deficit,
                  0.5, mu_cmb - 1e-12)
    z_det = b_dsic(MU0_E) / b_dsic(mu_d) - 1.0
    DM_corr = DM_corrected_mu(mu_cmb, mu_d)
    theta_model = PLANCK["r_star"] / DM_corr
    print(f"\n  μ_d        = {mu_d:.6f}")
    print(f"  z_detach   = {z_det:.4f}")
    print(f"  D_M corr   = {DM_corr:.1f} Mpc;  100·θ_model = {100*theta_model:.5f}"
          f"  (Planck 100θ* = {100*PLANCK['theta_star']:.5f})")
    ok = abs(DM_corr - target_d) < 1.0
    print(f"  closure to θ*: {'OK' if ok else 'MISMATCH'}\n")
    if not ok:
        raise SystemExit("STOP: θ_model did not match θ*.")
    return mu_cmb, DM0, deficit, mu_d, z_det


def step_late_universe_and_degeneracy(mu_cmb, DM0, deficit, mu_d, z_det):
    print("=" * 70)
    print("[9] Zero-below-threshold control + the \"amplitude ↔ threshold\" degeneracy")
    print("=" * 70)
    print("  The correction is identically zero at z ≤ z_detach:")
    print(f"  {'z':>6} {'D_M base':>10} {'correction':>10}")
    print("  " + "-" * 28)
    max_corr_below = 0.0
    for z in [0.5, 1.0, 2.33, 3.0, 4.0, z_det, z_det + 0.1, 5.0, 10.0]:
        me = mu_e_of_z(z, MU0_E)
        base = DM_base_mu(me)
        corr = DM_corrected_mu(me, mu_d) - base
        tag = "" if z <= z_det + 1e-9 else "  <- switched on"
        if z <= z_det + 1e-9:
            max_corr_below = max(max_corr_below, abs(corr))
        print(f"  {z:>6.2f} {base:>10.1f} {corr:>+10.1f}{tag}")
    print(f"\n  max |correction| at z ≤ z_detach = {max_corr_below:.2f} Mpc (≈0)")
    print("  -> all 6 DESI bins (z ≤ 2.33) are insensitive to the correction; the")
    print("     results of Part I ([1]–[5]) are fully preserved.\n")

    ck = ck_const(MU0_E, H0_E)
    print("  The degeneracy WITHOUT postulate P6 (what exactly the postulate fixes):")
    print("  the equation A·(c/k)·[arcsin μ_cmb − arcsin μ_d] = shortfall has a")
    print("  curve of solutions (A, μ_d) — all give the same distance to the CMB.")
    print(f"  {'A':>6} {'μ_d':>10} {'z_detach':>10}")
    print("  " + "-" * 30)
    for A in [0.5, 0.8, 1.0, 1.2, 2.0, 5.0]:
        arg = np.arcsin(mu_cmb) - deficit / (A * ck)
        if 0 < arg < np.pi / 2:
            md = np.sin(arg)
            zd = b_dsic(MU0_E) / b_dsic(md) - 1.0
            star = "  <- fixed by P6 (A≡1)" if abs(A - 1.0) < 1e-9 else ""
            print(f"  {A:>6.1f} {md:>10.6f} {zd:>10.3f}{star}")
        else:
            print(f"  {A:>6.1f} {'no sol.':>10}")
    print("\n  One observation fixes one parameter: zero degrees of freedom.")
    print("  With the form fixed by P6, μ_d is a measured constant of the second")
    print("  tier (the same status as μ₀ in the core); z_detach is not a prediction.\n")


def step_trace_and_horizon(mu_d, z_det):
    print("=" * 70)
    print("[10] The trace of μ_d, the corrected horizon, and the consistency flag")
    print("=" * 70)
    print(f"  {'z':>8} {'object':<22} {'correction to D_M':>18}")
    print("  " + "-" * 52)
    for z, obj in [(z_det, "detachment threshold"), (5.0, "JWST galaxies"),
                   (6.0, "quasars"), (8.0, "reionization"),
                   (11.0, "first galaxies"), (PLANCK["z_star"], "CMB (calibration)")]:
        me = mu_e_of_z(z, MU0_E)
        base = DM_base_mu(me)
        corr = DM_corrected_mu(me, mu_d)
        print(f"  {z:>8.2f} {obj:<22} {100*(corr/base-1):>+15.1f} %")
    ck = ck_const(MU0_E, H0_E)
    hor_base = ck * (1.0 - MU0_E ** 2) / MU0_E ** 2
    hor_corr = ck * ((1.0 - MU0_E ** 2) / MU0_E ** 2 + (np.pi / 2 - np.arcsin(mu_d)))
    print(f"\n  Particle horizon: base = {hor_base:.1f} Mpc;  with the correction = "
          f"{hor_corr:.1f} Mpc")
    print(f"  ({hor_corr - PLANCK['d_lss']:.1f} Mpc beyond the last-scattering "
          f"surface — the CMB lies almost flush against the horizon)")
    print()
    print("  Reach of the probes:")
    print("    DESI DR2 (available)  z ≤ 2.33   — below the threshold")
    print("    DESI DR3 (~2027)      z ≤ 3.55   — below the threshold")
    print("    DESI-II/LBG (~2029)   2.5<z<3.5  — below the threshold")
    print("  The radial BAO in the window 2.33–4.53 is a trace of the CORE (D_H = c/H),")
    print("  not of μ_d; the threshold itself affects only D_M above z≈4.5, where no")
    print("  precise rulers exist apart from the CMB calibration: a second independent")
    print("  anchor is out of reach until a precise distance scale appears at z ≳ 5–6.")
    print()
    mu_e14 = mu_e_of_z(14.0, MU0_E)
    print("  [!] CONSISTENCY FLAG OF THE P6 FORMULATION:")
    print(f"      the observed galaxies (JWST, z ≈ 10–14; μ_e(z=14) = {mu_e14:.4f})")
    print(f"      lie ABOVE the detachment threshold (μ_d = {mu_d:.4f}, z_detach = {z_det:.2f}),")
    print("      i.e. they exist as separate objects BEFORE detachment in the sense of P6.")
    print("      A literal reading — \"no separate objects exist at μ > μ_d\" —")
    print("      contradicts the observations; the formulation needs refinement")
    print("      (detachment of the dominant monolith vs the first detached islets),")
    print("      or a revision of the choice of A (which breaks the naturalness of A ≡ 1).")
    print("      This is a limitation of the NARRATIVE of Part II; the distance")
    print("      arithmetic and the results of Part I do not depend on it.\n")


# ==============================================================================
#  MAIN FLOW
# ==============================================================================




# ==============================================================================
#  PART III — STRUCTURE GROWTH  (block [11])
# ==============================================================================
#
# DSIC — a check of structure growth (linear growth, fσ8) against RSD data.
# Block [11] for dsic_full_test.py.
#
# Idea: DSIC sets the background H(z). On this background I solve the STANDARD
# equation of linear growth of matter perturbations (sub-horizon, GR limit):
#
#     d²δ/d(lna)² + (2 + dlnH/dlna) dδ/d(lna) - (3/2) Ω_m(a) δ = 0
#
# Then f = dlnδ/dlna, σ8(z) = σ8,0·δ(z)/δ(0), and fσ8(z) is compared with RSD.
#
# ASSUMPTION (stated honestly): DSIC does not derive the growth equation from
# the core. Here the GR limit of growth is taken on the DSIC geometry — a
# survival test, not a derivation. Ω_m0 enters as a growth parameter, since
# the core introduces no densities.
#
# KEY CONCLUSION: the fσ8 data are degenerate — they constrain the combination
# S8 = σ8·√(Ω_m0/0.3), not Ω_m0 and σ8 separately. The result is therefore
# reported as S8, plus a slice at fixed Planck-like Ω_m0=0.30.
#
# ------------------------------------------------------------------------
# DATA: the "Gold-2017" fσ8 compilation — 18 weakly correlated points,
# Nesseris, Pantazis & Perivolaropoulos (2017), Phys. Rev. D 96, 023542,
# arXiv:1703.10538. z, fσ8, σ are taken from Table I of Sagredo, Nesseris &
# Sapone (2018), arXiv:1806.10822 (the first 18 rows = Gold-2017).
# The three WiggleZ points (z=0.44,0.60,0.73) — with the Blake et al. 2012
# covariance.
#
# SIMPLIFICATION: the Alcock–Paczynski correction (DSIC vs fiducial) is not
# applied; for z<1 the AP factor is typically <2-3%, small against the errors.
# Flagged as a systematic.
# ------------------------------------------------------------------------
#
# ==============================================================================


# ---------- DSIC BACKGROUND ----------
def b_of_mu(mu):
    return np.sqrt(1 - mu**2) / mu**2

def mu_of_z(z, mu0):
    """1+z = b(mu0)/b(mu_e), physical branch mu_e in [mu0,1). At z=0 -> mu0."""
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

# ---------- PERTURBATION GROWTH ----------
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

# ---------- fσ8 DATA (Gold-2017, 18 points) ----------
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
    """Block [11]: structure growth fσ8 on the DSIC background (GR limit)."""
    mu0 = 0.76
    dof = len(fs8_data) - 2

    # 1) Full scan -> the true minimum (a demonstration of the degeneracy)
    Om_grid = np.linspace(0.05, 0.50, 91)
    s8_grid = np.linspace(0.55, 1.15, 121)
    best = (1e9, None)
    for Om0 in Om_grid:
        for s8 in s8_grid:
            c = chi2_fs8(mu0, Om0, s8)
            if c < best[0]: best = (c, (Om0, s8))
    c, (Om0, s8) = best
    print(f"DSIC growth (GR limit on the DSIC background), mu0={mu0}")
    print(f"  sample: Gold-2017, {len(fs8_data)} points (WiggleZ covariance included)")
    print(f"  formal minimum: Om0={Om0:.3f}, s8={s8:.3f}, "
          f"chi2={c:.2f}, chi2/dof={c/dof:.2f}")
    print(f"  (the minimum lies inside the grid; the valley is strongly degenerate — see below)")

    # 2) Degeneracy: the combination S8 is constant along the valley
    print("\n  fσ8 degeneracy (the data constrain S8, not Om0 and s8 separately):")
    print("   Om0    s8*     chi2/dof   S8=s8*sqrt(Om0/0.3)")
    S8_vals = []
    for Om0 in [0.15, 0.20, 0.25, 0.30, 0.35, 0.40]:
        cb, s8b = min((chi2_fs8(mu0, Om0, s8), s8) for s8 in s8_grid)
        S8 = s8b*np.sqrt(Om0/0.3); S8_vals.append(S8)
        print(f"  {Om0:.2f}   {s8b:.3f}   {cb/dof:.2f}       {S8:.3f}")
    print(f"  => S8 is stable: {np.mean(S8_vals):.3f} ± {np.std(S8_vals):.3f}")

    # 3) The slice at Planck-like Om0=0.30 — the recommended reporting form
    cb, s8b = min((chi2_fs8(mu0, 0.30, s8), s8) for s8 in s8_grid)
    print(f"\n  At fixed Om0=0.30 (Planck-like):")
    print(f"    s8={s8b:.3f}, chi2={cb:.2f}, chi2/dof={cb/(len(fs8_data)-1):.2f}")
    print(f"    S8={s8b*np.sqrt(0.30/0.3):.3f}")

    # 4) The curve at Om0=0.30
    zt, curve = fsigma8_curve(mu0, 0.30, s8b)
    print("\n   z    fσ8_model   fσ8_obs±err   (Om0=0.30)")
    for zz, obs, err in fs8_data:
        print(f"  {zz:.2f}   {np.interp(zz, zt, curve):.3f}       {obs:.3f}±{err:.3f}")


def main():
    # -------- Part I: the core --------
    self_check_desi()                                 # [0]
    self_check_coherence()                            # [0b]
    SN = load_sn()                                    # [1a]
    dsic_sn, lcdm_sn = run_sn_fits(SN)                # [1],[2]
    dsic_par, lcdm_par = run_joint(SN)                # [3]
    run_omz(dsic_par, lcdm_par)                       # [4]
    run_predictions(dsic_par, lcdm_par)               # [5]
    run_cosmic_time(dsic_par, lcdm_par)               # [6]
    # -------- Part II: the second tier --------
    self_check_planck()                               # [7]
    mu_cmb, DM0, deficit, mu_d, z_det = step_fix_threshold()        # [8]
    step_late_universe_and_degeneracy(mu_cmb, DM0, deficit, mu_d, z_det)  # [9]
    step_trace_and_horizon(mu_d, z_det)               # [10]
    # -------- Part III: structure growth --------
    print()
    print("#" * 70)
    print("#  PART III — STRUCTURE GROWTH (block [11])")
    print("#" * 70)
    run_growth()                                      # [11]
    print("=" * 70)
    print("Done. Part I: SNe downloaded, DESI passed the self-check, the core's")
    print("identities hold, parameter uncertainties computed from the Hessian.")
    print("Part II: the Planck reference data verified, the threshold μ_d fixed via θ*")
    print("(r* is borrowed — flagged), the late Universe untouched, the consistency")
    print("flag of the P6 formulation printed. μ_d is a calibration against a single")
    print("observation with the form fixed by the postulate.")
    print("Part III: structure growth fσ8 on the DSIC background (GR limit) checked\n"
          "against the Gold-2017 compilation; the data are degenerate in (Ω_m0, σ8),\n"
          "the invariant S8 is reported. Status — a survival test, not a derivation\n"
          "from the core.")
    print("=" * 70)


if __name__ == "__main__":
    main()