# DSIC — Empirical verification results (Part I + Part II + Part III)

**Model:** Dual Scale Inversion Cosmology (DSIC)
**Author:** Azamat Zaseev
**Script:** `dsic_full_test.py` (single, self-contained; blocks `[0]`–`[11]`)
**Comparison:** DSIC versus flat ΛCDM on identical data, with an equal number of free parameters.

> **In brief.** The internal identities of the core pass machine verification to `10⁻⁸` accuracy. On the expansion geometry of the late Universe (`z < 2.3`), DSIC with the single parameter `μ₀ = 0.7600 ± 0.0091` is statistically indistinguishable from ΛCDM on the supernovae and slightly preferred on the joint SN+BAO fit; the ΛCDM branch of the pipeline reproduces the published Pantheon+SH0ES values — a pipeline validation. The `Om(z)` diagnostic on current data **does not discriminate** between the models (the sign of Δχ² depends on the `H₀` normalization). The second tier (the detachment threshold `μ_d`) closes the distance to the CMB exactly, leaving the late Universe identically untouched; its status is a calibration with the form fixed by postulate P6. Structure growth (`fσ₈`) on the DSIC background in the GR limit describes the Gold-2017 compilation no worse than ΛCDM (`χ²/dof ≈ 0.7`); the data are degenerate in `(Ω_m0, σ₈)`, so the invariant `S8 = 0.761 ± 0.013` is reported — the status is a survival test, not a derivation from the core. The most stringent risk zone, stated as a number, is the compressed cosmic time at `z = 5–12` (at `z = 7.5`, 0.42 Gyr are available versus 0.71 in ΛCDM). The results are reproducible: the script downloads the public catalogs itself and computes χ² with the full covariance.

---

## 1. Summary of results

| Probe | Data / block | Metric | DSIC | ΛCDM | Verdict |
|---|---|---|---|---|---|
| Core consistency | identities, `[0b]` | max deviation | `≤ 5.7·10⁻⁸` | — | all identities hold |
| Supernovae (anchor + explicit `M`) | Pantheon+SH0ES (1657 SNe, 77 calib.), `[1]` | `χ²/dof` | 0.8783 | 0.8783 | indistinguishable (Δχ² = +0.05) |
| Parameter uncertainties | χ² Hessian, `[1]` | 1σ | `μ₀ = 0.7600 ± 0.0091` | `Ωm = 0.3308 ± 0.0180` | parameter measured (~1.2%) |
| Robustness to error treatment | STAT-ONLY, `[2]` | Δχ² | — | — | robust (Δχ² = −0.39) |
| SN + BAO joint | + DESI DR2 (6 bins), `[3]` | Δχ² (total) | — | — | DSIC better by −5.0 (weak) |
| BAO `D_M`–`D_H` | DESI DR2 (12×12 covariance), `[3]` | `χ²/dof` | 1.04 | 1.36 | DSIC no worse |
| `Om(z)` diagnostic | chronometers + DESI `D_H/r_d`, `[4]` | sign of Δχ² | — | — | **does not discriminate** (depends on `H₀`) |
| High `z` (`z > 3`) | prediction, `[5]` | `E(z)`, `D_H/r_d` | diverges | — | falsifiable |
| Cosmic time | `t(z)` at `z = 5–12`, `[6]` | fraction of ΛCDM | 0.68 → 0.48 | 1 | **the core's risk zone** (as a number) |
| Planck reference data | `d_LSS = r*/θ*`, `[7]` | consistency | — | — | verified (0.00%) |
| Distance to the CMB | closure via `θ*`, `[8]` | `100·θ_model` | 1.04110 | — | exact (by construction of `μ_d`) |
| Late Universe under `μ_d` | correction at `z ≤ 4.53`, `[9]` | Mpc | 0.00 | — | untouched (identically) |
| Degrees of freedom of `μ_d` | `A ↔ μ_d` degeneracy, `[9]` | dof | 0 | — | calibration with the form from P6 |
| P6 consistency | JWST galaxies vs the threshold, `[10]` | — | — | — | flag: requires the "average medium" formulation |
| Structure growth `fσ₈` | Gold-2017 (18 points), `[11]` | `χ²/dof` | 0.70 | ~0.7 | DSIC background no worse (survival test) |
| The `S8` invariant | `(Ω_m0, σ₈)` degeneracy, `[11]` | 1σ | `0.761 ± 0.013` | `0.76–0.83` (obs.) | within the observed range |

In all fits, DSIC and ΛCDM have **the same number of free parameters**, so Δχ² = ΔAIC = ΔBIC.

---

## 2. Core consistency: machine verification of the identities

Block `[0b]` verifies that the observable physics of the core is internally closed — distances, rate, and kinematics follow from one geometry without inconsistencies:

| Identity | Result |
|---|---|
| `b(μ)` strictly decreasing on `(μ₀, 1)` — the inversion `z → μ_e` is well defined | OK |
| `D_H = dχ/dz` (comoving `D_H = c/H`) | max deviation `5.7·10⁻⁸` |
| closed-form `χ(μ)` = numerical `∫ c dz′/H(z′)` | max deviation `4.4·10⁻¹²` |
| `q(μ_cross) = 0` at the exact root `μ_cross = √(3−√5) = 0.874032` | OK; `z_t = 0.769` |
| age: `t₀·H₀ = ln(1/μ₀)·(2−μ₀²)/(1−μ₀²) = 0.9241` | `t₀ = 13.41 Gyr` at `H₀ = 67.4` |

The closed-form formulas of the paper and the numerical geometry agree at the level of machine precision: the expansion rate, the distances, and the transition point are one and the same geometry, not independent fits.

---

## 3. Supernovae: the distance diagram

### 3.1 Anchored fit: calibrators + `M` as an explicit parameter

The full Pantheon+SH0ES catalog (1701 SNe; 1657 enter the fit, including 77 calibrators with Cepheid distances as the absolute anchor). The calibrators are not discarded — they pin the absolute scale; `M` is fitted explicitly; `H₀` is a genuine inference from the anchor. The full STAT+SYS covariance matrix. Three free parameters for each model (`μ₀`/`Ωm`, `H₀`, `M`).

```
DSIC : μ₀ = 0.7600 ± 0.0091,  H₀ = 73.35 ± 1.01,  M = −19.243 ± 0.030,  χ²/dof = 0.8783
ΛCDM : Ωm = 0.3308 ± 0.0180,  H₀ = 73.53 ± 1.02,  M = −19.244 ± 0.030,  χ²/dof = 0.8783

Δχ² = ΔAIC = ΔBIC = +0.05
```

The uncertainties are 1σ from the Hessian of χ² at the minimum (`C = 2H⁻¹`); the model parameter is measured to ~1.2%.

On fully calibrated data, DSIC is statistically **indistinguishable** from ΛCDM (`|Δχ²| ≪ 2`). The one-parameter geometry of `μ₀` reproduces the same distance–redshift relation as ΛCDM with the density parameter `Ωm` — without dark energy. `H₀ ≈ 73` inherits the familiar "Hubble tension" against the early Universe (`H₀ ≈ 67`) — identically in both models: DSIC creates no new tension and removes no existing one.

> **Pipeline validation.** The ΛCDM branch of the same pipeline reproduces the published values of the Pantheon+SH0ES analysis (`Ωm ≈ 0.33 ± 0.018`, `H₀ ≈ 73.5 ± 1.0`). Both models are computed by one code on one dataset; the agreement of the control branch with the literature confirms the correctness of the setup.

### 3.2 Robustness: STAT-ONLY versus STAT+SYS

Control: the same anchored fit with the matrix of statistical errors only (no systematics). If the parameters and the shape barely shift, the result does not rest on the systematics.

| | `μ₀` / `Ωm` | `H₀` | `χ²/dof` |
|---|---|---|---|
| DSIC, STAT+SYS | 0.7600 | 73.35 | 0.8783 |
| ΛCDM, STAT+SYS | 0.3308 | 73.53 | 0.8783 |
| DSIC, STAT-ONLY | 0.7661 | 73.04 | 0.9242 |
| ΛCDM, STAT-ONLY | 0.3466 | 73.22 | 0.9244 |

```
Δχ² (DSIC − ΛCDM), STAT-ONLY = −0.39
```

`μ₀` shifts only from 0.7600 to 0.7661 (within 1σ), the shape of the relation is preserved, and DSIC remains level with ΛCDM. The rise of `χ²/dof` when the systematics are removed is an expected effect, identical for both models.

---

## 4. Joint SN + BAO fit

Pantheon+SH0ES supernovae (1657 SNe, anchored) and DESI DR2 BAO (6 bins, `0.5 < z < 2.33`) with the **full 12×12 covariance** of `D_M`–`D_H` (the correlation between the transverse and radial scales within each bin is taken into account). The sound-horizon scale `r_d` is a free parameter (it does not cancel via a ratio). Four free parameters for each model.

```
DSIC : μ₀ = 0.7610,  H₀ = 73.33,  M = −19.243,  r_d = 136.5 Mpc,  total = 1461.0
ΛCDM : Ωm = 0.3029,  H₀ = 73.75,  M = −19.247,  r_d = 137.1 Mpc,  total = 1466.0

Δχ² total (DSIC − ΛCDM) = −5.0
BAO χ²/dof (dof = 12−4 = 8):  DSIC = 1.04,  ΛCDM = 1.36
```

On the joint data, DSIC describes the expansion geometry **no worse than, and formally slightly better than,** ΛCDM (Δχ² = −5.0; on the Jeffreys scale this is a weak, inconclusive preference, partly reflecting the known DESI–ΛCDM tension). The consistency on BAO (`χ²/dof ≈ 1`) shows that the same `μ₀` geometry that describes the supernovae simultaneously reproduces the baryon scale.

> **A remark on `r_d`.** In the joint fit, `r_d ≈ 136–137 Mpc` lies below the standard `≈ 147 Mpc`; this is a consequence of the high `H₀ ≈ 73` from the SH0ES anchor (the product `r_d·H₀` is fixed by the data). On the early scale `H₀ = 67.4`, the inverse problem from the same `D_M/r_d`, `D_H/r_d` yields a scale constant in `z`: `r_d = 147.6 Mpc` from `D_M` (scatter 1.7) and `r_d = 149.8 Mpc` from `D_H` (scatter 1.8) — consistent between the transverse and radial measurements and within ~1.5% of the standard value. The value of `r_d` shifts together with `H₀` through the "Hubble tension" identically in both models.

---

## 5. The `Om(z)` diagnostic — an honest result

The diagnostic `Om(z) = (E² − 1)/((1+z)³ − 1)`, `E = H/H₀`. For flat ΛCDM, `Om(z) = Ωm = const` (a horizontal line); DSIC predicts a U-shaped `Om(z)` with a minimum near `z ≈ 2`. `Om(z)` is built **from the data** (Moresco chronometer `H(z)` + DESI `D_H/r_d`) and compared with both curves. The chronometers are plugged in with two covariance modes: DIAG (diagonal errors) and FULL (the full Moresco et al. 2020 covariance for the 15 correlated points; the remaining 18 points of the compilation are diagonal).

```
  H0 |   Δχ² DIAG |   Δχ² FULL
--------------------------------
  64 |      −25.0 |      −24.2
  66 |      −15.6 |      −15.2
  68 |       −4.6 |       −4.7
  70 |       +8.3 |       +7.6
  72 |      +23.1 |      +21.6
  74 |      +39.9 |      +37.7
  76 |      +58.9 |      +55.8
```

**The sign of Δχ² changes with the choice of `H₀`** in both modes (the flip occurs near `H₀ ≈ 68.7`): at a low (Planck-like) normalization the data are formally closer to DSIC, at a high one (SH0ES) — to ΛCDM. The result is entirely determined by the `H₀` normalization, not by the shape of `Om(z)`; going from DIAG to FULL changes Δχ² by only a few units and does not change the verdict.

> **Conclusion.** With current data, the `Om(z)` test **does not discriminate** between DSIC and ΛCDM. Claiming a DSIC victory on this probe is not permissible. The test will become decisive only with an independently pinned `H₀` and more precise `H(z)`.

> **A caveat on the covariance.** The full matrix is a careful reconstruction of the Moresco et al. (2020) method on the official input tables (hard-coded into the script), not a verbatim run of the authors' code. The diagonal matches the publication exactly; the off-diagonal structure reproduces the prescribed recipe of correlated systematic components. This does not affect the qualitative conclusion (the sign depends on `H₀`).

---

## 6. Falsifiable predictions of DSIC

Where data are still scarce, DSIC and ΛCDM diverge far enough for measurements to deliver a verdict. The parameters are from the joint fit (`μ₀ = 0.7610`, `Ωm = 0.3029`). All signals in this section are a trace of the **core** (Part I) and are independent of the second-tier machinery: `D_H = c/H` is determined by `μ₀` alone, and the `μ_d` correction is identically zero at `z ≤ 4.53` (§8.3).

### 6.1 The shape of `Om(z)`

The values are determined by the parameters of the joint fit:

| `z` | `Om_DSIC` | `Om_ΛCDM` |
|---|---|---|
| 0.2 | 0.3566 | 0.3029 |
| 0.5 | 0.3297 | 0.3029 |
| 1.0 | 0.3064 | 0.3029 |
| 2.0 | 0.3032 | 0.3029 |
| 3.0 | 0.3266 | 0.3029 |
| 5.0 | 0.4043 | 0.3029 |

DSIC yields a U-shaped `Om(z)` (a minimum near `z ≈ 2`, a rise toward high `z`) — a qualitative signature that flat ΛCDM does not have. A flat `Om(z)`, measured precisely, would work against DSIC.

### 6.2 Radial BAO at high `z` (Lyman-α)

| `z` | `D_H/r_d` DSIC | `D_H/r_d` ΛCDM | difference |
|---|---|---|---|
| 2.33 | 8.612 | 8.599 | +0.1% |
| 2.5 | 7.978 | 8.013 | −0.4% |
| 3.0 | 6.449 | 6.615 | −2.5% |
| 3.5 | 5.300 | 5.572 | −4.9% |

By `z ≈ 3.5` the departure reaches ~5%. Precise Lyman-α BAO in this range (the window `2.33 < z < 4.53` lies entirely below the detachment threshold) is a direct test of the core.

### 6.3 The expansion rate beyond the data

```
z = 5  : E_DSIC / E_ΛCDM = 1.15   (+15%)
z = 10 : E_DSIC / E_ΛCDM = 1.47   (+47%)
```

DSIC predicts a substantially stiffer early expansion history. The direct consequence for cosmic time is Section 7.

---

## 7. Cosmic time `t(z)` — the core's risk zone, stated as a number

In DSIC, the observer's proper time from the Mush to the epoch `μ_e` is: `t(z)·H₀ = ln(1/μ_e(z))·(2−μ₀²)/(1−μ₀²)`; at `z = 0` it reproduces `t₀·H₀ ≈ 0.925` from block `[0b]`. The comparison is on the common scale `H₀ = 67.4` (`Ωm = 0.3029` from the joint fit, whence the ΛCDM `t₀` here is 13.95 Gyr — both models are put on an equal footing; the ≲1% radiation contribution is omitted, which does not skew the comparison in DSIC's favor):

```
Age today: t₀_DSIC = 13.38 Gyr,  t₀_ΛCDM = 13.95 Gyr
```

| `z` | `t_DSIC`, Gyr | `t_ΛCDM`, Gyr | ratio |
|---|---|---|---|
| 3.0 | 1.727 | 2.184 | 0.79 |
| 5.0 | 0.812 | 1.194 | 0.68 |
| 6.0 | 0.604 | 0.948 | 0.64 |
| 7.0 | 0.467 | 0.776 | 0.60 |
| 7.5 | 0.415 | 0.709 | 0.59 |
| 8.0 | 0.371 | 0.651 | 0.57 |
| 10.0 | 0.250 | 0.482 | 0.52 |
| 12.0 | 0.180 | 0.375 | 0.48 |

Context: quasars with `M_BH ~ 10⁹ M☉` are observed out to `z ≈ 7.5`; JWST galaxies are confirmed out to `z ≈ 10–14`; reionization completes by `z ≈ 6`.

**A Salpeter estimate.** At the Eddington limit, one e-folding of a black hole's mass takes ~45 Myr; growth from a stellar seed of `~10² M☉` to `~10⁹ M☉` requires ~16 e-foldings ≈ 0.72 Gyr. ΛCDM has 0.71 Gyr at its disposal by `z ≈ 7.5` — just barely enough (which already motivates heavy seeds); DSIC has 0.42 Gyr ≈ 9 e-foldings — either direct-collapse seeds of `~10⁵ M☉` with continuous accretion, or super-Eddington regimes, become mandatory. This is not a refutation (super-Eddington accretion is discussed in earnest, not least because of the "overmassive" JWST black holes), but it is **the narrowest point of the core** — the optical correction of Part II does not soften it, since it does not change `H(z)`.

**A caveat on the age today.** The oft-quoted `13.80 ± 0.02 Gyr` is an age derived within ΛCDM from the Planck data; measuring a non-ΛCDM model against it is partly circular. The value `t₀ = 13.38 Gyr` passes the model-independent stellar ages (the oldest globular clusters, `~13.5 ± 0.5 Gyr`).

---

## 8. Part II: the distance to the CMB (the detachment threshold `μ_d`)

> **Status.** `μ_d` is a constant of the second tier: the form of the correction is fixed by postulate P6 (amplitude `A ≡ 1`, a sharp threshold), and its value is calibrated against **one** observation. The verification below **does not prove** `μ_d`; it confirms that the arithmetic is reproducible, the Planck reference figure is correct, and the late Universe is untouched, and it honestly exhibits the degeneracy that is fixed precisely by the postulate, not by the data.

### 8.1 Planck reference data (self-check)

```
z* = 1089.92,  r* = 144.43 Mpc,  100θ* = 1.04110
d_LSS (hard-coded) = 13872.8 ± 25.0 Mpc
d_LSS = r*/θ*      = 13872.8 Mpc   (discrepancy 0.00%)
```

The directly measured quantity is the angle `θ*`. The sound horizon `r* = 144.43 Mpc` is a quantity **computed within ΛCDM** from the physics of the photon–baryon plasma, which DSIC does not have; it is used as an external input, and this is flagged explicitly. The origin of the acoustic scale `~147 Mpc` is an open question of the core.

### 8.2 The core's shortfall and the calibration of the threshold via `θ*`

The core parameters on the early scale are the same as for BAO, the chronometers, and the age: `μ₀ = 0.76`, `H₀ = 67.4`, `c/k = 14978.2 Mpc`.

```
μ_cmb(z* = 1089.92) = 0.9999994681
D_M base to the CMB  = 10936.2 Mpc
shortfall            = 2936.6 Mpc  (21.2%)

calibration:  θ_model(μ_d) = θ*  ⇔  D_M_corr(μ_cmb) = r*/θ* = 13872.8 Mpc
μ_d       = 0.980640
z_detach  = 4.526
100·θ_model = 1.04110   (Planck: 1.04110)  — exact closure
```

### 8.3 The late Universe is untouched; the degeneracy without P6

The correction is identically zero below the threshold:

| `z` | `D_M` base, Mpc | correction, Mpc |
|---|---|---|
| 0.50 | 1932.8 | +0.0 |
| 1.00 | 3384.2 | +0.0 |
| 2.33 | 5792.7 | +0.0 |
| 4.00 | 7335.6 | +0.0 |
| 4.53 | 7653.3 | +0.0 |
| 4.63 | 7707.7 | +49.3 ← switched on |
| 5.00 | 7896.9 | +220.4 ← switched on |
| 10.00 | 9247.2 | +1433.1 ← switched on |

All 6 DESI bins (`z ≤ 2.33`) are insensitive to the correction — the results of Part I (Sections 3–6) are fully preserved.

Without postulate P6 the calibration is degenerate — the equation `A·(c/k)·[arcsin μ_cmb − arcsin μ_d] = shortfall` has a curve of solutions:

| `A` | `μ_d` | `z_detach` |
|---|---|---|
| 0.5 | 0.923705 | 1.506 |
| 0.8 | 0.969868 | 3.344 |
| **1.0** | **0.980640** | **4.526** ← fixed by P6 |
| 1.2 | 0.986514 | 5.690 |
| 2.0 | 0.995097 | 10.266 |
| 5.0 | 0.999190 | 26.922 |

One observation fixes one parameter: **zero degrees of freedom**. With the form fixed by P6, `μ_d` is a measured constant of the second tier (the same status as `μ₀` in the core); `z_detach` is not a prediction.

### 8.4 The threshold's trace, the horizon, and the reach of the probes

| `z` | object | correction to `D_M` |
|---|---|---|
| 4.53 | detachment threshold | +0.0% |
| 5.00 | JWST galaxies | +2.8% |
| 6.00 | quasars | +7.1% |
| 8.00 | reionization | +12.4% |
| 11.00 | first galaxies | +16.6% |
| 1089.92 | CMB (calibration) | +26.9% |

The particle horizon: base `10953.6 Mpc`; with the correction `13905.7 Mpc` — `32.9 Mpc` beyond the last-scattering surface: in the corrected geometry the CMB lies almost flush against the horizon.

The reach of the probes: DESI DR2 (`z ≤ 2.33`), DESI DR3 (~2027, `z ≤ 3.55`), DESI-II/LBG (~2029, `2.5 < z < 3.5`) — all below the threshold. The radial BAO in the window `2.33–4.53` is a trace of the core, not of `μ_d`; the threshold itself affects only `D_M` above `z ≈ 4.5`, where no precise rulers exist apart from the CMB calibration itself. A second independent anchor for `μ_d` is out of reach until a precise distance scale appears at `z ≳ 5–6`.

### 8.5 Consistency flag of the P6 formulation

The observed galaxies (JWST, `z ≈ 10–14`; `μ_e(z=14) = 0.9972`) lie **above** the detachment threshold (`μ_d = 0.9806`, `z_detach = 4.53`) — they exist before detachment. A literal reading — "no separate objects exist at `μ > μ_d`" — would contradict the observations; the P6 formulation adopted in the paper refers the threshold to the **average medium along the light path**: local clumps with maximal binding fall into themselves and collapse into objects before the threshold, without changing the average availability of space. The script prints this flag explicitly; the distance arithmetic and the results of Part I do not depend on it.

---

## 9. Part III: structure growth (`fσ₈`)

The DSIC core sets only the **background** — the expansion rate `H(z)`. A separate question is whether this background is compatible with the observed **growth of inhomogeneities** (galaxies, clusters). Block `[11]` takes a minimal, honestly bounded step: on the DSIC background, the **standard** equation of linear growth of matter perturbations in the sub-horizon GR limit is solved,

```
d²δ/d(ln a)² + [2 + dln H/dln a]·dδ/d(ln a) − (3/2)·Ω_m(a)·δ = 0,
```

and the observable combination `fσ₈(z) = f·σ₈·δ(z)/δ(0)`, `f = dln δ/dln a`, is computed. The comparison is against the **Gold-2017** compilation (Nesseris, Pantazis & Perivolaropoulos 2017; 18 weakly correlated points, the covariance of the three WiggleZ points included).

> **Status — a survival test, not a derivation from the core.** DSIC does not derive the growth equation from first principles: the core introduces no densities, so `Ω_m0` enters here as a growth parameter, and the equation itself is taken in the GR limit on the DSIC geometry. The result removes the objection of gross incompatibility but is not a prediction from the core.

### 9.1 The result and the `(Ω_m0, σ₈)` degeneracy

The `fσ₈` data constrain not `Ω_m0` and `σ₈` separately but their combination `S8 = σ₈·√(Ω_m0/0.3)`. A scan along the degeneracy valley (`μ₀ = 0.76` fixed from the background):

| `Ω_m0` | `σ₈*` | `χ²/dof` | `S8 = σ₈·√(Ω_m0/0.3)` |
|---|---|---|---|
| 0.15 | 1.040 | 0.72 | 0.735 |
| 0.20 | 0.925 | 0.70 | 0.755 |
| 0.25 | 0.835 | 0.73 | 0.762 |
| 0.30 | 0.770 | 0.78 | 0.770 |
| 0.35 | 0.715 | 0.84 | 0.772 |
| 0.40 | 0.670 | 0.91 | 0.774 |

```
Formal minimum:  Ω_m0 = 0.185,  σ₈ = 0.955,  χ²/dof = 0.70  (18 points)
Invariant along the valley:  S8 = 0.761 ± 0.013
Slice at Planck-like Ω_m0 = 0.30:  σ₈ = 0.770,  χ²/dof = 0.73,  S8 = 0.770
```

Along the entire valley, `χ²/dof ≈ 0.7–0.9` — the DSIC background describes structure growth no worse than ΛCDM. The invariant `S8 = 0.761 ± 0.013` is stable and agrees with independent estimates (Planck `S8 ≈ 0.83`, weak lensing `S8 ≈ 0.76`); the DSIC background lands in the observed range without any stretching.

### 9.2 The `fσ₈(z)` curve (the slice `Ω_m0 = 0.30`)

The DSIC `fσ₈(z)` curve is shallow, peaking around `z ≈ 0.4–0.6`, and passes through the data cloud:

| `z` | `fσ₈` model | `fσ₈` obs. ± error |
|---|---|---|
| 0.02 | 0.399 | 0.428 ± 0.046 |
| 0.15 | 0.421 | 0.490 ± 0.145 |
| 0.38 | 0.440 | 0.440 ± 0.060 |
| 0.59 | 0.440 | 0.488 ± 0.060 |
| 0.73 | 0.434 | 0.437 ± 0.072 |
| 0.86 | 0.425 | 0.400 ± 0.110 |
| 1.40 | 0.375 | 0.482 ± 0.116 |

*(the table is thinned; the full output for all 18 points is in the log of block `[11]`).*

> **Simplifications (stated honestly).** (1) The Alcock–Paczynski correction for the difference between the DSIC background distances and the fiducial cosmology of each RSD measurement was not applied; at `z < 1` the AP factor is typically `< 2–3%`, small against the errors. (2) The growth equation is the GR-limit one; deriving growth from binding (§12 of the paper) remains a direction for development. Both simplifications are flagged as sources of systematics and do not affect the qualitative conclusion (compatibility).

---

## 10. Bottom line

**What has been checked and holds:**

- the internal identities of the core hold to machine precision (`≤ 5.7·10⁻⁸`);
- supernovae (anchored fit): Δχ² = +0.05 — indistinguishable from ΛCDM; the parameter is measured: `μ₀ = 0.7600 ± 0.0091`;
- pipeline validation: the ΛCDM branch reproduces the published Pantheon+SH0ES values;
- robustness to the treatment of errors (STAT-ONLY): Δχ² = −0.39, the shift of `μ₀` within 1σ;
- SN + BAO joint: Δχ² = −5.0 (a weak formal preference for DSIC); BAO `χ²/dof` = 1.04 versus 1.36;
- the ruler `r_d ≈ 148 Mpc` on the early scale, constant in `z`;
- the second tier: the closure to `θ*` is exact, the late Universe is identically untouched, the horizon is finite and lies `33 Mpc` beyond the CMB;
- structure growth: the DSIC background describes `fσ₈` (Gold-2017) no worse than ΛCDM (`χ²/dof ≈ 0.7`), the invariant `S8 = 0.761 ± 0.013` is within the observed range.

**What is not confirmed and is stated honestly:**

- the `Om(z)` test with current data does not discriminate between the models — the sign of Δχ² depends on `H₀` (in both covariance modes);
- `μ_d` is a calibration against a single observation (0 degrees of freedom) with the form fixed by the postulate; `r*` is borrowed from ΛCDM;
- the P6 formulation requires the "average medium" reading (the flag of §8.5);
- structure growth is a survival test, not a derivation from the core: the growth equation is taken in the GR limit, and `Ω_m0` enters as a growth parameter (§9).

**Where the model will meet its decisive tests:**

- the shape of `Om(z)`; Lyman-α BAO at `z ≈ 3–3.5` (up to ~5%, within reach of DESI DR3);
- **cosmic time at `z = 5–12`**: at `z = 7.5`, 0.42 Gyr are available versus 0.71 in ΛCDM — the Salpeter budget of black hole growth shrinks from ~16 to ~9 e-foldings; this is the core's narrowest point, stated as a number.

DSIC is a viable one-parameter geometric alternative on late-time data, with sharp falsifiable predictions on early ones. Current data neither reject it nor confirm it more strongly than ΛCDM; the decisive factors will be high `z`, an independently pinned `H₀`, and the cosmic time budget.

---

## 11. Reproducibility

All results are produced by a **single script**, `en_dsic_test.py`; it runs as a single cell in Colab or via `python en_dsic_test.py`. The script downloads the public catalogs itself and computes χ² with the full covariance; the provenance of every hard-coded figure is documented in the file itself.

> **On run time.** A full run takes roughly 3–5 minutes: downloading the supernova catalog (~65 MB), the Pantheon+SH0ES fits (a few minutes), and block `[11]` (structure growth) — the longest stage (a grid scan with numerical ODE integration at every node, ~2 minutes). Between the printout of block `[10]` and the results of `[11]` there may be a long silent pause — this is expected; the script is computing, no need to interrupt it.

Data: **Pantheon+SH0ES** (`Pantheon+SH0ES.dat`, `STAT+SYS.cov`, `STATONLY.cov`) from the `PantheonPlusSH0ES/DataRelease` repository (auto-download); **DESI DR2** BAO (Table IV, 6 bins, hard-coded, self-check `[0]`); cosmic chronometers — a compilation of 33 points with the full Moresco et al. (2020) covariance for the 15 Moresco points; **Planck 2018** (`z*`, `r*`, `θ*`, `d_LSS`, self-check `[7]`); **Gold-2017** `fσ₈` (18 points, Nesseris et al. 2017; from Table I of Sagredo et al. 2018; the WiggleZ covariance — Blake et al. 2012), hard-coded into block `[11]`.

| Block | What it does | Section |
|---|---|---|
| `[0]` | Self-check of the hard-coded DESI against Table IV | §11 |
| `[0b]` | Machine verification of the core's internal identities | §2 |
| `[1]` | Anchored SN fit (STAT+SYS), explicit `M`, 1σ from the Hessian | §3.1 |
| `[2]` | STAT-ONLY control | §3.2 |
| `[3]` | Joint SN+BAO, free `r_d`, the `r_d` ruler | §4 |
| `[4]` | `Om(z)`: robustness to `H₀`, DIAG + FULL covariances | §5 |
| `[5]` | Falsifiable predictions at high `z` | §6 |
| `[6]` | Cosmic time `t(z)`, DSIC vs ΛCDM | §7 |
| `[7]` | Self-check of the Planck reference data | §8.1 |
| `[8]` | The core's shortfall, fixing `μ_d` via `θ*` | §8.2 |
| `[9]` | Zero-below-threshold control; the `A ↔ μ_d` degeneracy | §8.3 |
| `[10]` | The trace of `μ_d`, the horizon, reach, the P6 flag | §8.4–8.5 |
| `[11]` | Structure growth `fσ₈` on the DSIC background (GR limit), the `S8` invariant | §9 |

Dependencies: `numpy`, `scipy` (pre-installed in Colab).
