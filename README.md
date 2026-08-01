# DSIC

**Dual Scale Inversion Cosmology (DSIC)** is a single-parameter cosmological model without dark energy.

Space and matter as two sides of one flow.

The model treats space and objectness as two poles of a single conserved quantity moving along the circle `S² + O² = R²` under one uniform flow. An internal observer, itself built of objectness, accesses only the modulus projection; its rod and clock uniquely fix the observable scale factor `b(μ) = √(1−μ²)/μ²`, from which redshift, kinematics and distances follow. Late-time acceleration emerges as a kinematic consequence of the divergence of the two scales — no dark energy is introduced. The entire late-Universe cosmology is generated from a single observable phase parameter, with the present epoch at `μ₀ = 0.7600 ± 0.0091`.

The shape of `b(μ)` has no functional freedom, so the parameter measured at low redshift forcibly determines the high-`z` behaviour: the divergences from ΛCDM at `z > 3` are predictions, not tuning.

---

## Repository structure

```
README.md
en DSIC 2 - Azamat Zaseev.pdf
en DSIC 2 - Azamat Zaseev.md
ru DSIC 2 - Azamat Zaseev.md
dsic2_formulas.mmd
dsic2_test.py                    <- current validation / reproducibility script (v2)
dsic2_test_results.txt           <- full output of a reference run of dsic2_test.py
```

---

## DSIC 2 — Azamat Zaseev.md (current paper)

A single self-contained document, organized in two tiers plus a growth check:

**Part I — the core (late Universe, z < 2.3).** Postulates; the dual-scale circle `S ↔ O`; the spinor (Möbius) form of the core `ψ = R·e^{iφ/2}` — the sign of objectness as the second sheet of a two-fold covering forced by the postulates; the modulus projection `μ`; the observer's rod and clock; the scale factor, redshift, kinematics with the built-in deceleration→acceleration crossover at `μ_cross = √(3−√5)` (`z_t = 0.769`); closed-form distances and a finite particle horizon; the mapping to the Friedmann form via `w_eff(z)`; machine verification of the internal identities; the confrontation with data; falsifiable high-`z` predictions.

**Part II — the second tier (early Universe).** One new entity (binding of objects), one postulate (P6) and one measured constant — the detachment horizon `μ_d = 0.9806` (`z ≈ 4.53`), calibrated against the directly measured angle `θ*` (calibration, not a prediction). The correction is identically zero below the threshold, so the late Universe is untouched. The same threshold sets the domain of the clock factor: clocks are normal inside the connected monolith and relative after detachment. The first objects are patches of the monolith that never uncoupled — expected early islands, not a contradiction with JWST galaxies at `z ≈ 10–14`. Island time makes the budget of every early object individually testable.

**Part III — structure growth.** Linear perturbation growth on the DSIC background in the GR limit against the Gold-2017 `fσ₈` compilation; the data constrain the invariant `S₈ = 0.761 ± 0.013`. A survival test, not a derivation from the core.

The status of every claim is separated explicitly: the core is self-contained; the extensions are not derived from it.

### Headline results (all reproduced by `dsic2_test.py`)

| probe | data | DSIC vs ΛCDM |
|---|---|---|
| supernovae | Pantheon+SH0ES, 1657 SN, STAT+SYS | indistinguishable, `Δχ² = +0.05` |
| parameter | χ² Hessian | `μ₀ = 0.7600 ± 0.0091` (~1.2 %) |
| joint SN+BAO | + DESI DR2, full 12×12 cov. | DSIC mildly preferred, `Δχ² = −5.0` |
| BAO `D_M`–`D_H` | DESI DR2 | `χ²/dof` 1.04 vs 1.36 |
| `Om(z)` | chronometers + DESI | does **not** discriminate (sign set by `H₀`) |
| high `z` | prediction | `D_H/r_d` diverges up to −5 % by `z = 3.5` |
| CMB distance | closure to `θ*` | exact by construction (`μ_d` calibrated) |
| cosmic time | detachment horizon | `t₀ = 13.5–13.6` Gyr (stellar anchor) |
| SMBH test | 6 early black holes | 5/6 accommodated; UHZ1 → heavy seed, same as its discoverers conclude |
| structure growth | Gold-2017, 18 points | `χ²/dof` 0.73 vs 0.74; `S₈ = 0.761 ± 0.013` |

In every fit DSIC and ΛCDM carry the same number of free parameters, so `Δχ² = ΔAIC = ΔBIC`.

---

## dsic2_test.py (current script)

One self-contained empirical-validation pipeline (`numpy`, `scipy`; runtime 10–15 minutes; the Pantheon+SH0ES catalog is downloaded automatically, ~65 MB). All hard-coded figures are documented in place and verified against the originals:

- **Pantheon+SH0ES** supernovae (auto-download, full STAT+SYS covariance);
- **DESI DR2 BAO** (Table IV, with a built-in self-check against the official values);
- **cosmic chronometers** (33 points; full Moresco et al. 2020 covariance for the 15-point subset);
- **Planck 2018** reference calibration (column `base_plikHM_TTTEEE_lowl_lowE_lensing`, Planck Legacy Archive tables);
- **early SMBHs** with spectroscopic redshifts and masses verified against the discovery papers (Bogdán 2024; Natarajan 2024; Goulding 2023; Maiolino 2024; Larson 2023; Wang 2021; Bañados 2018; Mortlock 2011);
- **Gold-2017 `fσ₈`** growth compilation (WiggleZ covariance included).

Blocks `[0]–[11]`: DESI self-check; machine verification of the core identities; the spinor (Möbius) core — 8 identities plus the forcedness test of the two-sheet covering; anchored SN fit with Hessian errors; STAT-ONLY robustness; joint SN+BAO with free `r_d`; the `Om(z)` diagnostic in two covariance modes; high-`z` predictions; piecewise cosmic time across the detachment horizon; the island-time test on real SMBHs; the Planck self-check; fixing `μ_d` via `θ*`; the zero-below-threshold control and the amplitude↔threshold degeneracy; the trace of `μ_d` and probe reach; structure growth and `S₈`.

`dsic2_test_results.txt` is the complete console output of a reference run — every number quoted in the paper can be located there by its block tag.

---

## Citation

If you use the DSIC (Dual Scale Inversion Cosmology) model, its equations, derivations, figures, or ideas, please cite the original work.

Azamat Zaseev
Dual Scale Inversion Cosmology (DSIC)
https://github.com/az-zase/DSIC.git

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21198493.svg)](https://doi.org/10.5281/zenodo.21198493)

---

## A personal note

If you've reached this page, thank you for taking the time.

I am not a physicist or affiliated with a scientific institution. I am just a guy who asked questions and searched for meanings. DSIC is the result of that path.

If these ideas help someone see things differently, or become a small step toward something greater — then it was perfect.

I hope this work proves useful, even at the very edge of the Universe...

...which does not exist.
