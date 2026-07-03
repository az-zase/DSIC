# DSIC

**Dual Scale Inversion Cosmology (DSIC)** is a single-parameter cosmological model without dark energy.

Space and matter as two sides of one flow.

The model describes the Universe as a continuous transfer between two complementary scales:
- Space (S)
- Objectness (O)

Rather than introducing dark energy, DSIC attributes the observed late-time acceleration to the geometry of this dual-scale evolution. The entire late-Universe cosmology is generated from a single observable phase parameter μ, with the present epoch corresponding to μ₀ ≈ 0.76.

The core implementation reproduces the expansion history of the mid-to-late Universe using one free cosmological parameter and matches ΛCDM on the full Pantheon+SH0ES supernova sample (Δχ² ≈ +0.05), while also providing consistent fits to BAO, cosmic chronometers, the transition redshift, and the age of the Universe. On the joint SN+BAO fit DSIC is marginally preferred (Δχ² ≈ −5.0, a weak formal preference). The paper additionally carries a **second floor** (the early Universe, via the detachment threshold μ_d) and a **structure-growth survival check** (fσ₈ on the DSIC background), both kept explicitly separate from the self-contained core. All empirical claims are reproduced by a single self-contained script included in the repository.

---

## Repository structure
```
README.md
article/
  DSIC - Azamat Zaseev.md
  dsic_formulas.mmd
  formulas/
    (25 formula images, .png)
  dsic_test/
    dsic_test.py
    dsic_results.md
    dsic_graph.png
```

---

## article/DSIC - Azamat Zaseev.md

The main DSIC paper — a single self-contained document that now merges the former core article and the CMB/gravity plugin into one text. It is organized in three parts, each with its status stated explicitly.

**Part I — the core** (late Universe, approximately z < 2.3; internally self-contained):
- dual-scale geometry (S ↔ O)
- observable projection (μ)
- scale factor
- Hubble expansion
- luminosity distance
- cosmological observables
- comparison with observational data, including a translation of the geometry into the standard Friedmann language (an effective equation of state w_eff(z))

**Part II — the second floor** (the early Universe; not derived from the core):
- the detachment threshold μ_d, introduced to extend the model toward the CMB distance while leaving the late-Universe fit unchanged (calibration, zero degrees of freedom)
- a tentative gravitational space-redistribution module (kept in an appendix as a development sketch, not a result)

**Part III — structure growth** (a survival check, not a derivation from the core):
- linear growth of density perturbations on the DSIC background in the GR limit, compared with the fσ₈ compilation

The boundary of applicability is stated throughout: the core is the accomplished part; the second floor and the growth check are flagged as exploratory and are not claimed as consequences of the core equations.

---

## article/dsic_formulas.mmd

Mermaid source of the complete DSIC formula diagram.

The diagram presents the logical structure of the model from the action principle through the observable cosmological equations and indicates the current boundary of applicability.

The rendered formula images are in `article/formulas/`.

---

## article/dsic_test/

Empirical test of the model and its reproducibility — now a single unified script and a single write-up covering all three parts.

- `dsic_test.py` — one self-contained script (blocks `[0]`–`[11]`). It downloads the public Pantheon+SH0ES catalog automatically; has DESI DR2 BAO and the cosmic-chronometer compilation built in (with a self-check against the official tables); carries the Planck 2018 reference data for the μ_d calibration (self-checked); and includes the Gold-2017 fσ₈ compilation for the growth test. It computes χ² with the full covariance for both DSIC and flat ΛCDM at an equal number of free parameters. Part II (the μ_d threshold) and Part III (structure growth) run in the same pass.
- `dsic_results.md` — the write-up of the results: supernovae (anchored fit, Δχ² = +0.05), STAT-ONLY robustness, the joint SN+BAO fit (Δχ² = −5.0), the BAO `D_M`–`D_H` fit, the honest `Om(z)` outcome (does not resolve the models — the sign of Δχ² depends on the `H₀` normalization), the falsifiable high-z predictions, the cosmic-time risk zone (z = 5–12), the second-floor μ_d calibration, and the structure-growth survival check (S8 invariant).
- `dsic_graph.png` — nine-panel summary figure.

> **Note on runtime.** A full run takes roughly 15–25 minutes: the supernova download (~65 MB), the Pantheon+SH0ES fits (a few minutes), and the growth block `[11]` (a grid scan with numerical ODE integration — the longest step). A long quiet pause between blocks `[10]` and `[11]` is expected; the script is working, not stuck.

---

## Current status

### Core DSIC
- one cosmological parameter (μ₀)
- analytic closed-form cosmology
- no dark energy
- late-Universe observational tests (see `article/dsic_test/`)

### Second floor (in the paper, not derived from the core)
- CMB detachment threshold (μ_d) — calibration only, zero degrees of freedom
- phenomenological gravity module (development sketch, kept in an appendix)

### Structure growth
- survival check on the DSIC background (GR-limit growth vs fσ₈, Gold-2017); the data constrain the invariant S8 rather than Ω_m0 and σ₈ separately — not a derivation from the core

The second floor and the growth check remain flagged as separate from the core until they can be independently derived or observationally validated.

---

## Citation

If you use the DSIC (Dual Scale Inversion Cosmology) model, its equations, derivations, figures, or ideas in academic work, presentations, articles, software, or derivative research, please cite the original work and acknowledge its author.

### Suggested citation

**Azamat Zaseev**  
*Dual Scale Inversion Cosmology (DSIC)*  
https://github.com/az-zase/DSIC.git

This repository is the original source of the DSIC model and is maintained by its author.

---

## A personal note

If you've reached this page, thank you for taking the time.

I am not a physicist or a member of any scientific institution. I am just a guy who asked questions and searched for meanings. DSIC is the result of that path.

If these ideas help someone see things differently, or become a small step toward something greater — then it was perfect.

I hope this work proves useful, even at the very edge of the Universe…

…which does not exist.
