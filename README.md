# DSIC

**Dual Scale Inversion Cosmology (DSIC)** is a single-parameter cosmological model without dark energy.

Space and matter as two sides of one flow.

The model treats space and matter as complementary manifestations of a single continuous flow.

Rather than introducing dark energy, DSIC interprets the observed late-time acceleration through the geometry of this dual-scale evolution. The entire late-Universe cosmology is generated from a single observable phase parameter μ, with the present epoch corresponding to μ₀ ≈ 0.76.

The core implementation reproduces the expansion history of the mid-to-late Universe using one free cosmological parameter. It matches ΛCDM on the full Pantheon+SH0ES supernova sample (Δχ² ≈ +0.05), while also providing consistent fits to BAO, cosmic chronometers, the transition redshift, and the age of the Universe. On the joint SN+BAO fit DSIC is marginally preferred (Δχ² ≈ −5.0, a weak formal preference).

The paper additionally carries a **second floor** (the early Universe, via the detachment threshold μ_d) and a **structure-growth survival check** (fσ₈ on the DSIC background). Both are kept explicitly separate from the self-contained core. All empirical claims are reproduced by a single self-contained script included in the 

## Repository structure

```
README.md
article/
  DSIC - Azamat Zaseev.md
  dsic_formulas.mmd
  dsic_test.py
  dsic_results.md
  dsic_graph.png
  EN/
    en DSIC - Azamat Zaseev.md
    en dsic_formulas.mmd
    en_dsic_test.py
    en dsic_results.md
    en dsic_graph.png
```
---

## article/DSIC - Azamat Zaseev.md

The main DSIC paper is a single self-contained document.

It is organized in three parts:

Part I — core model (late Universe, z < 2.3):
- dual-scale geometry (S ↔ O)
- observable projection (μ)
- scale factor
- Hubble expansion
- luminosity distance
- cosmological observables
- mapping to standard Friedmann form via an effective equation of state w_eff(z)

Part II — second floor (early Universe):
- detachment threshold μ_d extending the model toward the CMB distance (calibration only)
- gravitational space-redistribution module (appendix, exploratory)

Part III — structure growth:
- linear growth of density perturbations on the DSIC background in the GR limit
- comparison with the fσ₈ compilation

Applicability of each part is explicitly separated: the core is self-contained; extensions are not derived from it.

---

## article/EN/

Complete English-language mirror of the repository.

This directory contains English versions of:
- main paper
- formula diagram (Mermaid)
- reproducibility script
- results summary
- figures

All files correspond one-to-one with the original version.

---

## article/dsic_formulas.mmd

Mermaid source of the DSIC formula diagram.

The diagram describes the structure of the model from the action principle through observable cosmological relations and marks the boundary of applicability.

The rendered output is stored as article/dsic_graph.png.

---

## article/dsic_test.py

Empirical validation and reproducibility pipeline.

The script includes:
- Pantheon+SH0ES supernova compilation
- DESI DR2 BAO dataset
- cosmic chronometers compilation
- Planck 2018 reference calibration
- Gold-2017 fσ₈ growth dataset

It computes χ² for DSIC and ΛCDM using identical parameter counts.

Includes:
- late-time expansion fit
- BAO consistency checks
- joint SN+BAO analysis
- high-redshift extrapolation tests
- structure growth evolution

Runtime: 3–5 minutes.

---

## article/dsic_results.md

Summary of numerical results:
- Supernova fit: Δχ² ≈ +0.05
- Joint SN+BAO: Δχ² ≈ −5.0
- BAO consistency
- Cosmic chronometers agreement
- Transition redshift consistency
- High-redshift behavior (z ≈ 5–12)
- Structure growth consistency

---

## article/dsic_graph.png

Nine-panel visualization of cosmological fits and comparisons.

---

## Core model

- one cosmological parameter (μ₀)
- analytic closed-form evolution
- no dark energy
- direct observational comparison

---

## Second floor

- CMB extension via detachment threshold μ_d
- calibration-only parameterization
- exploratory gravitational module

---

## Structure growth

- linear perturbation evolution on DSIC background
- comparison with fσ₈ data
- invariant constraint on S₈

---

## Citation

If you use the DSIC (Dual Scale Inversion Cosmology) model, its equations, derivations, figures, or ideas, please cite the original work.

Azamat Zaseev  
Dual Scale Inversion Cosmology (DSIC)  
https://github.com/az-zase/DSIC.git

[![DOI](https://zenodo.org/badge/1281640542.svg)](https://doi.org/10.5281/zenodo.21198493)

---

## A personal note

If you've reached this page, thank you for taking the time.

I am not a physicist or affiliated with a scientific institution. I am just a guy who asked questions and searched for meanings. DSIC is the result of that path.

If these ideas help someone see things differently, or become a small step toward something greater — then it was perfect.

I hope this work proves useful, even at the very edge of the Universe...

...which does not exist.
