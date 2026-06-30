# DSIC

**Dual Scale Inversion Cosmology (DSIC)** is a single-parameter cosmological model without dark energy. 

Space and matter as two sides of one flow. 

The model describes the Universe as a continuous transfer between two complementary scales:
- Space (S)
- Objectness (O)

Rather than introducing dark energy, DSIC attributes the observed late-time acceleration to the geometry of this dual-scale evolution. The entire late-Universe cosmology is generated from a single observable phase parameter μ, with the present epoch corresponding to μ₀ ≈ 0.76.

Current implementation reproduces the expansion history of the mid-to-late Universe using one free cosmological parameter and matches ΛCDM on the full Pantheon+SH0ES supernova sample (Δχ² ≈ +0.05), while also providing consistent fits to BAO, cosmic chronometers, the transition redshift, and the age of the Universe. On the joint SN+BAO fit DSIC is marginally preferred (Δχ² ≈ −5.0, a weak formal preference). All empirical claims are reproduced by self-contained scripts included in the repository.

---

## Repository structure
```
README.md
article/
  dsic-azamat-zaseev.md
  dsic-formulas-azamat-zaseev.mmd
  formulas/
    (25 formula images, .png)
  dsic_test/
    dsic_test.py
    dsic-test-results.md
    dsic_graph.png
  dsic_cmb_plugin/
    dsic-cmb-plugin-azamat-zaseev.md
    cmb_plugin_test/
      dsic_cmb_test.py
      dsic-cmb-test-results.md
      dsic_cmb_graph.png
```

---

## article/dsic-azamat-zaseev.md

The main DSIC paper.

Contains the complete derivation of the core model:
- dual-scale geometry (S ↔ O)
- observable projection (μ)
- scale factor
- Hubble expansion
- luminosity distance
- cosmological observables
- comparison with observational data

The core model is intended for the late Universe (approximately z < 2.3) and is internally self-contained.

---

## article/dsic-formulas-azamat-zaseev.mmd

Mermaid source of the complete DSIC formula diagram.

The diagram presents the logical structure of the model from the action principle through the observable cosmological equations and indicates the current boundary of applicability.

The rendered formula images are in `article/formulas/`.

---

## article/dsic_test/

Empirical test of the core model and its reproducibility.

- `dsic_test.py` — a single self-contained script. It downloads the public Pantheon+SH0ES catalog automatically, has DESI DR2 BAO and the cosmic-chronometer compilation built in (with a self-check against the official tables), and computes χ² with the full covariance for both DSIC and flat ΛCDM at an equal number of free parameters.
- `dsic-test-results.md` — the write-up of the results: supernovae (anchored fit, Δχ² = +0.05), STAT-ONLY robustness, the joint SN+BAO fit (Δχ² = −5.0), the BAO `D_M`–`D_H` fit, the honest `Om(z)` outcome (does not resolve the models — the sign of Δχ² depends on the `H₀` normalization), and the falsifiable high-z predictions.
- `dsic_graph.png` — summary figure.

---

## article/cmb_plugin/dsic-cmb-plugin-azamat-zaseev.md

An optional extension to the core model.

This document collects ideas that are **not** derived from the DSIC core equations and are therefore intentionally kept separate.

It contains two independent phenomenological modules:
- **Part A:** the detachment parameter μ_d, introduced to extend the model toward the CMB distance while leaving the late-Universe fit unchanged
- **Part B:** a tentative gravitational space-redistribution module describing local anisotropic redistribution of space

These extensions are exploratory and are not considered part of the core DSIC model.

---

## article/cmb_plugin/cmb_plugin_test/

Separate reproducibility test for the CMB plugin (the detachment parameter μ_d).

- `dsic_cmb_test.py` — a single self-contained script (no external files; the Planck reference data are built in and self-checked).
- `dsic-cmb-test-results.md` — the write-up. It does **not** prove μ_d; it confirms that the plugin arithmetic is reproducible, that the Planck reference number is correct, that the late Universe (z ≤ z_detach ≈ 4.53) is untouched by the correction, and it states honestly that the parameter has zero degrees of freedom (the amplitude ↔ threshold degeneracy). The value `z_detach ≈ 4.53` is a calibration, not a prediction.
- `dsic_cmb_graph.png` — summary figure.

---

## Current status

### Core DSIC
- one cosmological parameter (μ₀)
- analytic closed-form cosmology
- no dark energy
- late-Universe observational tests (see `article/dsic_test/`)

### Extensions
- CMB detachment parameter (μ_d) — calibration only, zero degrees of freedom (see `article/cmb_plugin/`)
- phenomenological gravity module

These extensions remain separate until they can be independently derived or observationally validated.

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
