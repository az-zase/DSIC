# DSIC

**Dual Scale Inversion Cosmology (DSIC)** is a single-parameter cosmological model. 

Space and matter as two sides of one flow. 

The model describes the Universe as a continuous transfer between two complementary scales:
- Space (S)
- Objectness (O)

Rather than introducing dark energy, DSIC attributes the observed late-time acceleration to the geometry of this dual-scale evolution. The entire late-Universe cosmology is generated from a single observable phase parameter μ, with the present epoch corresponding to μ₀ ≈ 0.76.

Current implementation reproduces the expansion history of the mid-to-late Universe using one free cosmological parameter and matches ΛCDM on the full Pantheon+SH0ES supernova sample (Δχ² ≈ +0.05), while also providing consistent fits to BAO, cosmic chronometers, the transition redshift, and the age of the Universe.

---

## Repository structure
```
README.md
article/
  dsic-azamat-zaseev.md
  dsic-cmb-plugin-azamat-zaseev.md
  dsic-formulas-azamat-zaseev.mmd
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

## article/dsic-cmb-plugin-azamat-zaseev.md

An optional extension to the core model.

This document collects ideas that are not derived from the DSIC core equations and therefore are intentionally kept separate.

It contains two independent phenomenological modules:
- **Part A:** the detachment parameter μ_d, introduced to extend the model toward the CMB distance while leaving the late-Universe fit unchanged
- **Part B:** a tentative gravitational space-redistribution module describing local anisotropic redistribution of space

These extensions are exploratory and are not considered part of the core DSIC model.

---

## article/dsic-formulas-azamat-zaseev.mmd

Mermaid source of the complete DSIC formula diagram.

The diagram presents the logical structure of the model from the action principle through the observable cosmological equations and indicates the current boundary of applicability.

---

## Current status

### Core DSIC
- one cosmological parameter (μ₀)
- analytic closed-form cosmology
- no dark energy
- late-Universe observational tests

### Extensions
- CMB detachment parameter (μ_d)
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
