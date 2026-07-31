# Dual Scale Inversion Cosmology (DSIC): a one-parameter background cosmology from the relativity of scales

**A. Zaseev**

Independent Researcher, Russia
ORCID: 0009-0004-9345-1069
DOI: https://doi.org/10.5281/zenodo.21330687
Code and data: https://github.com/az-zase/DSIC

---

## Abstract

A geometric model of the cosmological background is presented in which space and objectness are two poles of a single conserved quantity and no absolute scale exists: only the ratio of the poles carries physical content. The state of the system moves along the circle `S² + O² = R²` driven by a single uniform flow; an internal observer, itself built of objectness, has access only to the modulus projection. The observer's ruler and clock uniquely fix the observable scale factor `b(μ) = √(1−μ²)/μ²`, from which redshift, kinematics and distances follow. The acceleration of the expansion turns out to be a kinematic consequence of the divergence of the two scales; dark energy and analogues of `Ω_m, Ω_Λ` are not introduced.

The shape of `b(μ)` is fixed by the construction and contains no functional freedom, so the single parameter `μ₀ = 0.7600 ± 0.0091`, measured with late-time data, uniquely determines the high-`z` behaviour: the divergences from ΛCDM at `z > 3` are forced predictions, not the result of tuning. On the Pantheon+SH0ES catalog (1657 SNe, STAT+SYS covariance) the model is statistically indistinguishable from flat ΛCDM (`Δχ² = +0.05` at an equal number of parameters); the joint SN+BAO fit formally mildly prefers DSIC (`Δχ² = −5.0`); the `Om(z)` diagnostic does not discriminate between the models on current data. The early Universe is placed on a second tier: the binding postulate fixes the form of the correction to the light path, and the single threshold `μ_d = 0.9806` (`z ≈ 4.53`) is calibrated against the angle `θ*`; its status is a calibration under a postulated form, not a prediction. Cosmic time is tied to the same threshold, read as a detachment horizon: clocks are normal in the connected monolith and relative after detachment, which turns the time budget of early objects into an individually testable requirement — 5 of the 6 observed supermassive black holes are accommodated at reasonable assembly epochs. The background of the model is compatible with standard (GR-limit) growth of linear perturbations: `S₈ = 0.761 ± 0.013`. All empirical claims are reproduced by one open script.

**Keywords:** cosmology; dark energy; accelerated expansion; alternative background models; type Ia supernovae; baryon acoustic oscillations.

---

# 1. Introduction

## 1.1 Statement of the problem

The transition of the cosmological expansion from deceleration to acceleration at `z_t ≈ 0.6–0.8` [1, 2] is a firmly established fact. The standard model describes it with a component of negative pressure: two density parameters set the shape of the expansion curve, but the constant Λ itself remains a postulate without a derivation. The question is posed: can the same curve arise without an energy ingredient — from the geometry of what is called "scale"?

DSIC answers in the affirmative, in the form of an explicit construction. The starting premise: no absolute scale exists, neither for objects nor for space — only their ratio is observable. Space and objectness are treated as two poles of a single conserved norm, flowing along a circle under a single uniform flow. The observer sits inside the system and is itself made of objectness, so neither the sign of quantities nor their absolute values are accessible to it — only moduli and ratios. The requirement that the speed of light be constant in its units uniquely fixes the rate of its clock, and from the ruler and the clock follows

```
b(μ) = √(1 − μ²) / μ² ,
```

where `μ ∈ [0, 1]` is the observable fraction of objectness. This single function generates the redshift, the expansion rate, the deceleration parameter and all distances. The sign change of the deceleration parameter is built into the curvature of `b(μ)` and requires no Λ: acceleration arises as a kinematic consequence of the growing divergence of the two scale poles.

The key property of the model is not economy but **rigidity**: the shape of `b(μ)` contains no functional freedom, so the value of `μ₀`, measured with late-time data to ~1 % precision, uniquely determines the high-`z` behaviour with no room for tuning.

## 1.2 Organization of the paper and summary of results

The exposition is built up in order of increasing empirical risk and is split into two tiers of different status.

**Part I (the core, §2–§8)** is a closed construction: the postulates (§2), the geometry of the core (§3), the projection for the internal observer (§4), the observable physics (§5), machine verification of the internal identities (§6), the confrontation with late-Universe data at a number of free parameters equal to ΛCDM's (§7), and falsifiable high-`z` predictions (§8). The core is self-contained: none of its results depends on Part II.

**Part II (the second tier, §9–§11)** adds exactly one new entity (the binding of objects), one postulate and one measured constant — the detachment threshold `μ_d`. It closes the distance to the relic radiation (§10) and the cosmic-time budget of the early epochs (§11), leaving the late Universe identically untouched.

The paper then checks the compatibility of the background with structure growth (§12), outlines the limits of applicability (§13) and concludes (§14). All numerical results are produced by one open script (the "Reproducibility" section); the block tags `[0]–[11]` in the text refer to its output.

Summary of results:

| probe | data | metric | DSIC | ΛCDM | verdict |
|---|---|---|---|---|---|
| core consistency | identities | max deviation | `≤ 5.7·10⁻⁸` | — | hold |
| supernovae | Pantheon+SH0ES, 1657 SNe | `χ²/dof` | 0.8783 | 0.8783 | indistinguishable (`Δχ² = +0.05`) |
| parameter | `χ²` Hessian | 1σ | `μ₀ = 0.7600 ± 0.0091` | `Ω_m = 0.3308 ± 0.0180` | measured (~1.2 %) |
| robustness | STAT-ONLY | `Δχ²` | — | — | stable (`−0.39`) |
| SN + BAO | + DESI DR2 (12×12) | total `Δχ²` | — | — | DSIC by `−5.0` (weak) |
| BAO `D_M`–`D_H` | DESI DR2 | `χ²/dof` | 1.04 | 1.36 | no worse |
| `Om(z)` | chronometers + DESI | sign of `Δχ²` | — | — | **does not discriminate** |
| high `z` | prediction | `E(z)`, `D_H/r_d` | diverges | — | falsifiable |
| distance to the CMB | closure to `θ*` | `100·θ` | 1.04110 | — | exact (by construction) |
| cosmic time | detachment horizon | `t₀` | 13.5–13.6 Gyr | 13.95 | inside the stellar anchor |
| SMBH test | 6 early black holes | accommodation | 5/6 | — | falsifiable |
| structure growth | Gold-2017, 18 points | `χ²/dof` | 0.73 | 0.74 | no worse |
| invariant `S₈` | `(Ω_m0, σ₈)` degeneracy | 1σ | `0.761 ± 0.013` | 0.760 | consistent |

In every fit DSIC and ΛCDM carry the same number of free parameters, so `Δχ² = ΔAIC = ΔBIC`.

---
---

# Part I — The core and the observable cosmology

---

# 2. Postulates and the language of the model

## 2.1 Postulates

**P1 (relativity of scale).** No absolute scale exists. Only the ratio of scales carries physical content; any construction that opens access to an absolute scale is forbidden.

**P2 (two poles, one norm).** Objectness and spatiality are two poles of a single relation. An object is a local convexity: it has an interior and an outer boundary, and occupies a place relative to other objects. Space is the antipode of the object, a universal concavity: it is located nowhere and is contained in nothing; the category "where" arises only inside it. The scales of the poles obey a conservation law of the norm:

```
S² + O² = R² = const .
```

One pole can grow only at the expense of the other; the physics is carried by the fixedness of the norm during the flow between the poles, not by the absolute values of the poles.

**P3 (a single flow).** The system has one fundamental degree of freedom — a continuous monotonic flow `φ` that never reverses. Its dynamics is free: the action contains no potential, no forces, no cosmological constant,

```
𝒮[φ] = ∫ dt · ½ φ̇² ,     δ𝒮 = 0  ⇒  φ̈ = 0  ⇒  φ(t) = ω·t ,   ω > 0 .
```

Cyclicity, if it is observed, is born not of the flow but of the way it is observed.

**P4 (the internal observer).** Any observer is inside the system and is itself an object: its instruments transform together with it. Neither the sign of quantities nor the absolute scale is accessible to it; only moduli are.

**P5 (constancy of the speed of light).** In the internal observer's units `c = const`; together with P1 this uniquely fixes the rate of its clock (§5.1). This is an exchange rate fixed by convention, linking the ruler and the clock — one more ratio, not an absolute quantity; the numerical value of `c` does not enter the derivation of the geometry.

## 2.2 Layers of the exposition

So that the boundary between the mechanism and its observable shadow is never blurred, the exposition is split into layers; the layer of a symbol is visible from its form:

- **layer 0 — the flow**: the single monotonic parameter `φ` (time);
- **layer 1 — the core**: signed quantities `S, O, σ`; no forces, no cycles, no reversals; the modulus is not applied;
- **layer 2 — the projection**: unsigned quantities `P, Q, μ`, accessible to the observer; born by the operation `|·|`; cycles exist only here;
- **layer 3 — the observable physics**: `b, z, H_τ, q, χ, d_L`.

The boundary between the layers is absolute: no sign passes from the core into the projection, and no projection quantity enters the equations of the core.

## 2.3 Glossary

The model uses its own terminology; each term is given a standard analogue. The terms are defined rigorously in §3–§5; the table serves as a map read in advance.

| DSIC term | what it is | closest familiar image |
|---|---|---|
| objectness (`O`) | the pole of "thingness" | matter density / matter |
| space (`S`) | the pole of "extension" | spatial scale |
| norm `R` | `S² + O² = R²` | a conservation law |
| flow `φ` | the single monotonic degree of freedom | time as the evolution parameter |
| the Mush | a physical state of maximal objectness — the zone on both sides of the zero of space (`μ → 1`) | a regular analogue of the initial singularity |
| object bounce | the limit `μ → 0` | de Sitter asymptotics |
| eversion | a sign-phase change in the core | a mirror phase, unobservable from inside |
| bounce | the observable reversal at a limit | an artifact of the projection, not an event |
| binding (Part II) | a measure of the cohesion of objects | gravitational binding |
| detachment (`μ_d`) | the detachment horizon: the threshold at which the second pole seeps through for a fraction of the medium | the epoch after which the core geometry works without correction; the clock rate switches at the same threshold |

Two statuses must be kept apart. The observable **reversals** ("bounces") at the limits are not events: they are how the internal observer sees the smooth monotonic motion of the core through the modulus operation (§4); beneath the visible oscillation between the two limits lies one forward motion. **The Mush**, however, is a physical state: objectness at its maximum. It occupies the entire zone near the zero of space, on both its sides — both when space approaches zero, squeezing the objects into the Mush, and when it seeps out of zero, gradually releasing them; on the outgoing side the state of the Mush holds until monolithicity begins to break apart (§3.2, Part II).

---

# 3. The core: the circle, the inversions, the sign phases

## 3.1 The circle and the quadrature

Both poles are built from one smooth carrier, shifted by a quarter turn:

```
S(φ) = R·sin(φ/2) ,   O(φ) = R·cos(φ/2) ,   S² + O² = R² ;
S(−φ) = −S(φ) ,   O(−φ) = O(φ) ,
W = S·(dO/dφ) − O·(dS/dφ) = −R²/2 = const .
```

The carrier has a period of `4π` in the flow phase — a full sign cycle. The constancy of the Wronskian `W` is the formal record of the quadrature: when one pole is at an extremum, the other passes through zero, so every discontinuous derivative of a modulus (arising later, in the projection) is multiplied by a paired factor that vanishes.

In the `(S, O)` axes the invariant draws a circle of radius `R`: every admissible state is a point on it, and the evolution is a glide at constant angular velocity. Leaving the circle inward or outward is impossible; passing a point on an axis, the state neither stops nor bounces. The rate of the flow between the poles is nevertheless non-uniform, `dμ/dφ ∝ sin(φ/2)` — slow near the poles, fast in the middle of the arc. This non-linearity is a geometric property of projecting uniform rotation, not the action of a force.

## 3.2 The inversions: the Mush and the object bounce

Two events strictly alternate along the circle — smooth zero crossings of the poles:

```
inversion of space:    S(φ) = 0  ⇒  φ = 2kπ ,
inversion of objects:  O(φ) = 0  ⇒  φ = π + 2kπ .
```

The zero of `S` coincides with the maximum `|O| = R`: this is the state of **the Mush** — the limiting density of objectness. The Mush is the entire zone near the zero of space, on both its sides: the half before the zero is the finale of the contraction, when the vanishing space squeezes the objects together; the half after it is the start of the expansion, when the seeping space releases them. Only the second, outgoing half deserves the name **Big Bang**; on it the state of the Mush persists while the medium remains monolithic (its break-up is described in Part II). The zero of `O` coincides with the maximum `|S| = R`: the minimum of objectness — the **object bounce** in the observer's terms.

The order of alternation is strict: the Mush → the object bounce → the next Mush.

## 3.3 The spinor form: the sign as a second sheet

The sign mechanics of the core — eversion, inheritance, mirror phases — needs no separate rules: it folds into one complex object. The poles are combined into

```
ψ(φ) = O(φ) + i·S(φ) = R·e^{iφ/2} ,      O = Re ψ ,   S = Im ψ .
```

Then the conservation of the norm `|ψ|² = R²` and the quadrature `W = −R²/2` hold identically rather than being postulated separately (machine check — §6). The entire content of §3.1–3.2 is the motion of a single complex number along a circle at the constant angular velocity `φ̇/2`.

**The two-sheet covering (the Möbius-strip structure).** The phase of `ψ` is `φ/2`, hence

```
ψ(φ + 2π) = −ψ(φ) ,        ψ(φ + 4π) = +ψ(φ) .
```

One turn of the observable cycle (`φ: 0 → 2π`, over which `μ` runs a full period) flips `ψ` into its negative; the return to the initial state comes only after **two** turns (`4π`). This is exactly the property of a Möbius strip — the side flips over one circuit and is restored after two — and it is at the same time the spinor property of the `SU(2)` covering over `SO(3)`: a `360°` rotation gives a minus, `720°` restores. The observer, who has access only to the modulus `μ = |cos(φ/2)|` with period `2π` (§4), lives on one sheet; the core lives on the double one. **The minus phase of objectness is not a separate rule but the second sheet of this covering.**

**The asymmetry of the two zeros** falls out of which pole crosses the real/imaginary axis:

- At `O = 0` (`ψ` crosses the imaginary axis) only objectness passes through its own zero — it everts; space is at its maximum and performs merely an observable bounce, the sign of its phase unchanged.
- At `S = 0` (`ψ` crosses the real axis) space itself passes through its zero — the common container everts, and objectness, being at its maximum, **inherits** the sign inversion synchronously, including the observer's instruments.

The inherited sign of objectness is a choice of sheet by the sign of space, and it is written in closed form, without manual `σ`:

```
O_full(φ) = R·cos(φ/2)·sign[sin(φ/2)]     ( = Re ψ on the sheet set by the sign of Im ψ ) .
```

The quantity `O_full` is continuous at object inversions and discontinuous at space inversions (a sheet flip). The discontinuity is unobservable from inside: `|O_full| = |O|` is continuous everywhere, and the synchronous sign change of the system and the instruments cancels. It is precisely the continuity/discontinuity of `O_full` that determines, in §5.3, the observed asymmetry of the two scale limits (objectness is **subordinate** to space).

> **Status of the form (forcedness, not a choice).** The two-sheet covering is neither a postulate nor a fit but a consequence of two postulates already accepted. From P2 (a smooth core with a conserved norm, motion along a circle in the `(O,S)` plane) and P4 (only the modulus is accessible to the observer) it follows that the observable quantity `μ` identifies `ψ` and `−ψ`, so the sign cycle of the core is exactly twice the observable one: `ψ(φ+T_obs) = −ψ(φ)`, `ψ(φ+2T_obs) = +ψ(φ)`. This is the Möbius-strip (spinor `SU(2)`-covering) structure, and it is forced — verified by scanning the generalized core `ψ = R·e^{iαφ}`: the cycle ratio equals 2 for any `α > 0` (§6, block `[0c]`). The phase exponent (`α = 1/2` in `sin(φ/2), cos(φ/2)`) is a choice of normalization of the flow `φ` — an unobservable quantity: the redefinition `φ → φ/(2α)` maps any `α` into `1/2` without changing a single observable. The spinor form is therefore identically equivalent to the trigonometric core of §3.1 (eight identities, §6); the novelty is not in the physics but in the fact that the sign becomes the geometry of a sheet, forced by P2+P4, rather than a set of rules.

| zero | geometry of `ψ` | what changes its sign phase | kink for the observer |
|---|---|---|---|
| `O = 0` | crossing the imaginary axis | objectness itself | on `Q = \|O\|` |
| `S = 0` | crossing the real axis | space + objectness by inheritance (sheet flip) | on `P = \|S\|` |

Terminological convention: **eversion** is a sheet change at the level of the core; **bounce** is what the internal observer sees. These are descriptions of one mechanism from two vantage points.

## 3.4 Falling into itself

The universal concavity of space (P2) sets the character of the collapse of objectness: falling is not a motion toward a point in space but the convergence of an object's scale toward its own center, proceeding inward from everywhere at once, with no preferred direction. The gravitational development of this idea (a common center of fall for bound groups) is placed in Appendix B, as not derived from the core.

---

# 4. The projection: the world of the internal observer

An external observer (an imagined vantage point outside space; in reality there is none) sees the flow whole, with its sign: the history is continuous and monotonic, the quantities pass through their zeros and keep changing in the same direction. To the internal observer, by P4, only moduli are accessible: where the external observer sees a smooth zero crossing, the internal one sees a touch of zero and a reversal.

The transition from the core to the observable is performed by a single operation:

```
π : (S, O) ↦ (|S|, |O|) .
```

It cuts off the sign, folds the smooth zero crossing into a kink, and identifies adjacent sign branches, so that the observable period turns out to be half the sign period. The observable quantities:

```
P = |S| = R·|sin(φ/2)| ≥ 0 ,   Q = |O| = R·|cos(φ/2)| ≥ 0 ,   P² + Q² = R² ,
μ = Q/R = |cos(φ/2)| ∈ [0, 1] ,   P = R·√(1 − μ²) ,   Q = R·μ .
```

The pair `(P, Q)` always lies in the first quadrant. Over one observable cycle the point traverses this quadrant there and back: from the Mush (`P = 0, Q = R`) to the object bounce (`P = R, Q = 0`) and back. The back-and-forth motion arises from the symmetry of `Q = R·|cos(φ/2)|` about the bounce point, not from a superposition of arcs.

The carrier is smooth and has no vertices of its own, so both observable kinks of `μ` are of one nature — they are created by the modulus. The reversal requires no force, no braking, no cause: at the level of the mechanism there is no reversal at all, only a change in the way of seeing. What repeats is the projection, not the current: beneath it lies one monotonic forward motion in which nothing returns.

---

# 5. The observable physics

Everything in this section is built from `μ` and the internal observer's scales; it is an interpretation of the projection in physical units, not a separate mechanism.

## 5.1 The ruler, the clock, and the scale factor

The observer is made of objectness, so its unit of length is tied to the object pole: the ruler `L(μ) ∝ μ`. The metric is the ratio of the poles `P/Q = √(1 − μ²)/μ`.

The clock rate is fixed by postulate P5. Let the period of the proper clock be `T(μ) ∝ μᵖ`; then `c = L/T ∝ μ^{1−p}`, and the requirement `c = const` yields the unique `p = 1`:

```
dτ = dt / μ .
```

Any other exponent would drive `c` across the range and open access to the absolute scale forbidden by P1. The observable scale is composed of the metric and the clock rate:

```
b(μ) = (P/Q)·(1/μ) = √(1 − μ²) / μ² .
```

This is the central formula of the model: one function of one variable, with no free shape parameters.

**The domain of the clock factor.** The stretching `dτ = dt/μ` arises because time is measured *relative to* the second pole: an observer built of objectness reads the passage of time off the space available to it. The factor is therefore defined only where the second pole has already seeped through — below the detachment horizon (§9). In the connected monolithic objectness above the threshold there is nothing to measure against, and the clock rate there is normal (`dτ = dt`). For the entire late Universe, where the confrontations of §7 are performed, the factor acts in full force; the consequences for the early epochs are worked out in §11. This is the domain of a factor already introduced, not a separate postulate.

## 5.2 Redshift

Light propagates through space (the pole `P`) but is measured by an observer built of objectness. Propagation and measurement run along different poles, so both enter the observed wavelength, and each contribution enters exactly once:

```
(1) stretching of the path through space:  ∝ √(1 − μ²)    [pole P]
(2) the source clock, dτ = dt/μ:           ∝ 1/μ          [time]
(3) the object ruler, L ∝ μ:               ∝ 1/μ          [length]

λ_obs(μ) ∝ √(1 − μ²) / μ² = b(μ) .
```

The observed wavelength is identically equal to the scale factor, and the redshift drops out of the bookkeeping automatically:

```
1 + z = b(μ₀)/b(μ_e) = (μ_e²/μ₀²) · √(1 − μ₀²)/√(1 − μ_e²) ,
```

where `μ₀` is the phase of reception and `μ_e` the phase of emission. As `μ_e → 1`, `z → ∞`; as `μ_e → μ₀`, `z → 0`. The formula decomposes into the two poles:

```
1 + z = [ √(1 − μ₀²)/√(1 − μ_e²) ] · [ μ_e²/μ₀² ]
        └───── pole P (space) ─────┘  └ object scale ┘
```

The alternating dominance of these factors along the cycle is the very mechanism of the observed kinematics (§5.4). Near the Mush the pole `P` is squeezed almost to nothing — light has almost nothing to propagate through; the observed picture is already carried by the object scale.

The alternative "redshift by one ruler", `1 + z = μ₀/μ_e`, is rejected: it accounts only for contribution (3), whereas the photon lives through all three.

`b(μ)` is strictly decreasing on `(0, 1)`, so the correspondence `z ↔ μ_e` on the physical reception branch `μ_e ∈ (μ₀, 1)` is one-to-one and the inversion is well defined. Only modulus-valued quantities enter `1 + z` — the layer boundary is not violated.

## 5.3 The limits of the cycle and their asymmetry

```
μ → 0  (object bounce):  b → ∞ ;      μ → 1  (the Mush):  b → 0 .
```

At the object bounce space is maximal, objects are vanishingly small relative to it, and the observable scale diverges; this is a limit of approach, not an attainable singularity — the observer glides past `μ = 0` without reaching it. The `∞ vs 0` asymmetry follows from the sign mechanics (§3.3): objectness passes through its own zero by itself and smoothly (`O_full` continuous → divergence), while through the zero of space the container passes, everting the sign of its contents (`O_full` discontinuous → vanishing). The divergence of `b` is a direct expression of the relativity of scales: bounded forms of `b` do not satisfy postulate P1.

## 5.4 Kinematics

Introduce the branch direction `s = ±1` (`s = −1` — expansion) and the characteristic phase rate `k = (2/π)·ω > 0`, so that along a branch `dμ/dt = s·k = const`. All derivatives are taken with respect to proper time, `d/dτ = μ·d/dt`. The constancy of `k` is the definition of a branch, not an approximation: with a non-constant rate a term `∝ dk/dt` would enter the deceleration parameter, and here it is identically zero.

```
H_τ(μ) = (db/dτ)/b = s·k·(μ² − 2)/(1 − μ²) ,   |H_τ| = (2 − μ²)/(1 − μ²)  [units of k] ;
q(μ) = −(d²b/dτ²)·b/(db/dτ)² = (−μ⁴ + 6μ² − 4)/(μ⁴ − 4μ² + 4) .
```

The deceleration parameter depends on the direction `s` only through `s² = 1`, and not at all on the rate `k`: `q` is a function of the phase alone. The key property: `q(μ)` changes sign, with an exact crossover at `μ_cross = √(3−√5) = 0.874032`.

| state | `μ` | `q` | `\|H_τ\|` (units of `k`) | `w_eff` |
|---|---|---|---|---|
| the Mush | `→ 1` | `+1` | `→ ∞` | `+1/3` |
| transition | `0.874` | `0` | `≈ 5.2` | `−1/3` |
| today | `≈ 0.76` | `−0.43` | `≈ 3.4` | `−0.62` |
| midpoint | `0.5` | `−0.84` | `2.33` | `−0.89` |
| object bounce | `→ 0` | `−1` | `2` | `−1` |

The deceleration → acceleration transition arises from the curvature of `b(μ)`, without an introduced constant, and its mechanism can be read off the two-pole decomposition of the redshift. Near the Mush the space factor dominates — the deceleration regime (`q → +1`). In the late Universe the observer's object scale dominates — the acceleration regime, with de Sitter asymptotics as `μ → 0`. Acceleration here is not the result of something "pushing" the Universe: the object ruler and the scale of space diverge ever faster, and their relative divergence is read from inside as accelerated expansion.

The clock factor changes the kinematics substantially across the entire range, not only near the edges: at `μ = 0.5` the deceleration parameter without the clock equals `−1.27` versus `−0.84` with it.[^noclock]

[^noclock]: The formula without the clock factor: `q_noclock(μ) = (−2μ⁴ + 9μ² − 6)/(μ⁴ − 4μ² + 4)`.

## 5.5 Distances

Light travels with `c = const` in the local units of each epoch; the comoving distance accumulates as `dχ = (c/k)·[b(μ₀)/b(μ)]·dμ/μ`, and the integral is taken in closed form:

```
χ(μ_e; μ₀) = (c/k) · [√(1 − μ₀²)/μ₀²] · [ √(1 − μ₀²) − √(1 − μ_e²) ] ,
d_L(z) = (1 + z)·χ ,      DM(z) = 5·log₁₀(d_L / 10 pc) .
```

Two consequences. **The particle horizon is finite by itself**, with no cutoff: `χ(μ_e → 1) = (c/k)·(1 − μ₀²)/μ₀²`. **Self-consistency**: the identity `D_H = dχ/dz` holds exactly — the rate and the distances are tied by one geometry (machine check in §6).

## 5.6 Translation into the language of the Friedmann metric

The core can also be read in the standard language of GR without changing the formulas. The observable scale factor plays the role of the scale factor of the flat FRW metric, with the observer's proper time serving as cosmic time:

```
ds² = −c² dτ² + b(μ)² · dx² .
```

The question of the effective medium is then settled by the standard identity of flat FRW cosmology `w_eff = (2q − 1)/3`, which yields the closed expression

```
w_eff(μ) = (−μ⁴ + (16/3)μ² − 4) / (μ⁴ − 4μ² + 4) ,
```

whose values are listed in the table of §5.4. No new entity is introduced here: `w_eff` is entirely determined by the same `q(μ)`. In the language of GR, DSIC is a flat Universe with a single medium whose effective equation of state slides smoothly from radiation (`+1/3`) to a cosmological constant (`−1`). The model does not postulate a medium with such a `w(z)` — it obtains this curve from the geometry of the flow between scales; `w_eff` is merely its shadow in the language of GR.

> The coincidence of the background kinematics with some FRW model does not imply coincidence at the level of perturbations: a dynamical `w(z)` fixes the background but does not uniquely set the evolution of inhomogeneities (§12).

---

# 6. Internal consistency of the core: machine verification

Before any confrontation with data, it is checked that the observable physics of the core is internally closed — distances, rate and kinematics follow from one geometry without mismatches (script blocks `[0b]`, `[0c]`):

| identity | result |
|---|---|
| `b(μ)` strictly decreasing on `(μ₀, 1)`; the inversion `z → μ_e` well defined | OK |
| `D_H = dχ/dz` (comoving `D_H = c/H`) | max deviation `5.7·10⁻⁸` |
| closed-form `χ(μ)` = numerical `∫ c dz′/H(z′)` | max deviation `4.4·10⁻¹²` |
| `q(μ_cross) = 0` at the exact root `μ_cross = √(3−√5) = 0.874032` | OK; `z_t = 0.769` |
| age: `t₀·H₀ = ln(1/μ₀)·(2−μ₀²)/(1−μ₀²) = 0.9241` | `t₀ = 13.41 Gyr` at `H₀ = 67.4` |

**The spinor form `ψ = R·e^{iφ/2}` (§3.3)** — eight identities confirming the equivalence to the trigonometric core and the two-sheet covering structure:

| identity | deviation |
|---|---|
| `Re ψ = O = R cos(φ/2)` | `0` |
| `Im ψ = S = R sin(φ/2)` | `0` |
| `\|ψ\|² = R²` (norm conservation) | `4.4·10⁻¹⁶` |
| `W = S·O′ − O·S′ = −R²/2` (quadrature) | `1.1·10⁻¹⁶` |
| `O_full = R cos(φ/2)·sign[sin(φ/2)]` (inherited sign without `σ`) | `0` |
| `μ = \|cos(φ/2)\|` has period `2π` (the observable cycle) | `1.9·10⁻¹⁵` |
| `ψ(φ+2π) = −ψ(φ)` (the spinor minus, sheet flip) | `2.0·10⁻¹⁵` |
| `ψ(φ+4π) = +ψ(φ)` (return after two sheets) | `2.2·10⁻¹⁵` |

The forcedness of the covering is checked separately: for the generalized core `ψ = R·e^{iαφ}` the ratio of the sign cycle to the observable one equals 2 for every `α ∈ {0.25, 0.333, 0.5, 0.75, 1, 2}`.

The closed formulas and the numerical geometry coincide at machine precision: the expansion rate, the distances and the transition point are one geometry, not independent fits. The sign mechanics of §3.3 folds into a single complex number with no hand-set sign functions.

---

# 7. Confrontation with observations (`z < 2.3`)

The shape of every observable dependence is rigidly fixed by the function `b(μ)`; analogues of `Ω_m, Ω_Λ` are absent. The free quantities are the initial condition `μ₀` and the scale `H₀` (for the supernovae also the absolute magnitude `M`, as in ΛCDM).

**On scales.** BAO, the chronometers and the age are computed at `μ₀ = 0.76` on the early scale `H₀ = 67.4`; the geometry is converted into megaparsecs via `c/k = (c/H₀)·(2−μ₀²)/(1−μ₀²) = 14978 Mpc`. A direct fit to the supernovae gives `H₀ ≈ 73` — the familiar "Hubble tension" [8], reproduced identically by DSIC and ΛCDM. The model creates no new tension and removes no existing one.

## 7.1 Supernovae: Pantheon+SH0ES

The full catalog (1701 SNe; 1657 enter the fit, including 77 calibrators [5] as the anchor of the absolute scale), with the full STAT+SYS covariance [3, 4]. Three free parameters for each model (block `[1]`):

```
DSIC : μ₀ = 0.7600 ± 0.0091 ,  H₀ = 73.35 ± 1.01 ,  M = −19.2433 ± 0.0295 ,  χ²/dof = 0.8783
ΛCDM : Ω_m = 0.3308 ± 0.0180 , H₀ = 73.53 ± 1.02 ,  M = −19.2442 ± 0.0295 ,  χ²/dof = 0.8783

Δχ² (DSIC − ΛCDM) = +0.05        (ΔAIC = ΔBIC = +0.05)
```

The uncertainties are 1σ from the `χ²` Hessian at the minimum (`C = 2H⁻¹`); the model parameter is measured to ~1.2 %. On fully calibrated data DSIC is statistically indistinguishable from ΛCDM (`|Δχ²| ≪ 2`): a one-parameter geometry reproduces the same distance–redshift relation that the standard model sets with two density parameters.

**Pipeline validation.** The ΛCDM branch of the same pipeline reproduces the published values of the Pantheon+SH0ES analysis (`Ω_m ≈ 0.33 ± 0.018`, `H₀ ≈ 73.5 ± 1.0`); both models are computed by one code on the same data.

**Robustness.** A control run with the statistical-errors-only matrix (block `[2]`):

| | `μ₀` / `Ω_m` | `H₀` | `χ²/dof` |
|---|---|---|---|
| DSIC, STAT+SYS | 0.7600 | 73.35 | 0.8783 |
| ΛCDM, STAT+SYS | 0.3308 | 73.53 | 0.8783 |
| DSIC, STAT-ONLY | 0.7661 | 73.04 | 0.9242 |
| ΛCDM, STAT-ONLY | 0.3466 | 73.22 | 0.9244 |

`Δχ²` (STAT-ONLY) `= −0.39`; the shift of `μ₀` is within 1σ, the shape of the dependence is preserved. The growth of `χ²/dof` upon removing the systematics is the expected effect, identical for both models.

## 7.2 The baryon scale: DESI DR2

The DESI DR2 data (Table IV [6], `0.5 < z < 2.33`) give `D_M/r_d` and `D_H/r_d` separately in six bins; the values hard-coded in the script pass a self-check against the official table (block `[0]`). The dimensionless ratio `D_M/D_H` depends neither on `r_d` nor on `H₀` and tests the pure geometry of the expansion: `χ²/dof ≈ 0.78`, on a par with the standard model. The point `z = 0.934` gives the largest tension in ΛCDM as well (≈ 2.5σ) — a feature of the data, not of the model.

The scale `r_d` is recovered by an inverse problem, independently from the transverse and the radial measurement, on the early scale `H₀ = 67.4`:

```
r_d (from D_M) = 147.6 Mpc  (scatter 1.7) ;   r_d (from D_H) = 149.8 Mpc  (scatter 1.8) ,
```

consistent with the standard `≈ 147 Mpc` to within ~1.5 %. There is no dynamical `r_d(z)` in the data: the poles change with `z`, but the recovered physical scale is constant — the DSIC ruler is static in the observed window.

## 7.3 The expansion rate and the age

The cosmic chronometers [9] on `0.07 < z < 2` give `χ²/dof ≈ 0.37`. The integral over the expansion history gives `t₀·H₀ = 0.9241`, i.e. at `H₀ = 67.4` the age is `t₀ = 13.41 Gyr`.

A caveat on the comparison: the frequently quoted `13.80 ± 0.02 Gyr` is an age *derived within ΛCDM* from the Planck data, and measuring a non-ΛCDM model with it is partly circular; the model-independent stellar ages (the oldest globular clusters, `13.5 ± 0.5 Gyr`) pass this value. The quoted `t₀` is computed with the relative clock along the entire path; with the detachment horizon taken into account (§11), the age of an early-detached observer, `13.5–13.6 Gyr`, is the same quantity within the stellar anchor. The comparison of `H(z)` at `z < 2` (entirely below the threshold) is unaffected by the distinction.

## 7.4 The deceleration parameter and the transition redshift

Through the relation `μ(z)` at `μ₀ = 0.76`, the theoretical curve is compared with model-independent kinematic reconstructions [12]:

| `z` | `q_DSIC` | reconstruction |
|---|---|---|
| `0` | `−0.43` | `−0.55` |
| `0.3` | `−0.25` | `−0.34` |
| `0.65` | `−0.06` | `0.0` |
| `1.0` | `+0.11` | `+0.18` |
| `1.5` | `+0.30` | `+0.33` |
| `2.0` | `+0.44` | `+0.42` |

The shape of `q(z)` is S-shaped, reproducing the reconstructions within their scatter; the transition redshift `z_t = 0.769` agrees with the observed `0.6–0.8` [10].

## 7.5 The joint SN + BAO fit

Supernovae (the anchor) and the DESI DR2 BAO (6 bins) with the full 12×12 `D_M`–`D_H` covariance; `r_d` is a free parameter. Four free parameters for each model (block `[3]`):

```
DSIC : μ₀ = 0.7610 , H₀ = 73.33 , M = −19.243 , r_d = 136.5 Mpc , total = 1461.0
ΛCDM : Ω_m = 0.3029 , H₀ = 73.75 , M = −19.247 , r_d = 137.1 Mpc , total = 1466.0

Δχ² total (DSIC − ΛCDM) = −5.0 ;    BAO χ²/dof (dof = 8):  DSIC = 1.04 ,  ΛCDM = 1.36
```

DSIC describes the expansion geometry no worse, and formally slightly better, than ΛCDM; on the Jeffreys scale this is a weak, inconclusive preference, partly reflecting the known DESI–ΛCDM tension. The value `r_d ≈ 136–137 Mpc` is below the standard one here as a consequence of the high `H₀` from the SH0ES anchor (the product `r_d·H₀` is what is fixed); on the early scale the inverse problem gives `≈ 148 Mpc` (§7.2).

## 7.6 The `Om(z)` diagnostic: a null result

The diagnostic `Om(z) = (E² − 1)/((1+z)³ − 1)` [11] is a horizontal line for flat ΛCDM; DSIC predicts a U-shaped curve with a minimum near `z ≈ 2`. `Om(z)` is built from the data (chronometers + DESI `D_H/r_d`) in two covariance modes: diagonal and with the full Moresco et al. 2020 matrix (block `[4]`).

```
  H₀ |  Δχ² DIAG |  Δχ² FULL
-----------------------------
  64 |     −25.0 |     −24.2
  66 |     −15.6 |     −15.2
  68 |      −4.6 |      −4.7
  70 |      +8.3 |      +7.6
  72 |     +23.1 |     +21.6
  74 |     +39.9 |     +37.7
  76 |     +58.9 |     +55.8
```

**The sign of `Δχ²` flips with the choice of `H₀`** in both modes (the flip near `H₀ ≈ 68.7`): at a low normalization the data are formally closer to DSIC, at a high one — to ΛCDM. The result is set by the normalization, not by the shape of `Om(z)`; going from DIAG to FULL changes `Δχ²` by mere units. **Conclusion: the test does not discriminate between the models.** It becomes decisive only with an independently pinned `H₀`.

> The full matrix is a reconstruction of the Moresco et al. (2020) method on the official input tables, not a verbatim run of the authors' code; the diagonal coincides with the publication exactly.

## 7.7 Summary of the confrontation

| probe | data | result |
|---|---|---|
| supernovae `d_L(z)` | Pantheon+SH0ES, 1657 SNe | `Δχ² = +0.05` — indistinguishable |
| BAO `D_M/D_H` | DESI DR2, 6 points | `χ²/dof ≈ 0.78` |
| the `r_d` ruler | DESI DR2 | `≈ 148 Mpc`, static |
| chronometers `H(z)` | `z < 2` | `χ²/dof ≈ 0.37` |
| age `t₀` | integral over the history | `13.41 Gyr` (stellar `13.5 ± 0.5`) |
| transition `z_t` | `q(z)` reconstructions | `0.769` (observed `0.6–0.8`) |

A single initial condition `μ₀ ≈ 0.76` consistently describes six independent probes: one geometry, one parameter, no dark energy.

---

# 8. Distinguishing predictions and falsifiability

On late-time data DSIC and ΛCDM are practically degenerate; high `z` separates them. The parameters are taken from the joint fit (`μ₀ = 0.7610`, `Ω_m = 0.3029`), block `[5]`.

**The shape of `Om(z)`.** DSIC gives a U-shaped curve — a qualitative signature that flat ΛCDM lacks:

| `z` | 0.2 | 0.5 | 1.0 | 2.0 | 3.0 | 5.0 |
|---|---|---|---|---|---|---|
| `Om_DSIC` | 0.3566 | 0.3297 | 0.3064 | 0.3032 | 0.3266 | 0.4043 |
| `Om_ΛCDM` | 0.3029 | 0.3029 | 0.3029 | 0.3029 | 0.3029 | 0.3029 |

A flat `Om(z)`, measured precisely, would work against DSIC.

**Radial BAO in the Lyman-α region.**

| `z` | 2.33 | 2.5 | 3.0 | 3.5 |
|---|---|---|---|---|
| `D_H/r_d` DSIC | 8.612 | 7.978 | 6.449 | 5.300 |
| `D_H/r_d` ΛCDM | 8.599 | 8.013 | 6.615 | 5.572 |
| difference | +0.1 % | −0.4 % | −2.5 % | −4.9 % |

By `z ≈ 3.5` the divergence reaches ~5 %; the window `2.33 < z < 4.53` lies entirely below the Part II threshold, so this is a direct test of the **core**.

**The expansion rate beyond the data:** `E_DSIC/E_ΛCDM = 1.15` at `z = 5` and `1.47` at `z = 10`.

It is important to separate the sources of the signals: the `D_H` divergence in the window `2.33 < z < 4.53` is a trace of the core (`D_H = c/H` depends only on `μ₀`) and is independent of the Part II mechanics; the trace of the threshold `μ_d` lies in `D_M` above `z ≈ 4.5` (§10.5).

---
---

# Part II — The second tier: the binding of objects

> **Status.** Everything in this part rests on the core but is not derived from it: the poles `S, O` contain only the flow of scale, while the constructions below require a second, independent entity — the binding of objects. The core is self-contained: all six probes of §7 are independent of Part II, and the second-tier correction is identically zero throughout the range of their data. Part II introduces one new entity, one additional postulate and one measured constant.

---

# 9. The second entity and the detachment postulate

The core runs into one and the same missing entity from two sides. From the side of the early Universe: the geometry of the core falls short of the distance to the relic radiation by ~21 % — what is missing is the physics of the phase in which matter was still stuck together, and being stuck together is a state of binding. From the late side: the picture in which space accrues predominantly between weakly bound objects (Appendix B) requires a measure of that bond. Both superstructures are two manifestations of one entity; it is introduced once.

**P6 (the postulate of binding and detachment).**

- The early epoch is almost monolithic objectness, and the threshold refers to the **dominant (path-averaged) fraction**, not to every point. Until the monolith has on average parted (`μ > μ_d`), the internal observer's measuring chain is undefined: the ruler `L ∝ μ` is assembled from separate objects, while the medium is ruled by the stuck-together component. Light has almost nowhere to propagate — the path is counted not by the observer's metric but by the pure availability of space: an addition `∝ dμ/P = dμ/√(1−μ²)`. The weighting `b(μ₀)/b(μ)` of the late-time formula does not apply here — there is nothing to weight with.
- Local clumps in which the binding is strongest **do not uncouple** at the detachment of the medium: they retain their original stuck-together state and fall into themselves (§3.4). The first objects are not newly assembled but never-disintegrated patches of the monolith; their existence does not change the mean availability of space and does not shift the threshold.
- The detachment of the medium is a mechanical event: a **hard threshold** `μ_d`. At `μ ≤ μ_d` the correction is identically zero and the late-time geometry of the core works unchanged.
- The form is fixed by principle: the amplitude is natural (`A ≡ 1`, the form `1/P` with no prefactor), the window is sharp. Exactly one number remains free — the threshold `μ_d`.
- The threshold is a **detachment horizon** — the boundary at which the second pole first seeps through for a given fraction of the medium. It has two consequences: a *spatial* one — the optical correction to the light path (§10), and a *temporal* one — the domain of the clock factor `dτ = dt/μ` (§5.1, §11). Detachment is local: each detached patch switches to the relative clock at the moment of its own detachment, not at the common threshold.

**Consistency with observations of early galaxies.** JWST galaxies are observed out to `z ≈ 10–14` (`μ_e(z=14) ≈ 0.9972`), i.e. above the threshold `μ_d = 0.9806`. A literal reading — "no separate objects exist at `μ > μ_d`" — would contradict these observations; in the adopted formulation `μ_d` is the detachment horizon of the *dominant* fraction, while the first objects are patches of the monolith that never uncoupled: where the binding is maximal, the original stuck-together state persists as the surrounding medium detaches, and such a patch separates from it as an island earlier than the mean threshold. Their existence above the mean threshold is not a contradiction but the expected early islands: detachment is a process distributed over objects and epochs, not one global event, and objects in it are not assembled anew — they merely un-stick at different times.

**A physical analogy.** A set of objects is compressed in a press; when the press releases, the objects first push against one another, and at some moment they come unstuck — after which there is no more pressure on the walls. That point of unsticking is `μ_d`; individual objects deep inside may never come unstuck at all — remaining stuck together when everything else has already detached — and this does not affect the threshold. The sharpness of the threshold is motivated mechanically: a smooth window would introduce an extra parameter without grounds.

**Naturalness of the form.** The base integrand of the core path equals `(√(1−μ₀²)/μ₀²)·μ/P`, and the integrand of the correction is `1/P`; at the Mush (`μ → 1`) they coincide in form. The addition in the stuck-together phase is asymptotically comparable to the base measure of the path.

**Status of `μ_d`.** The value of the threshold is not derived from the core: there is no number in the poles `S, O` singling out `0.9806`. But with the form fixed by P6, the "amplitude ↔ threshold" degeneracy is broken by principle, and `μ_d` becomes a **measured constant of the second tier** — of the same status as `μ₀` in the core. A theory is allowed to have measured constants; it is not allowed arbitrary functions.

---

# 10. The early Universe: the detachment threshold

## 10.1 The core's shortfall and the calibration of the threshold

The Planck 2018 reference data [7] (column `base_plikHM_TTTEEE_lowl_lowE_lensing`, parameter tables from the Planck Legacy Archive) pass a self-check (block `[7]`): `z* = 1089.92 ± 0.25`, `r* = 144.43 ± 0.26 Mpc`, `100θ* = 1.04110 ± 0.00031`, whence `d_LSS = r*/θ* = 13872.8 ± 25 Mpc` (a 0.00 % discrepancy with the hard-coded value).

The core parameters on the early scale (the same as for BAO, the chronometers and the age): `μ₀ = 0.76`, `H₀ = 67.4`, `c/k = 14978.2 Mpc`. Then (block `[8]`):

```
μ_cmb(z* = 1089.92) = 0.9999994681
D_M base to the CMB = 10936.2 Mpc
shortfall           = 2936.6 Mpc  (21.2 %)
```

The cause of the shortfall lies on the core's side: the particle horizon is finite "for free" (§5.5), almost the whole path is accumulated already by `z ≈ 10`, and beyond that the distance integral runs out of steam.

By P6 the addition to the path is `1/P`, cut off by the threshold; its integral is `arcsin μ`:

```
D_M_corr(μ_e) = D_M(μ_e) + (c/k)·[arcsin μ_e − arcsin μ_d]   for μ_e > μ_d ,
D_M_corr(μ_e) = D_M(μ_e)                                      for μ_e ≤ μ_d .
```

Closure to `θ_model(μ_d) = θ*` gives

```
μ_d = 0.980640 ,   z_detach = 4.526 ,   100·θ_model = 1.04110   (Planck: 1.04110) .
```

## 10.2 What the postulate fixes: the degeneracy without P6

Without fixing the amplitude, the equation `A·(c/k)·[arcsin μ_cmb − arcsin μ_d] = shortfall` has a curve of solutions — any amplitude yields its own `μ_d` with the same distance to the CMB (block `[9]`):

| `A` | 0.5 | 0.8 | **1.0** | 1.2 | 2.0 | 5.0 |
|---|---|---|---|---|---|---|
| `μ_d` | 0.923705 | 0.969868 | **0.980640** | 0.986514 | 0.995097 | 0.999190 |
| `z_detach` | 1.506 | 3.344 | **4.526** | 5.690 | 10.266 | 26.922 |

One observation fixes one parameter: **zero degrees of freedom**. The value `z_detach ≈ 4.53` is a calibration under a form fixed by principle, not a prediction, and it is not presented as one.

## 10.3 Properties of the optical correction

The correction is **optical**, not dynamical: the light path is lengthened, the rate `H(z)` is unchanged. The consequences hold by construction:

- **Frequency is untouched.** The redshift remains `1+z = b(μ₀)/b(μ_e)`.
- **Distance duality is preserved.** `d_L` and `d_A` are built from the same `χ_corr`, so the Etherington relation holds automatically.
- **Achromaticity and non-dissipativity.** The blackbody character of the relic spectrum (FIRAS) excludes any medium that is frequency-dependent or scattering; the correction must be a geometry of the path — the form `1/P` satisfies this.
- **The particle horizon remains finite:** `χ_corr(μ_e → 1) = 13905.7 Mpc` against the base `10953.6` — only 32.9 Mpc beyond the last-scattering surface: the relic radiation lies almost flush against the horizon.

## 10.4 The late Universe is untouched

The correction is strictly zero below the threshold — this is an identity, not an approximation (block `[9]`):

| `z` | 0.50 | 1.00 | 2.33 | 4.00 | 4.53 | 4.63 | 5.00 | 10.00 |
|---|---|---|---|---|---|---|---|---|
| `D_M` base, Mpc | 1932.8 | 3384.2 | 5792.7 | 7335.6 | 7653.3 | 7707.8 | 7896.9 | 9247.2 |
| correction, Mpc | +0.0 | +0.0 | +0.0 | +0.0 | +0.0 | +49.3 | +220.4 | +1433.0 |

All six DESI DR2 bins (`z ≤ 2.33`) are insensitive to the correction, the BAO `χ²/dof` is unchanged — the results of §7 are fully preserved.

## 10.5 The trace of the threshold and the reach of the probes

The relative trace of the threshold in `D_M` (block `[10]`):

| `z` | 4.53 | 5.00 | 6.00 | 8.00 | 11.00 | 1089.92 |
|---|---|---|---|---|---|---|
| object | threshold | JWST galaxies | quasars | reionization | first galaxies | CMB (calibration) |
| correction to `D_M` | +0.0 % | +2.8 % | +7.1 % | +12.4 % | +16.6 % | +26.9 % |

The reach of the probes: DESI DR2 (`z ≤ 2.33`), DESI DR3 (~2027, `z ≤ 3.55`), DESI-II/LBG (~2029) — all below the threshold. The `D_H/r_d` divergence in the window `2.33 < z < 4.53` is a trace of the core, not of `μ_d`. The threshold itself affects only `D_M` above `z ≈ 4.5`, where no precise standard rulers exist apart from the CMB calibration: a second independent anchor for `μ_d` is out of reach in accessible data until a precise distance scale appears at `z ≳ 5–6`.

## 10.6 Reading of the relic imprint

The imprint is a freeze-frame of the shape of objectness at a late stage of the Mush, when space had already begun to seep through and light could escape the monolith for the first time: the contrast on the map is a contrast in the density of the seeping space. The Mush at that point is not a strictly zero state but a step back from zero (`μ_cmb ≈ 0.9999995`). The sign of the contrast: cold patches ↔ dense objectness, hot ones ↔ pores of space that have seeped through; the smallness of the contrast, `~10⁻⁵`, agrees qualitatively with the pores being only nascent at the threshold. This sign coincides with the standard large-scale Sachs–Wolfe effect — a pleasant consistency, but not evidence, since neither the amplitude nor the spectrum is derived from the model.

An important distinction: `μ_cmb` is the observable "floor" of the early zone, `μ_d` its late edge.

## 10.7 What the second tier settles and what it does not

Exactly one number is closed — the distance to the last-scattering surface, precisely and with no consequences for the late Universe. The following remains open.

**The borrowed sound horizon.** The closure uses `r* = 144.43 Mpc` — a quantity not measured directly but computed within ΛCDM from the physics of the photon–baryon plasma, which DSIC does not have. The directly measured quantity is the angle `θ*`, and the calibration is performed against it; the origin of the acoustic scale `~147 Mpc` itself — consistently showing through in the late-time ruler `r_d ≈ 148 Mpc` as well — is declared an open question of the core.

**The relic radiation is not one number.** Planck constrains the whole structure of the acoustic peaks, the damping tail and the polarization — thousands of multipoles. The second tier closes one parameter; everything else requires a theory of perturbations and of composition, which the model does not have.

---

# 11. Cosmic time and the detachment horizon

## 11.1 Statement

The optical correction of §10 changes the path but not `H(z)`: at `z = 5–10` the rate remains above ΛCDM by +15…+47 %. If proper time is counted with the relative clock along the whole path from the Mush, the time budget of the early epochs turns out compressed (0.42 Gyr by `z ≈ 7.5` against 0.71 in ΛCDM). Such an application of the clock factor is, however, inconsistent: `dτ = dt/μ` is defined only below the detachment horizon (§5.1, §9), while in the connected monolith above the threshold the clock rate is normal. Below is the computation with the domain of the factor taken into account.

## 11.2 Piecewise background time

Below the horizon (`μ < μ_d`) the clock is relative; above it, normal. The threshold is taken from the same closure to `θ*` (§10.1) and is not a new parameter. Matching the rate at the threshold makes `t(z)` continuous (`Δ ≈ 6.5·10⁻⁴` Gyr); the late Universe (`z ≤ z_detach`) is identically untouched (block `[6]`):

| `z` | `t` relative | `t` piecewise | `t_ΛCDM` | zone |
|---|---|---|---|---|
| 3.0 | 1.73 | 2.70 | 2.18 | below the threshold (untouched) |
| 4.0 | 1.15 | 2.12 | 1.57 | below the threshold |
| 5.0 | 0.81 | 1.78 | 1.19 | above (normal clock) |
| 6.0 | 0.60 | 1.53 | 0.95 | above |
| 7.5 | 0.42 | 1.27 | 0.71 | above |
| 10.0 | 0.25 | 0.99 | 0.48 | above |
| 12.0 | 0.18 | 0.84 | 0.37 | above |

## 11.3 The island time of the observer

The piecewise curve above is the clock of the **background**. But the relativity of time is a property of each detached patch, switched on at the moment of its own detachment `z_form`: before it the matter is in the connected monolith (normal clock), after it a free island (relative clock). The proper time of a given object therefore depends on the `z_form` of its host structure.

Here `z_form` is an **input from observations** (the assembly epoch of a galaxy), not a parameter of the model: the core sets the rule for the clock rate after detachment, but not the epoch of detachment. The latter is a matter for structure-formation theory, exactly as the Friedmann background does not derive the assembly epoch of a given galaxy.

The construction has a built-in self-regulation: the earlier an island detached, the longer it lived on the compressed relative clock and the less proper time it accumulated; the later it detached, the longer its monolithic phase with the normal clock. The two effects pull in opposite directions and stabilize the budget at about 0.5–0.8 Gyr for almost any island, with no tuning.

The age today for an early-detached observer (`z_form = 10–20`) is `t₀ = 13.5–13.6 Gyr` — falling inside the model-independent stellar anchor `13.5 ± 0.5 Gyr`, without a circular tie to the Planck value of `13.80`.

## 11.4 A falsifiable test on observed SMBHs

For each early supermassive black hole, `M_BH` and `z_obs` are known. Salpeter growth (`0.045` Gyr per e-fold, Eddington-limited, a light seed of `10² M☉`) sets the required proper time `t_need`; the island picture gives the maximum available `t_avail` (detachment immediately before the observation) and the limiting detachment epoch. Since earlier detachment means less proper time, the constraint reads: the host structure must have detached **late enough** — at `z_form ≲ X`, i.e. no earlier than the epoch of redshift `X`. If even the maximum `t_avail < t_need`, the object cannot be accommodated with a light seed. The redshifts of the objects are spectroscopic and the masses are the virial/accretion estimates of the discovery publications [16–23] (block `[6b]`):

| object | `z_obs` | `M_BH`, M☉ | e-folds | `t_need`, Gyr | max `t_avail` | verdict |
|---|---|---|---|---|---|---|
| UHZ1 [16–18] | 10.073 | 4·10⁷ | 12.9 | 0.58 | 0.49 | a heavy seed is required |
| GN-z11 [19] | 10.603 | 1.6·10⁶ | 9.7 | 0.44 | 0.45 | OK if `z_form ≲ 11.0` |
| CEERS 1019 [20] | 8.679 | 1·10⁷ | 11.5 | 0.52 | 0.64 | OK if `z_form ≲ 11.4` |
| J0313−1806 [21] | 7.642 | 1.6·10⁹ | 16.6 | 0.75 | 0.80 | OK if `z_form ≲ 8.4` |
| J1342+0928 [22] | 7.541 | 7.8·10⁸ | 15.9 | 0.71 | 0.82 | OK if `z_form ≲ 9.0` |
| J1120+0641 [23] | 7.085 | 2·10⁹ | 16.8 | 0.76 | 0.91 | OK if `z_form ≲ 9.1` |

Five of the six objects are accommodated with detachment at `z_form ≲ 8–11` — observationally reasonable assembly epochs. UHZ1 does not pass with a light seed; however, it is precisely this object that serves in standard cosmology as the main argument for direct-collapse seeds of `10⁴–10⁵ M☉`, for which the required number of e-folds drops and the object passes — the interpretation of UHZ1 as evidence for a heavy seed is given by the discoverers themselves [16, 17]. The model arrives at the same conclusion as standard astrophysics.

> **A caveat on the data.** The values are verified against the discovery publications. Quasar masses are virial estimates that depend on the method (Mg II or C IV) and the calibration and differ by factors of a few between works, which can shift the limiting `z_form`. For J1342+0928 the discovery value is `~8·10⁸ M☉` [22], while `7.8·10⁸` is the refined virial estimate of later works; for J1120+0641 the conservative discovery value [23] is used, refined by JWST observations. This precision suffices to test the logic of the model; for quantitative claims one should adopt the value of the specific publication cited.

**The statuses are separated explicitly.** (1) *Derived from the core and falsifiable:* the rule for the clock rate — the relativity of time switches on at the detachment horizon. (2) *An input, not a derivation:* the detachment epoch of a specific object. A systematic failure — objects requiring detachment earlier than the observed assembly epochs for any reasonable seed — would falsify the island picture of time.

---

# 12. Structure growth: a check on the DSIC background

Structure growth is not part of the core. A minimal step is taken here: it is checked whether the DSIC **background** is compatible with the standard growth of linear perturbations — not by deriving the growth equation from the core, but by taking it in the GR limit on the DSIC geometry (block `[11]`).

On a background with rate `H(z)`, the matter contrast in the sub-horizon GR limit obeys the standard equation

```
d²δ/d(ln a)² + [2 + dln H/dln a] · dδ/d(ln a) − (3/2) Ω_m(a) · δ = 0 ,
```

where `Ω_m(a) = Ω_m0·a⁻³/E(a)²`. Since the core introduces no densities, `Ω_m0` enters as a growth parameter. The comparison is with the "Gold-2017" compilation (18 points [13, 14]; for the three WiggleZ points [15] their covariance is included).

The `fσ₈` data are degenerate: they constrain not `Ω_m0` and `σ₈` separately but the combination `S₈ = σ₈·√(Ω_m0/0.3)`. A scan along the degeneracy valley:

| `Ω_m0` | 0.15 | 0.20 | 0.25 | 0.30 | 0.35 | 0.40 |
|---|---|---|---|---|---|---|
| `σ₈*` | 1.040 | 0.925 | 0.835 | 0.770 | 0.715 | 0.670 |
| `χ²/dof` | 0.72 | 0.70 | 0.73 | 0.78 | 0.84 | 0.91 |
| `S₈` | 0.735 | 0.755 | 0.762 | 0.770 | 0.772 | 0.774 |

```
Invariant along the valley:  S₈ = 0.761 ± 0.013
Slice at Ω_m0 = 0.30:        σ₈ = 0.770 ,  χ²/dof = 0.73
ΛCDM on the same points:     S₈ = 0.760 ,  χ²/dof ≈ 0.74   (Δχ² ≈ +0.1)
```

The DSIC background describes structure growth no worse than ΛCDM; the recovered `S₈ ≈ 0.76` agrees with independent estimates (Planck `≈ 0.83`, weak lensing `≈ 0.76`). This does not prove that growth in DSIC must follow the GR-limit equation — the core does not dictate it — but it removes the objection of gross incompatibility.

> **Simplifications.** (1) The Alcock–Paczynski correction for the difference between the DSIC background distances and the fiducial cosmology of each RSD measurement was not applied; at `z < 1` the factor is usually a few per cent and small against the errors. (2) The growth equation is the GR-limit one; deriving growth from binding (Appendix B) remains a direction of development.

---

# 13. Limits of applicability and the relation to ΛCDM

**What has been done.** The background kinematics of the late Universe (`z < 2.3`): one parameter `μ₀` describes six independent probes statistically no worse than flat ΛCDM at an equal number of free parameters. The second tier closes the distance to the relic radiation with one measured constant, without touching the late Universe.

**What has not been done.**

- *Structure growth* is not derived from the core; only the compatibility of the background with standard GR-limit growth has been checked (§12).
- *The origin of the acoustic scale* `~147 Mpc` is an open question; `r*` is borrowed from ΛCDM (§10.7).
- *The detachment epoch* `z_form` of a specific object is an input from observations, not a derivation from the core (§11.3); a physical criterion for detachment is a subject for separate study.
- *Gravity* (Appendix B) is a phenomenological sketch with a free parameter; no comparison with data has been made.
- *Ontology versus ansatz.* The form of `b(μ)` is derived here from the circular geometry of the poles, but it can equally be adopted simply as a one-parameter ansatz: the formulas work on their own, and the ontological picture is not required for them to work.

**Relation to ΛCDM.** DSIC does not refute ΛCDM. On late-time data the models are statistically indistinguishable; DSIC offers a different ontology of the same background — geometry instead of dark energy — at an equal number of free parameters. Its value lies not in the quality of the fit (which is no better) but in the origin of the curve: the same distance–redshift relation follows from a single geometric starting point, without `Ω_Λ`, and moreover the deceleration → acceleration transition and the finiteness of the particle horizon are derived rather than postulated.

**Where the model will meet decisive tests.**

- the shape of `Om(z)` with an independently pinned `H₀`;
- radial BAO in the Lyman-α region at `z ≈ 3–3.5` (up to ~5 %, within reach of DESI DR3);
- the test of §11.4 on a growing sample of early SMBHs as masses and host-galaxy assembly epochs are refined.

---

# 14. Conclusion

From four structural elements — the relativity of scale, the conservation of the norm of two poles, a single uniform flow, and the modulus-only access of the internal observer — the scale factor `b(μ) = √(1−μ²)/μ²` follows uniquely, and from it the entire background kinematics: the redshift `1+z = b(μ₀)/b(μ_e)`, the rate `|H_τ| = (2−μ²)/(1−μ²)`, a sign-changing deceleration parameter with the exact crossover `μ_cross = √(3−√5)`, closed-form distances and a finite particle horizon. The single parameter `μ₀ = 0.7600 ± 0.0091` consistently passes six independent probes of the late Universe, remaining statistically indistinguishable from ΛCDM on its best data and sharply distinguishable at high `z`. Acceleration of the expansion in this picture is not the action of a hidden component but the shadow of the divergence of two scales, read off by an observer who is itself made of one of them.

The second tier adds one entity — the binding of objects — and one measured constant `μ_d = 0.9806`. The late Universe is thereby left identically untouched, the distance to the last-scattering surface is closed exactly, and distance duality and achromaticity hold by construction. The same threshold, read as a detachment horizon, also sets the domain of the clock factor: clocks are normal in the connected monolith and relative after detachment. This turns the cosmic-time budget of the early epochs from a global deficit into an individually testable requirement on the assembly epoch of each early structure.

Current data neither reject the model nor support it more strongly than ΛCDM. The decisive evidence will come from high `z`, an independently pinned `H₀`, and a growing sample of early supermassive black holes.

---

# Reproducibility

All empirical claims are reproduced by a single self-contained script, `dsic2_test.py`, from the repository https://github.com/az-zase/DSIC. The script downloads the public catalogs itself and computes `χ²` with the full covariance; the provenance of every hard-coded figure is documented in place. A full run takes 10–15 minutes.

**Data.** Pantheon+SH0ES (`Pantheon+SH0ES.dat`, `STAT+SYS.cov`, `STATONLY.cov`) from the `PantheonPlusSH0ES/DataRelease` repository (auto-download); DESI DR2 BAO (Table IV, 6 bins, hard-coded with a self-check); cosmic chronometers — a compilation of 33 points with the full Moresco et al. (2020) covariance for 15 of them; Planck 2018 (`z*`, `r*`, `θ*`); Gold-2017 `fσ₈` (18 points).

**Script blocks.**

| block | what it does | section |
|---|---|---|
| `[0]` | self-check of the hard-coded DESI against Table IV | §7.2 |
| `[0b]` | machine verification of the core's internal identities | §6 |
| `[0c]` | the spinor form `ψ = R·e^{iαφ}`: 8 identities at α=1/2 + the forcedness test of the covering (cycle ratio = 2 for any α) | §3.3, §6 |
| `[1]` | the anchored SN fit (STAT+SYS), explicit `M`, 1σ from the Hessian | §7.1 |
| `[2]` | the STAT-ONLY control | §7.1 |
| `[3]` | joint SN+BAO, free `r_d`, the `r_d` ruler | §7.2, §7.5 |
| `[4]` | `Om(z)`: robustness to `H₀`, DIAG and FULL covariances | §7.6 |
| `[5]` | falsifiable predictions at high `z` | §8 |
| `[6]` | cosmic time `t(z)`, the piecewise clock | §11.2 |
| `[6b]` | the island test on observed SMBHs | §11.4 |
| `[7]` | self-check of the Planck reference data | §10.1 |
| `[8]` | the core's shortfall, fixing `μ_d` via `θ*` | §10.1 |
| `[9]` | the zero-below-threshold control; the `A ↔ μ_d` degeneracy | §10.2, §10.4 |
| `[10]` | the trace of `μ_d`, the horizon, the reach of the probes | §10.3, §10.5 |
| `[11]` | structure growth `fσ₈` on the DSIC background (GR limit), `S₈` | §12 |

Dependencies: `numpy`, `scipy`.

---

# Appendix A. Summary of notation

| layer | symbol | meaning | sign / range |
|---|---|---|---|
| 0 | `φ = ωt` | the phase of the flow (time) | grows monotonically |
| 0 | `ω` | the rate of the flow | `> 0` |
| 1 | `R` | the norm (the scale budget) | `> 0` |
| 1 | `S = R·sin(φ/2)` | the scale of space | signed (odd) |
| 1 | `O = R·cos(φ/2)` | the scale of objectness | signed (even) |
| 1 | `σ_S, σ_O` | the sign phase of the poles (= the choice of covering sheet) | `±1` |
| 1 | `ψ = R·e^{iφ/2}` | the spinor-complex core, `O = Re ψ`, `S = Im ψ` (§3.3) | `\|ψ\| = R`, period `4π` |
| 1 | `O_full = O·σ_S` | objectness with the inherited sign (`Re ψ` on the sheet set by the sign of `Im ψ`) | signed |
| 1 | `W = −R²/2` | the Wronskian (quadrature) | const |
| 2 | `P = \|S\| = R√(1−μ²)` | the observable scale of space | `≥ 0` |
| 2 | `Q = \|O\| = R·μ` | the observable scale of objectness | `≥ 0` |
| 2 | `μ = Q/R = \|cos(φ/2)\|` | the observable phase (the fraction of objectness) | `∈ [0,1]`, period `2π` |
| 3 | `L(μ) ∝ μ` | the observer's ruler | — |
| 3 | `τ`, `dτ = dt/μ` | proper time (below the detachment horizon) | — |
| 3 | `b = √(1−μ²)/μ²` | the observable scale factor | — |
| 3 | `z` | redshift, `1+z = b(μ₀)/b(μ_e)` | `≥ 0` |
| 3 | `s = ±1` | the branch direction | — |
| 3 | `k = (2/π)ω` | the rate of `μ` along a branch | `> 0` |
| 3 | `H_τ`, `q` | the expansion rate, the deceleration parameter | — |
| 3 | `w_eff = (2q−1)/3` | the effective equation of state (§5.6) | — |
| 3 | `χ`, `d_L`, `DM` | distances and the distance modulus | — |
| 3 | `μ₀, μ_e` | the phase today / at emission | physical branch `μ_e ∈ (μ₀, 1)` |
| 3 | `μ_cross = √(3−√5)` | the zero of `q` | `0.874032`, `z_t = 0.769` |

**Notation of Part II:**

| symbol | meaning | status |
|---|---|---|
| binding | the second entity of the second tier | postulated (P6) |
| `μ_d = 0.980640` | the detachment horizon (`z_detach = 4.526`) | a measured constant (calibrated against `θ*`) |
| `A ≡ 1` | the amplitude of the `1/P` correction | fixed by postulate P6 |
| `z_form` | the detachment epoch of a specific structure | an input from observations (§11.3) |
| `D_M_corr` | the distance with the `arcsin` correction added | — |
| `r_ij`, `M_i`, `G`, `κ` | quantities of the module of Appendix B | not from the core; `κ` is free |

---

# Appendix B. A development programme: gravity and the distribution of space

> **Status.** This is not a result but a sketch of a direction. The module is phenomenological: it uses the Newtonian quantities `G, M_i`, which are absent from the flow `S ↔ O`; the parameter `κ` is free; no comparison with data has been carried out. The value of the module is architectural: it shows that binding may also have a late, gravitational face, closing the early (the detachment threshold) and the late (the distribution of space) onto one entity.

When objects fall into themselves (§3.4), gravitationally bound groups fall in concert — a common center of fall shows through for them: the more strongly the objects are bound, the less space is released between them. The loop is self-reinforcing: dispersal weakens the bond → the region receives more space → the dispersal accelerates. Dense structures stay compact, rarefied ones receive more of the space that seeps through.

```
binding:                 g_ij = G·M_i·M_j / r_ij²
potential gradient:      |∇Φ_ij| = | Σ_k G·M_k·(mid_ij − pos_k)/|mid_ij − pos_k|³ |
weight of space:         w_ij = (1/|∇Φ_ij|) / Z ,     Z = Σ_{i≠j} 1/|∇Φ_ij|
increment:               Δr_ij = κ·|ΔP|·w_ij / (M_i + M_j)
conservation:            Σ_{i≠j} Δr_ij·(M_i + M_j)/κ = |ΔP(φ)|
```

Properties of the weights: `w_ij ≥ 0`, `Σ w_ij = 1`; a pair in a void receives a large weight, a pair in a cluster a small one.

The minimum programme: derive `μ_d` as a condition on binding from the same root as the weights `w_ij` — then both constructions become consequences of one quantity rather than two add-ons.

---

# References

1. A. G. Riess et al., *Observational Evidence from Supernovae for an Accelerating Universe and a Cosmological Constant*, AJ **116**, 1009 (1998).
2. S. Perlmutter et al., *Measurements of Ω and Λ from 42 High-Redshift Supernovae*, ApJ **517**, 565 (1999).
3. D. M. Scolnic et al., *The Pantheon+ Analysis: The Full Dataset and Light-Curve Release*, ApJ **938**, 113 (2022); arXiv:2112.03863.
4. D. Brout et al., *The Pantheon+ Analysis: Cosmological Constraints*, ApJ **938**, 110 (2022); arXiv:2202.04077.
5. A. G. Riess et al., *A Comprehensive Measurement of the Local Value of the Hubble Constant* (SH0ES), ApJL **934**, L7 (2022); arXiv:2112.04510.
6. DESI Collaboration (M. Abdul-Karim et al.), *DESI DR2 Results II: Measurements of Baryon Acoustic Oscillations and Cosmological Constraints*, Phys. Rev. D **112**, 083515 (2025); arXiv:2503.14738.
7. Planck Collaboration (N. Aghanim et al.), *Planck 2018 results. VI. Cosmological parameters*, A&A **641**, A6 (2020); erratum A&A **652**, C4 (2021); arXiv:1807.06209. (The column `base_plikHM_TTTEEE_lowl_lowE_lensing` is used; parameter tables: Planck Legacy Archive, `Baseline_params_table_2018_68pc`.)
8. E. Di Valentino et al., *In the Realm of the Hubble Tension — a Review of Solutions*, Class. Quantum Grav. **38**, 153001 (2021); arXiv:2103.01183.
9. M. Moresco et al., *Setting the Stage for Cosmic Chronometers*, ApJ **898**, 82 (2020); arXiv:2003.07362. (With the compilation of Moresco et al. 2012; Moresco 2015; Moresco et al. 2016; primary sources of the other `H(z)` points: Simon et al. 2005; Stern et al. 2010; Zhang et al. 2014; Ratsimbazafy et al. 2017; Borghi et al. 2022.)
10. O. Farooq & B. Ratra, *Hubble Parameter Measurement Constraints on the Cosmological Deceleration–Acceleration Transition Redshift*, ApJ **766**, L7 (2013); arXiv:1301.5243.
11. V. Sahni, A. Shafieloo & A. A. Starobinsky, *Two New Diagnostics of Dark Energy*, Phys. Rev. D **78**, 103502 (2008); arXiv:0807.3548.
12. V. Sahni, T. D. Saini, A. A. Starobinsky & U. Alam, *Statefinder — a New Geometrical Diagnostic of Dark Energy*, JETP Lett. **77**, 201 (2003); arXiv:astro-ph/0201498.
13. S. Nesseris, G. Pantazis & L. Perivolaropoulos, *Tension and Constraints on Modified Gravity Parametrizations of G_eff(z) from Growth Rate and Planck Data*, Phys. Rev. D **96**, 023542 (2017); arXiv:1703.10538.
14. B. Sagredo, S. Nesseris & D. Sapone, *Internal Robustness of Growth Rate Data*, Phys. Rev. D **98**, 083543 (2018); arXiv:1806.10822.
15. C. Blake et al., *The WiggleZ Dark Energy Survey: Joint Measurements of the Expansion and Growth History at z < 1*, MNRAS **425**, 405 (2012); arXiv:1204.3674.
16. Á. Bogdán et al., *Evidence for Heavy-Seed Origin of Early Supermassive Black Holes from a z ≈ 10 X-ray Quasar* (UHZ1), Nat. Astron. **8**, 126 (2024); arXiv:2305.15458.
17. P. Natarajan et al., *First Detection of an Over-massive Black Hole Galaxy Candidate: UHZ1*, ApJL **960**, L1 (2024); arXiv:2308.02654.
18. A. D. Goulding et al., *UNCOVER: The Growth of the First Massive Black Holes from JWST/NIRSpec — Spectroscopic Redshift of UHZ1, z = 10.073*, ApJL **955**, L24 (2023); arXiv:2308.02750.
19. R. Maiolino et al., *A Small and Vigorous Black Hole in the Early Universe* (GN-z11), Nature **627**, 59 (2024); arXiv:2305.12492.
20. R. L. Larson et al., *A CEERS Discovery of an Accreting Supermassive Black Hole 570 Myr after the Big Bang* (CEERS 1019), ApJL **953**, L29 (2023); arXiv:2303.08918.
21. F. Wang et al., *A Luminous Quasar at Redshift 7.642* (J0313−1806), ApJL **907**, L1 (2021); arXiv:2101.03179.
22. E. Bañados et al., *An 800-Million-Solar-Mass Black Hole in a Significantly Neutral Universe at a Redshift of 7.5* (J1342+0928), Nature **553**, 473 (2018); arXiv:1712.01860.
23. D. J. Mortlock et al., *A Luminous Quasar at a Redshift of z = 7.085* (J1120+0641), Nature **474**, 616 (2011).
