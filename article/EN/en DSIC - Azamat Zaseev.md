# Dual Scale Inversion Cosmology (DSIC): a one-parameter background cosmology from the relativity of scale

**Azamat Zaseev**

Independent Researcher, Russia

ORCID: 0009-0004-9345-1069

https://doi.org/10.5281/zenodo.21198494

https://github.com/az-zase/DSIC

---

## Abstract

I present a geometric model of the cosmological background — Dual Scale Inversion Cosmology (DSIC) — in which space and objectness are treated as two poles of a single conserved quantity, and no absolute scale exists: only the ratio of the poles carries physical content. The state of the system moves along the circle `S² + O² = R²` under a single uniform flow; an internal observer, itself built of objectness, has access only to the projection. From the ruler and clock of such an observer, the observable scale factor `b(μ) = √(1−μ²)/μ²` is derived uniquely, and from it — redshift, kinematics, and distances. The late-time acceleration of the Universe turns out to be a kinematic consequence of the relative scale drift of the two poles; dark energy and analogues of `Ω_m, Ω_Λ` are not introduced into the model.

A single parameter `μ₀ = 0.7600 ± 0.0091` (the observer's phase "today"; 1σ from the χ² Hessian) consistently describes six independent probes of the background kinematics at `z < 2.3`. The form of `b(μ)` is fixed by construction and contains no functional freedom, so the value of `μ₀` measured from late-time data uniquely determines the model's behavior at high `z`: the departures from ΛCDM at `z > 3` are forced predictions, not the result of tuning. On the full Pantheon+SH0ES catalog [3, 4] with the STAT+SYS covariance, the model is statistically indistinguishable from flat ΛCDM (`Δχ² ≈ +0.05` with an equal number of free parameters); the joint SN+BAO fit formally gives a weak preference to DSIC (`Δχ² ≈ −5.0`, partly reflecting the known DESI–ΛCDM tension). The `Om(z)` diagnostic [11] does not discriminate between the models on current data (a null result: the sign of `Δχ²` depends on the `H₀` normalization). The most stringent tests of the core: the radial BAO scale in the Lyman-α region (a departure from ΛCDM of up to ~5% at `z ≈ 3.5`, within reach of DESI DR3) and the cosmic time budget at `z = 5–12` (at `z = 7.5`, 0.42 Gyr are available versus 0.71 in ΛCDM), analyzed quantitatively in §11.6.

The early Universe is delegated to the **second tier** of the model (Part II): the form of the correction to the light path is derived from an additional postulate on the binding of objects, and the single threshold `μ_d ≈ 0.9806` (`z ≈ 4.5`) is calibrated against a single observation — the distance to the last-scattering surface; the status of this superstructure — a calibration with the form fixed by the postulate, not a prediction — is stated explicitly. Structure growth is not derived from the core; the model's background has been checked for compatibility with the standard (GR-limit) growth of linear perturbations: on the `fσ₈` compilation ("Gold-2017" [13, 14]) one recovers `S₈ ≈ 0.76`, statistically indistinguishable from ΛCDM (§12); deriving growth from binding remains a direction for development. All empirical claims are reproducible with scripts from the open repository.

---

## 1. Introduction

The transition of the cosmological expansion from deceleration to acceleration at `z_t ≈ 0.6–0.8` [1, 2] is a firmly established observational fact. The standard model describes it by introducing a component with negative pressure: two density parameters `Ω_m` and `Ω_Λ` set the shape of the expansion curve, but the constant Λ itself remains a postulate without derivation. A natural question arises: can the same expansion curve emerge without an energy ingredient — from nothing but the geometry of what I call "scale"?

DSIC answers this question in the affirmative, in the form of an explicit construction. The model starts from the premise that no absolute scale exists, either for objects or for space: only their ratio is observable. Space and objectness are treated as two poles of a single conserved norm, flowing along a circle under a single uniform flow. An observer located inside the system and itself made of objectness has, in principle, no access to the sign of quantities or to their absolute values — only to moduli and ratios. The requirement that the speed of light be constant in the observer's units uniquely fixes their clock, and from the ruler and the clock the observable scale factor is derived:

```
b(μ) = √(1 − μ²) / μ² ,
```

where `μ ∈ [0, 1]` is the observable fraction of objectness. This single function generates the redshift, the expansion rate, the deceleration parameter, and all distances. The sign change of the deceleration parameter — the transition from deceleration to acceleration — is built into the curvature of `b(μ)` and requires neither Λ nor any new component: acceleration arises as a kinematic consequence of the growing scale drift of the two poles.

The claimed domain of applicability is stated at once. What the core establishes is the **background kinematics of the late Universe** (`z < 2.3`): one parameter `μ₀` describes six independent probes statistically no worse than flat ΛCDM with an equal number of free parameters. The key property of the model is not economy but rigidity: the form of `b(μ)` is fixed by construction and contains no functional freedom, so the value of `μ₀`, measured from late-time data to ~1% precision, uniquely determines the model's behavior at high `z` — including `E(z)`, the radial BAO scale in the Lyman-α region, and the cosmic time budget at `z = 5–12` — with no possibility of tuning. These predictions depart from ΛCDM at a level within reach of DESI DR3 and constitute a direct test of the core (§11.6). The early Universe is reached by the model's **second tier** (Part II): one new entity — the binding of objects — fixes the form of the correction to the light path by postulate, and the single threshold `μ_d` is calibrated against a single observation — the distance to the CMB; the status of this superstructure (a calibration, not a prediction) and its limitations are examined in §11.6. Structure growth is not part of the core; the model's background has been checked only for compatibility with the standard growth of linear perturbations (§12).

The exposition follows a strict layer discipline, and the layer a symbol belongs to is always visible from its form:

- **layer 0 — the flow**: the single monotonic parameter `φ` (time);
- **layer 1 — the core**: signed quantities `S, O, σ`; there are no forces, cycles, or reversals here, and the modulus operation is not applied;
- **layer 2 — the projection**: unsigned quantities `P, Q, μ`, accessible to the internal observer; they are born of the single operation `|·|`; bounces and cycles exist only here;
- **layer 3 — observable physics**: `b, z, H_τ, q, χ, d_L` — the interpretation of the projection in physical units through the observer's scales.

Since the model is built in its own language, I provide a dictionary of correspondences between its terms and familiar physical notions. The terms are deliberately retained (they carry the model's ontology), but each is given a standard analogue:

| DSIC term | what it is | closest familiar image |
|---|---|---|
| objectness (`O`) | the pole of "thingness", the scale of objects | matter density / matter |
| space (`S`) | the pole of "extension", the scale of space | volume / spatial scale |
| norm `R` | the conserved sum of the poles `S² + O² = R²` | a conservation law (the link between the two poles) |
| flow `φ` | the single monotonic degree of freedom | time as the evolution parameter |
| the Mush | the limit `μ → 1`: maximum objectness, minimum space | a regular analogue close to the initial singularity |
| object bounce | the limit `μ → 0`: maximum space | de Sitter asymptotics |
| eversion | a sign-phase change in the core (a pole passing through its zero) | a mirror phase, unobservable from inside |
| bounce | the observable reversal at a limit | an artifact of the projection (of taking the modulus), not an event |
| projection (`μ`, `P`, `Q`) | modulus-valued, observable quantities | what is accessible to measurement |
| binding (Part II) | a measure of the cohesion of objects | an analogue of gravitational binding |
| detachment (`μ_d`) | the threshold beyond which monolithic matter lets enough space seep through for light to propagate freely | the epoch after which the core geometry works without correction |

The key thing to keep in mind: the "bounce" and the "Mush" are not physical events but the way the internal observer sees the smooth monotonic motion of the core through the modulus operation (§4.3). Beneath the observed oscillation between two limits lies a single forward motion in which nothing ever returns.

The boundary between the layers is absolute: no sign passes from the core into the projection (the modulus cuts it off), and no projection quantity enters the equations of the core.

Plan of the paper. Part I — the core: §2 states the postulates; §3 constructs the core (the circle, inversions, sign phases); §4 — the projection for the internal observer; §5 derives the observable physics; §6 confronts the model with data; §7 lists falsifiable predictions; §8 outlines the boundaries; §9 positions the model relative to ΛCDM. Part II — the second tier: §10 introduces the binding of objects and the detachment postulate; §11 constructs the early Universe (the threshold `μ_d`); Appendix B — the space-distribution module (gravity). §13 — general conclusion. Appendix A contains a summary of notation.

---
---

# Part I — The core and observable cosmology

---

## 2. Postulates

**P1 (relativity of scale).** No absolute scale exists. Only the ratio of scales carries physical content; any construction that opens access to an absolute scale is forbidden.

**P2 (two poles, one norm).** Objectness and spatiality are not two independent entities but two poles of a single relation, linked by a common degree of freedom. An object is a local convexity: it has an interior and an outer boundary, and occupies a place relative to other objects (a ball, a stone, a planet, an atom — variants of one structure). Space is the antipode of the object, a universal concavity: it is located nowhere and is contained in nothing; the very category of "where" arises only within it. The scales of the poles `S` (space) and `O` (objectness) obey the norm conservation law:

```
S² + O² = R² = const .
```

Growth of one pole is possible only at the expense of the other; the physics is carried by the fixedness of the norm under the flow between the poles, not by the absolute values of the poles.

**P3 (single flow).** The system has one fundamental degree of freedom — a continuous monotonic flow `φ` that never reverses. Its dynamics is free: the action contains no potential, no forces, no cosmological constant,

```
𝒮[φ] = ∫ dt · ½ φ̇² ,     δ𝒮 = 0  ⇒  φ̈ = 0  ⇒  φ(t) = ω·t ,   ω > 0 .
```

`φ` has no special points; everything else in the model is a function of `φ`. Cyclicity, if observed, is born not of the flow but of the way it is observed.

**P4 (internal observer).** Any observer is inside the system and is itself an object: its instruments transform along with it. Neither the sign of quantities nor absolute scale is accessible to it; only moduli are — in practice, only the approach and recession of objects.

**P5 (constancy of the speed of light).** In the units of the internal observer, `c = const`. As shown in §5.1, together with P1 this uniquely fixes the rate of the observer's clock.
(This is a conventionally fixed exchange rate linking the observer's ruler and clock — that is, one more ratio, not an absolute quantity. The numerical value of `c` enters nowhere in the derivation of the geometry and appears only when translating dimensionless predictions into observational units.)

---

## 3. The core: circle, inversions, sign phases

### 3.1 The carrier and the circle

Both poles are built from a single smooth carrier, shifted by a quarter turn (quadrature):

```
S(φ) = R·sin(φ/2) ,      O(φ) = R·cos(φ/2) ,      S² + O² = R² .
```

The carrier has a period of `4π` in the flow phase — a full sign cycle. The full description has definite parity and a constant Wronskian:

```
S(−φ) = −S(φ) ,    O(−φ) = O(φ) ,
W = S·(dO/dφ) − O·(dS/dφ) = −R²/2 = const .
```

The constancy of `W` is the formal record of the quadrature: when one pole is at an extremum, the other passes through zero; every discontinuous derivative of a modulus (arising later, in the projection) is multiplied by a paired factor that vanishes, so the product in the full signed description remains smooth.

If `S` and `O` are plotted along two axes, the invariant draws a circle of radius `R`: every admissible state of the system is a point on it. Evolution is the sliding of this point at constant angular velocity (the uniform flow `φ`); leaving the circle inward or outward is impossible (the norm is fixed), and when passing a point on an axis, the state does not stop or bounce — it keeps moving in the same direction. The rate of the flow between the poles, however, is not constant:

```
dμ/dφ ∝ sin(φ/2)   — slow near the poles, fast in the middle of the arc.
```

This nonlinearity is a purely geometric property of projecting uniform circular motion onto the axes, not the action of a force: the flow itself remains strictly uniform (P3).

### 3.2 Inversions

Along the circle, two events strictly alternate — smooth zero crossings of the poles:

```
space inversion:   S(φ) = 0  ⇒  φ = 2kπ        (k ∈ ℤ) ,
object inversion:  O(φ) = 0  ⇒  φ = π + 2kπ .
```

The zero of `S` coincides with the maximum `|O| = R`: this is the state of the **Mush** — the limiting density of objectness, in which almost no space remains and objects are compressed to the limit. An important distinction: the Mush is the entire zone near the zero of space, **on both sides** of it; the half before the zero is the finale of contraction, the half after it is the beginning of expansion. Only the second, outgoing half — where the modulus of space grows again and objectness is released — properly deserves the name **Big Bang**. One and the same zero of space separates the finale of the previous contraction from the start of the new expansion.

The zero of `O` coincides with the maximum `|S| = R`: this is the minimum of objectness — the **object bounce** in the observer's terms (see §4).

Because of the nonlinearity of the flow between the poles (§3.1), the intervals between inversions in `φ` are uneven, but the order of alternation is strict: Mush → object bounce → next Mush.

### 3.3 Sign phases and inheritance

Having passed its zero, a pole enters a **mirror (sign) phase**: it continues the same motion with the opposite sign — an eversion, a turning of the coordinates inside out, smooth at the level of the carrier. Denote the sign phases `σ_S = sign S(φ)`, `σ_O = sign O(φ)`.

The key structural asymmetry of the model: objectness is a **child** of space, and the two zeros are arranged differently.

- At `O = 0` (the object bounce), **only objectness** passes through its own zero — it is what everts, changing its sign phase. Space at that moment is at its maximum and undergoes only an observable bounce: its scale stops growing and begins to shrink, but the sign phase of space does not change.
- At `S = 0` (the Mush), **space itself** passes through its own zero — the common container everts, and the contents evert together with it: objectness, being at its maximum and not passing through its own zero, **inherits** the sign inversion forcibly, in sync with its parent — including the observer's measuring instruments.

The events are not mirror images of each other: at the zero of objectness one pole changes sign; at the zero of space, both do (one by itself, the other by inheritance). Formally, inheritance is expressed by the quantity

```
O_full(φ) = O(φ)·σ_S(φ) ,
```

which is continuous at the object-inversion points (`O` passes through zero smoothly, `σ_S` is constant) and discontinuous at the space-inversion points (an instantaneous sign flip). The discontinuity is unobservable from inside: `|O_full| = |O|` is continuous everywhere, and the synchronous sign change of the system and the instruments cancels out. It is precisely the continuity/discontinuity of `O_full` that will determine, in §5.3, the observed asymmetry of the two scale limits.

A terminological convention the paper adheres to throughout: **eversion** is a sign-phase change at the level of the core (it happens to the pole passing through its own zero); **bounce** is what the internal observer sees: a touching of the limit and a reversal. These are not synonyms but descriptions of one mechanism from two vantage points.

| zero | what everts (changes sign phase) | the other pole | the observer's "V" kink |
|---|---|---|---|
| `O = 0` (object bounce) | objectness itself | space at its peak: a mere bounce, sign unchanged | in `Q = \|O\|` |
| `S = 0` (Mush) | space + objectness by inheritance | objectness at its peak, everted forcibly | in `P = \|S\|` |

### 3.4 Falling into itself

The universal concavity of space (P2) also sets the character of the "collapse" of objectness: falling, in this picture, is not motion toward a point in space but the convergence of an object's scale toward its own center, proceeding inward from everywhere at once, with no preferred direction. Objects do not fly into a common pit "somewhere out there" — they fall into themselves. The gravitational development of this idea (a common center of infall for bound groups) is discussed in §8 as a direction not derived from the core.

---

## 4. The projection: the world of the internal observer

### 4.1 Two observers

One and the same event looks different from two vantage points. The **external observer** is an imaginary point of view outside space (in reality there is none: space is fundamental and has no "outside", but it is convenient for the full description). It sees the flow whole, with its sign: the history is continuous and monotonic, quantities pass through zeros and keep changing in the same direction; no reversals whatsoever. The **internal observer** is anyone inside. By P4, only moduli are accessible to it: where the external observer sees a smooth zero crossing, the internal one sees a touching of zero and a reversal — a bounce. Both the quiet eversion inside out and the smooth zero crossing read identically from within — as a change in the regime of approach/recession of objects, observed by a being that is itself made of what is changing.

### 4.2 The modulus operation

The transition from the core to the observable is performed by a single operation — taking the modulus:

```
π : (S, O) ↦ (|S|, |O|) .
```

It does three things: it cuts off the sign (the mirror branch is inaccessible); it folds a smooth zero crossing into a kink — a bounce; and it identifies adjacent sign branches, whereby the observable period turns out to be half the sign period. The observable quantities:

```
P(φ) = |S| = R·|sin(φ/2)| ≥ 0 ,      Q(φ) = |O| = R·|cos(φ/2)| ≥ 0 ,
P² + Q² = R² ,

μ(φ) = Q/R = |cos(φ/2)| ∈ [0, 1] ,      P = R·√(1 − μ²) ,   Q = R·μ .
```

The observable scale phase `μ` (the fraction of objectness) has period `2π` versus the sign period `4π`: two events distinct in the full description merge, for the observer, into one.

```
full description (signed):  … +O → 0 → −O …   (monotonic crossing)
observable (modulus):       …  |O| → 0 → |O| …  (bounce)
```

The pair `(P, Q)` is nonnegative and always lies in the first quadrant. Over one observable cycle the point traverses this quadrant there and back: from the Mush (`P = 0, Q = R`) to the object bounce (`P = R, Q = 0`) and back. The back-and-forth motion arises not from a superposition of arcs but from the symmetry of `Q = R·|cos(φ/2)|` about the bounce point: on reaching `Q = 0`, the point reflects and returns along the same arc. The monotonic circular sliding of the core looks from inside like an oscillation between two limits; the number of sign beats of the core is thereby unobservable — a property of the layer boundary, not an indeterminacy of the model.

### 4.3 The bounce is not a physical event

The carrier is smooth and has no vertices of its own, so **both** observed kinks of `μ` have the same nature — they are created by the modulus as the corresponding pole passes through zero: at `O = 0`, `Q = |O|` folds into a "V" kink (at this point `P` smoothly passes its maximum); at `S = 0`, so does `P = |S|` (the maximum of `Q` is smooth). The bounce requires no force, no braking, no cause: there is no reversal at the level of the mechanism at all — only a change in the way of seeing. "Passing through zero" and "bounce" are one event from two vantage points.

Over one cycle the internal observer sees four regimes: (i) the recession of objects reaches its limit (`μ` is minimal) — the object scale begins to grow; (ii) objects approach one another, space shrinks; (iii) the Mush — objects collide and acquire a superdense state, squeezed tight by space; (iv) the Mush disperses (Big Bang) — objects recede. All the turns look like changes of regime of the same kind, but it is the projection that repeats, not the flow: beneath it lies a single monotonic forward motion in which nothing ever returns.

---

## 5. Observable physics

Everything in this section is built from `μ` and the internal observer's scales; it is an interpretation of the projection in physical units, not a separate mechanism.

### 5.1 The scale factor

The observer is made of objectness, so its unit of length is tied to the object pole: the ruler `L(μ) ∝ μ`. The metric is the ratio of the poles `P/Q = √(1 − μ²)/μ`; near the zero of objectness the ruler contracts to a point, and relative measurements diverge.

The clock rate is fixed by postulate P5. The period of the proper clock is `T(μ) ∝ μᵖ`; then the speed of light from inside is `c = L/T ∝ μ^{1−p}`. The requirement `c = const` yields the unique value `p = 1`, i.e. the observer's proper time

```
dτ = dt / μ .
```

Any other exponent `p` would drive `c` across a range and open access to the absolute scale forbidden by P1. The observable scale is composed of the metric and the clock rate:

```
b(μ) = (P/Q)·(1/μ) = √(1 − μ²) / μ² .
```

This is the central formula of the model: one function of one variable, with no free shape parameters.

### 5.2 Redshift

Light **propagates through space** (the pole `P`) but is **measured** by an observer built of objectness (`L ∝ μ`, `dτ = dt/μ`, `c = const`). Propagation and measurement run along different poles, so both enter the observed wavelength. At an epoch with phase `μ` it is composed of three contributions, each entering exactly once:

```
(1) stretching of the path through space:  ∝ √(1 − μ²)    [pole P]
(2) the source's clock, dτ = dt/μ:         ∝ 1/μ          [time]
(3) the object ruler, L ∝ μ:               ∝ 1/μ          [length]

λ_obs(μ) ∝ √(1 − μ²) / μ² = b(μ) .
```

The observed wavelength is identically equal to the scale factor, and the redshift falls out of the accounting of contributions automatically — it is derived, not introduced:

```
1 + z = λ_obs(μ₀)/λ_obs(μ_e) = b(μ₀)/b(μ_e)
      = (μ_e²/μ₀²) · √(1 − μ₀²)/√(1 − μ_e²) ,
```

where `μ₀` is the phase of reception ("today") and `μ_e` the phase of emission. As `μ_e → 1` (the Mush), `z → ∞`; as `μ_e → μ₀`, `z → 0`. The formula admits a decomposition into the two poles:

```
1 + z = [ √(1 − μ₀²)/√(1 − μ_e²) ] · [ μ_e²/μ₀² ]
       └───── pole P (space) ─────┘ └ object scale ┘
```

— the first factor accounts for the stretching of space along the path, the second for the observer's object scale (clock + ruler). Their alternating dominance along the cycle is precisely the mechanism of the observed kinematics (§5.4). Note also the physics of the limit: near the Mush the pole `P` contracts to almost nothing — light simply has almost nothing left to propagate through; the path through space collapses, and the observed picture rests on the object scale alone. Light does not "slow down" — the stage it travelled on vanishes.

The alternative — a "redshift by a single ruler", `1 + z = μ₀/μ_e` — is rejected: it accounts only for contribution (3) and ignores the path (1) and the clock (2), whereas the photon lives through all three. The observer's full metric is `b`, so the redshift goes by `b`, not by `L`.

Two properties of the formula. First, `b(μ)` is strictly decreasing on `(0, 1)`, so the correspondence `z ↔ μ_e` on the physical reception branch `μ_e ∈ (μ₀, 1)` is one-to-one and the inversion `z → μ_e` is well defined. Second, only the modulus quantity `μ` enters `1 + z` — no sign from the core is involved, and the layer boundary is not violated. By itself, `z` in isolation is unobservable: it enters the data only through `d_L(z)`.

### 5.3 The limits of the cycle and their asymmetry

```
μ → 0  (object bounce):  b → ∞ ,
μ → 1  (Mush):           b → 0 .
```

At the object bounce, space is maximal, objects are vanishingly small relative to it, and the observable scale diverges; this is **a pole of approach, not an attainable singularity** — the observer slides past `μ = 0` without ever landing on it. In the Mush, no ratio of scales remains at all. The asymmetry `∞ vs 0` does not contradict the symmetry of the core but follows directly from the mechanics of sign (§3.3): objectness passes through its own zero by itself and smoothly (`O_full` continuous → divergence), whereas through the zero of space it is the container that passes, everting the sign of its contents (`O_full` discontinuous → vanishing). The divergence of `b` near the zero of objectness is a direct expression of the relativity of scale (bounded forms of `b` fail postulate P1) and an artifact of the projection, not an event of the core.

### 5.4 Kinematics

Introduce the branch direction `s = ±1` (`s = −1` — expansion, `μ` decreasing; `s = +1` — contraction) and the characteristic phase speed `k = (2/π)·ω > 0`, so that along a branch `dμ/dt = s·k = const`; the factor `2/π` fixes the time scale so that the observable cycle in `μ` has period `2π`. All derivatives are taken with respect to proper time `dτ = dt/μ`, i.e. `d/dτ = μ·d/dt`. The constancy of `k` is the definition of a branch, not an approximation: with a non-constant rate, a term `∝ dk/dt` would enter the deceleration parameter, and here it is identically zero. The kinematics therefore depends only on the phase `μ`, not on how fast the system traverses it.

The expansion rate:

```
H_τ(μ) = (db/dτ)/b = s·k·(μ² − 2)/(1 − μ²) ,
|H_τ(μ)| = (2 − μ²)/(1 − μ²)                    [units of k] ;
|H_τ| → 2 as μ → 0 ,     |H_τ| → ∞ as μ → 1 .
```

The deceleration parameter:

```
q(μ) = −(d²b/dτ²)·b/(db/dτ)² = (−μ⁴ + 6μ² − 4)/(μ⁴ − 4μ² + 4) .
```

It depends on the direction `s` only through `s² = 1`, and not on the speed `k` at all (the factor `k²s²` cancels): `q` is a function of the phase alone. The key property: `q(μ)` **changes sign**, with a crossover at `μ_cross ≈ 0.874`:

| state | `μ` | `q` | `\|H_τ\|` (units of `k`) |
|---|---|---|---|
| Mush | `→ 1` | `+1` | `→ ∞` |
| transition | `0.874` | `0` | `≈ 5.2` |
| today | `≈ 0.76` | `−0.43` | `≈ 3.4` |
| midpoint | `0.5` | `−0.84` | `2.33` |
| object bounce | `→ 0` | `−1` (de Sitter) | `2` |

The deceleration → acceleration transition arises from the very curvature of `b(μ)`, without an introduced constant, and its mechanism can be read off the two-pole decomposition of the redshift (§5.2). Near the Mush (`μ_e → 1`) the spatial factor dominates: the expansion is set by the stretching of space seeping through the monolithic objectness — the deceleration regime (`q → +1`). In the late Universe (`μ_e → μ₀`) the observer's object scale dominates — the acceleration regime (`q₀ ≈ −0.43` at `μ₀ ≈ 0.76`), with the de Sitter asymptotics `q → −1` as `μ → 0`. Acceleration here is not the result of something "pushing" the Universe: the observer's object ruler and the scale of space drift apart ever faster, and their relative drift is read from inside as accelerated expansion. Acceleration is the shadow of the drift, not a force; de Sitter is the asymptotics of the approach to the object bounce, not an introduced constant.

The clock factor `dτ = dt/μ` is not cosmetic: it changes the kinematics substantially across the whole range, not only near the edges. In the middle of the range (`μ = 0.5`), the deceleration parameter without the clock is `≈ −1.27` versus `≈ −0.84` with it.[^noclock]

[^noclock]: The formula without the clock factor (derivatives taken directly with respect to `dt`, where `dμ/dt = s·k`): `q_noclock(μ) = (−2μ⁴ + 9μ² − 6)/(μ⁴ − 4μ² + 4)`; at `μ = 0.5` it gives `−1.27` versus `−0.84` from the full formula of §5.4.

### 5.5 Distances

Light travels with `c = const` in the local units of each epoch; the comoving distance accumulates along the path as

```
dχ = (c/k) · [b(μ₀)/b(μ)] · dμ/μ ,
```

and with `b(μ) = √(1 − μ²)/μ²` the integral is taken in closed form (`∫ μ/√(1 − μ²) dμ = −√(1 − μ²)`):

```
χ(μ_e; μ₀) = (c/k) · [√(1 − μ₀²)/μ₀²] · [ √(1 − μ₀²) − √(1 − μ_e²) ] ,

d_L(z) = (1 + z)·χ(μ_e; μ₀) ,        1 + z = b(μ₀)/b(μ_e) ,
DM(z) = 5·log₁₀(d_L / 10 pc) .
```

Two consequences. **The particle horizon is finite by itself**, without a cutoff: at the Mush

```
χ(μ_e → 1) = (c/k)·(1 − μ₀²)/μ₀²   — finite.
```

**Self-consistency:** the identity `D_H = dχ/dz` (where `D_H = c/H`) holds exactly — the rate `H(z)` of §5.4 and the distances `χ(z)` are tied by one geometry, and the distances computed both ways coincide.

---

### 5.6 Translation into the language of the Friedmann metric

The DSIC core can also be read in the standard language of general relativity — without changing a single formula, merely by matching quantities. The observable scale factor `b(μ)` plays the role of the scale factor `a` in a flat Friedmann–Robertson–Walker metric, with the observer's proper time `dτ = dt/μ` serving as cosmic time:

```
ds² = −c² dτ² + b(μ)² · dx² .
```

This is an identical translation: `b(μ)` is not an analogue but literally the scale factor of such a metric. The question "what effective medium does the DSIC geometry imitate" is then settled by the standard identity of flat FRW cosmology linking the equation of state `w_eff = p/ρ` to the deceleration parameter,

```
w_eff = (2q − 1)/3 ,
```

and, substituting `q(μ)` from §5.4, yields the closed expression

```
w_eff(μ) = (−μ⁴ + (16/3)μ² − 4) / (μ⁴ − 4μ² + 4) .
```

Not a single new entity is introduced here: `w_eff` is entirely determined by the same `q(μ)` already derived from the core. The limiting values show that DSIC, in GR language, is a flat Universe with a single medium whose effective equation of state glides smoothly from radiation to a cosmological constant:

| state | `μ` | `w_eff` | GR interpretation |
|---|---|---|---|
| Mush | `→ 1` | `+1/3` | radiation-like (deceleration) |
| transition | `0.874` | `−1/3` | the acceleration boundary (`q = 0`) |
| today | `≈ 0.76` | `−0.62` | dark-energy-like (not a constant) |
| midpoint | `0.5` | `−0.89` | close to de Sitter |
| object bounce | `→ 0` | `−1` | exactly a cosmological constant |

Let me stress the meaning of this translation. DSIC does **not postulate** a medium with such a `w(z)` — it obtains this curve from the pure geometry of the flow of scale; `w_eff(μ)` is merely its shadow in GR language. Where ΛCDM inserts `Ω_Λ` by hand, DSIC reproduces the equivalent dynamics with a single geometric parameter `μ₀`. The translation does not strip the model of its ontology — it merely makes it commensurable with the standard formalism and shows that DSIC does not conflict with general relativity at the level of the background metric, but fits into it as a universe with a dynamical `w(z)`.

> Remark. Agreement of the background kinematics with some FRW model does not imply agreement at the level of perturbations and structure growth: a dynamical `w(z)` fixes the background but does not uniquely determine the evolution of inhomogeneities (see §12). The translation concerns precisely the background.

---

## 6. Confrontation with observations (`z < 2.3`)

The form of all observable dependences is rigidly fixed by the function `b(μ) = √(1−μ²)/μ²`; analogues of `Ω_m, Ω_Λ` are absent. Free are the initial condition `μ₀` (the phase "today") and the scale `H₀` (for supernovae, also the absolute magnitude `M`, as in ΛCDM).

On scales. BAO, cosmic chronometers, and the age are computed at `μ₀ = 0.76` on the early scale `H₀ = 67.4`; the conversion of the geometry into megaparsecs goes through `c/k = (c/H₀)·(2−μ₀²)/(1−μ₀²) ≈ 14978 Mpc` on the same scale. A direct fit to the supernovae (§6.1) yields `H₀ ≈ 73` — this is the familiar "Hubble tension" [8], reproduced identically in DSIC and ΛCDM; the fitted value is not substituted into the early scale, since these are different scales, kept apart by the tension itself. The model creates no new tension and removes no existing one.

### 6.1 Supernovae: Pantheon+SH0ES

The model is confronted with the full **Pantheon+SH0ES** catalog (1701 SNe; 1657 enter the fit, including 77 calibrators [5] as the anchor of the absolute scale) with the full STAT+SYS covariance matrix. Three parameters are free (`μ₀, H₀, M`); flat ΛCDM is fitted to the same data with likewise three (`Ω_m, H₀, M`):

```
DSIC : μ₀ = 0.7600 ± 0.0091 ,  H₀ = 73.35 ± 1.01 ,  M = −19.243 ± 0.030 ,  χ²/dof = 0.8783
ΛCDM : Ω_m = 0.3308 ± 0.0180 , H₀ = 73.53 ± 1.02 ,  M = −19.244 ± 0.030 ,  χ²/dof = 0.8783

Δχ² (DSIC − ΛCDM) = +0.05        (ΔAIC = ΔBIC = +0.05)
```

The uncertainties are 1σ from the Hessian of `χ²` at the minimum (`C = 2H⁻¹`); the model parameter is measured to ~1.2%. A separate pipeline validation: the ΛCDM branch of the same pipeline reproduces the published values of the Pantheon+SH0ES analysis (`Ω_m ≈ 0.33 ± 0.018`, `H₀ ≈ 73.5 ± 1.0`) — both models are computed by one code on one dataset, and the agreement of the control branch with the literature confirms the correctness of the setup.

The difference in fit quality is five hundredths: on fully calibrated data, DSIC is **statistically indistinguishable** from ΛCDM (`|Δχ²| ≪ 2`). The one-parameter geometry of `μ₀` reproduces the same distance–redshift relation that the standard model sets with two density parameters — without dark energy. The Hubble constant thereby emerges as a genuine inference from the anchor (`H₀ ≈ 73`), not as an absorber of a shift in `M`. A run with the STAT-ONLY covariance gives `μ₀ = 0.7661`, `Δχ² = −0.39` — the result is robust to the treatment of errors.

### 6.2 The baryon scale: DESI DR2

The DESI DR2 data (Table IV [6], `0.5 < z < 2.33`) give `D_M/r_d` and `D_H/r_d` separately in six bins. The dimensionless ratio `D_M/D_H` depends neither on `r_d` nor on `H₀` and tests the pure geometry of the expansion:

```
D_M/D_H (DSIC vs DESI DR2, 6 points):  χ²/dof ≈ 0.78
```

— on a par with the standard model on the same points. The `z = 0.934` point gives the largest tension in ΛCDM as well (≈ 2.5σ) — a feature of the data, not of the model.

The scale `r_d` is recovered as an inverse problem: `r_d = D_M^DSIC/(D_M/r_d)_DESI` and, independently, `r_d = D_H^DSIC/(D_H/r_d)_DESI` in each bin. The result is constant in `z` and consistent between the transverse and radial measurements:

```
r_d ≈ 148 Mpc    (scatter over z ≈ 1.4%; from D_M ≈ 147.7, from D_H ≈ 150.2) ,
```

which coincides with the standard scale `≈ 147 Mpc`. There is no dynamical `r_d(z)` in the data: the poles `P` and `Q` change with `z`, but the recovered physical scale stays constant — the DSIC ruler is static in the observed window.

### 6.3 Expansion rate and age

Cosmic chronometers [9] `H(z)` at `0.07 < z < 2` give `χ²/dof ≈ 0.37`, on a par with the standard model. The integral over the expansion history gives `t₀·H₀ ≈ 0.925`, i.e. at `H₀ = 67.4` an age `t₀ ≈ 13.4 Gyr`. A caveat on the comparison: the oft-quoted `13.80 ± 0.02 Gyr` is an age *derived within ΛCDM* from the Planck data, and measuring a non-ΛCDM model against it is partly circular; the value `t₀ ≈ 13.4 Gyr` passes the model-independent stellar ages (the oldest globular clusters, `~13.5 ± 0.5 Gyr`).

### 6.4 The deceleration parameter and the transition redshift

Through the relation `μ(z)` (§5.2) at `μ₀ ≈ 0.76`, the theoretical curve `q(z)` is compared with model-independent kinematic reconstructions [12]:

| `z` | `q_DSIC` | reconstruction (model-indep.) |
|---|---|---|
| `0` | `−0.43` | `−0.55` |
| `0.3` | `−0.25` | `−0.34` |
| `0.65` | `−0.06` | `0.0` |
| `1.0` | `+0.11` | `+0.18` |
| `1.5` | `+0.30` | `+0.33` |
| `2.0` | `+0.44` | `+0.42` |

The shape of `q(z)` — S-shaped, with deceleration at high `z` and acceleration at low `z` — reproduces the reconstructions within their scatter; the transition redshift `z_t ≈ 0.77` [10] agrees with the observed `z_t ≈ 0.6–0.8`.

### 6.5 Joint fit and `Om(z)`

The joint SN+BAO fit formally gives a weak preference to DSIC (`Δχ² ≈ −5.0`); the advantage partly reflects the known DESI–ΛCDM tension and is not interpreted as decisive. The `Om(z)` diagnostic **does not discriminate** between the models on current data: the sign of `Δχ²` depends on the `H₀` normalization. This is a null result.

### 6.6 Summary of the confrontation

A single initial condition `μ₀ ≈ 0.76` consistently describes six independent probes of the late Universe:

| probe | data | result |
|---|---|---|
| supernovae, `d_L(z)` | Pantheon+SH0ES (1657 SNe, anchor + covariance) | `Δχ² = +0.05` vs ΛCDM — indistinguishable |
| BAO `D_M/D_H` | DESI DR2, 6 points | `χ²/dof ≈ 0.78` |
| the `r_d` ruler | DESI DR2 `D_M/r_d`, `D_H/r_d` | `≈ 148 Mpc`, static |
| chronometers `H(z)` | `z < 2` | `χ²/dof ≈ 0.37` |
| age `t₀` | integral over the history | `≈ 13.4 Gyr` (stellar ages `~13.5 ± 0.5`) |
| transition `z_t` | `q(z)` reconstructions | `≈ 0.77` (obs. `0.6–0.8`) |

The observer's position on the cycle is fixed by these observations consistently: one geometry, one parameter, no dark energy.

---

## 7. Distinguishing predictions and falsifiability

On late-time data DSIC and ΛCDM are practically degenerate; what separates them is high `z`, where the model makes sharp, testable predictions:

- **the shape of `Om(z)`** over an extended redshift range;
- **the radial BAO in the Lyman-α region** (`z ≈ 3.5`): a departure from ΛCDM of up to ~5%;
- **the expansion rate `E(z)` at `z = 5–10`**: a departure of +15…+47%.

It is important to separate the sources of the signals: the `D_H` departure in the window `2.33 < z < 4.53` is a trace of the **core** (`D_H = c/H` depends only on `μ₀`) and is independent of the Part II machinery; the trace of the threshold `μ_d` lies in `D_M` above `z ≈ 4.5` and is discussed in §11.5.

---

## 8. Boundaries of applicability and open directions

The claims are deliberately kept narrow: only the background kinematics of the late Universe (`z < 2.3`) of §6 is considered established. The following has not been done.

**Structure growth.** Perturbations are not derived from the core, but the DSIC background has been checked for compatibility with the standard (GR-limit) growth of linear inhomogeneities: on the `fσ₈` compilation the model describes the data no worse than ΛCDM (§12). A full derivation of the growth equation from binding (Appendix B) remains open.

**The early Universe.** From the core alone, the distance to the CMB falls short by ~21%. This is closed by the model's second tier (Part II): the form of the correction is derived from the postulate on the binding of objects (P6, §10), and its single threshold `μ_d` is calibrated against a single observation. What this resolves, what it does not, and why — §11.6.

**Ontology versus ansatz.** The form `b(μ) = √(1−μ²)/μ²` is derived here from the circular geometry of the poles, yet the reader is free to accept it simply as a fortunate one-parameter ansatz: the formulas work on their own, and the ontological picture is not required for them to work.

**Gravity.** The local distribution of the seeping space requires the same entity the core lacks — the binding of objects — and is treated in Appendix B as the late manifestation of the second tier; this module has not yet been confronted with data.

---

## 9. Relation to ΛCDM

DSIC does not "refute" ΛCDM. On late-time data the models are statistically indistinguishable; DSIC offers a different **ontology** of the same background — geometry instead of dark energy — with an equal number of free parameters. The value lies not in the quality of the fit (it is no better) but in the origin of the curve: the same distance–redshift relation is obtained from a single geometric principle, without `Ω_Λ`, with the deceleration → acceleration transition and the finiteness of the particle horizon derived rather than postulated. The decisive test is the high-`z` measurements of §7.

---
---

# Part II — The second tier of the model: the binding of objects

> **Status.** Everything in this part rests on the core but is not derived from it: the poles `S, O` contain nothing but the flow of scale, while both constructions below require a second, independent entity — the **binding of objects**. The core (Part I) is self-sufficient: all six probes of §6 are independent of Part II, and the second-tier correction is identically zero throughout the domain of their data. Part II introduces one new entity, one additional postulate, and one measured constant; where derivation ends and calibration begins is stated at every step.

---

## 10. The second entity and the detachment postulate

The core runs up against one and the same missing entity from two sides. From the side of the early Universe: the core geometry falls short of the distance to the CMB by ~21% (§11.1) — what is missing is the physics of the phase in which matter is still stuck together, and "being stuck together" is a state of binding. From the late side: the qualitative picture in which space accrues predominantly between objects weakly bound by gravity (Appendix B) requires a measure of that bond. Both superstructures — the early detachment and the late distribution of space — are two manifestations of one new entity: **binding**. I introduce it once.

**P6 (postulate of binding and detachment).**

- The early epoch is nearly monolithic objectness, and the threshold refers to the **average medium along the light path**, not to every point. Until the dominant monolith has, on average, parted (`μ > μ_d`), the internal observer's measurement chain along the path is undefined: the ruler `L ∝ μ` is assembled from separate objects, while the medium is still ruled by the stuck-together component. Light in this phase has almost nowhere to propagate — the path is counted not by the observer's metric but by the sheer availability of space: the addition to the path is `∝ dμ/P = dμ/√(1−μ²)`. The weighting `b(μ₀)/b(μ)` that weights the path in the late formula (§5.5) is not applied here — there is nothing to weight with.
- Local clumps in which binding is strongest **fall into themselves** (§3.4) and collapse into objects long before the threshold — the first galaxies and quasars. Their appearance does not change the average availability of space along the light path and does not shift the threshold: `μ_d` marks the moment when the monolith parted on average, not the moment the first objects were born.
- The detachment of the medium is a mechanical event, not a process: a **hard threshold** `μ_d`. Once matter has admitted space to the degree at which it could detach from itself and part, light gains a growing ability to propagate: at `μ ≤ μ_d` the correction is identically zero, and the late core geometry operates unchanged.
- The form is fixed by principle: the amplitude is the natural one (`A ≡ 1`, the form `1/P` without a prefactor), the window is sharp (no smoothing). Exactly one number remains free — the threshold `μ_d`.

The distinction "average medium vs. first islets" is not cosmetics but a necessary consistency condition: JWST galaxies are observed out to `z ≈ 10–14` (`μ_e(z=14) ≈ 0.9972`), i.e. **above** the threshold `μ_d = 0.9806` (`z_detach ≈ 4.53`). A literal reading — "no separate objects exist at `μ > μ_d`" — would contradict these observations; in the adopted formulation, the first objects are clumps that have fallen into themselves inside the monolith where binding is maximal, and their existence is compatible with a threshold of the average medium. The verification script prints this consistency flag explicitly.

The physical picture is the analogy of a press. A set of items has been compressed; when the press releases, the items at first still push against one another, and at some moment they **detach** — after which there is no more pressure on the walls, though the press keeps opening. The detachment point is precisely `μ_d` — the moment the mass as a whole parts; individual items deep inside may have latched onto one another earlier, which does not affect the threshold. The sharpness of the threshold is motivated precisely by the mechanics: detachment is an event; a smooth window would introduce an extra parameter with no physical basis.

A detail supporting the naturalness of the form. The basic path integrand of the core is `(√(1−μ₀²)/μ₀²)·μ/P` (§5.5), the correction integrand is `1/P`; right at the Mush (`μ → 1`) they coincide in form: both are `∝ 1/P`. The path made up in the stuck-together phase is asymptotically commensurate with the basic path measure — light "picks its way" through objectness too, not only through the space that has already seeped through.

**The status of `μ_d`.** The threshold's value is not derivable from the core: the poles `S, O` contain no number that singles out `0.9806` — the threshold is a condition on binding, which the core does not have. But with the form fixed by P6, the "amplitude ↔ threshold" degeneracy is broken by principle, and `μ_d` becomes a **measured constant of the second tier** — of the same status as `μ₀` in the core. A theory is permitted to have measured constants; what it is not permitted to have is arbitrary functions. Without P6 the calibration would be degenerate, as is honestly shown by the table in §11.2.

---

## 11. The early Universe: the detachment threshold

### 11.1 The core's shortfall

The core parameters on the early scale are the same as in §6 for BAO, chronometers, and the age: `μ₀ = 0.76`, `H₀ = 67.4`, `c/k = 14978.2 Mpc`. The reference figures of Planck 2018 [7] (TT,TE,EE+lowE+lensing): `z* = 1089.92`, `r* = 144.43 Mpc`, `100θ* = 1.04110`, whence the comoving distance to the last-scattering surface `d_LSS = r*/θ* = 13872.8 ± 25 Mpc` (the internal consistency of `r*/θ*` is verified by a self-check of the script). The core geometry without the correction:

```
μ_cmb(z* = 1089.92) = 0.9999994681        (from 1+z = b(μ₀)/b(μ_e))
D_M base to the CMB  = 10936.2 Mpc
d_LSS (Planck)       = 13872.8 Mpc
shortfall            = 2936.6 Mpc  (21.2%)
```

The cause of the shortfall lies on the core's side: the particle horizon is finite "for free" (§5.5), almost the whole path is accumulated by `z ≈ 10`, beyond which the distance integral runs out of breath. The near and intermediate Universe (`z < 2.3`) is computed correctly all the while, so the detachment threshold affects only the far end.

### 11.2 The correction and the calibration of the threshold

By P6, the addition to the path in the stuck-together phase is `1/P`, cut off by the threshold; its integral is `arcsin μ`:

```
correction(μ) = 1/√(1−μ²)   for μ > μ_d ;      0   for μ ≤ μ_d

D_M_corr(μ_e) = D_M(μ_e) + (c/k)·[arcsin μ_e − arcsin μ_d]   for μ_e > μ_d
D_M_corr(μ_e) = D_M(μ_e)                                       for μ_e ≤ μ_d
```

The threshold is fixed by the requirement that the shortfall be closed exactly:

```
(c/k)·[arcsin μ_cmb − arcsin μ_d] = 2936.6 Mpc
⇒  μ_d = 0.980640 ,      z_detach = 4.526
D_M_corr(μ_cmb) = 13872.8 Mpc   (= d_LSS Planck, exact closure)
```

For honesty's sake, let me show what P6 resolves. Without the amplitude being fixed, the equation `A·(c/k)·[arcsin μ_cmb − arcsin μ_d] = shortfall` has a curve of solutions — any amplitude yields its own `μ_d` with the same distance to the CMB:

| `A` | `μ_d` | `z_detach` |
|---|---|---|
| 0.5 | 0.923706 | 1.506 |
| 0.8 | 0.969868 | 3.344 |
| **1.0** | **0.980640** | **4.526** ← fixed by P6 |
| 1.2 | 0.986514 | 5.691 |
| 2.0 | 0.995097 | 10.266 |
| 5.0 | 0.999190 | 26.922 |

One observation fixes one parameter: zero degrees of freedom. The value `z_detach ≈ 4.53` is a calibration with the form fixed by principle, not a prediction; it is not presented as a prediction.

### 11.3 Properties of the optical correction

The correction is **optical**, not dynamical: the light path is lengthened, the expansion rate `H(z)` is unchanged. Testable properties follow from this choice, and all of them hold here by construction:

- **Frequency is untouched.** The redshift remains `1+z = b(μ₀)/b(μ_e)` (§5.2); lengthening the path does not affect it.
- **Distance duality is preserved.** `d_L` and `d_A` are built from one and the same `χ_corr`, so the Etherington relation `d_L = (1+z)²·d_A` holds automatically.
- **Achromaticity and non-dissipativity are a hard requirement.** The blackbody nature of the CMB spectrum (FIRAS) rules out any "medium" that depends on frequency or scatters; the correction is obliged to be geometry of the path, not an interaction — and the form `1/P` satisfies this.
- **The particle horizon remains finite.** With the correction, `χ_corr(μ_e → 1) = (c/k)·[(1−μ₀²)/μ₀² + (π/2 − arcsin μ_d)] ≈ 13906 Mpc` — only ~33 Mpc beyond the last-scattering surface: in the corrected geometry the CMB lies almost flush against the particle horizon.

### 11.4 The late Universe is untouched

The correction is strictly zero below the threshold — an identity, not an approximation:

| `z` | `D_M` base, Mpc | correction, Mpc |
|---|---|---|
| 0.50 | 1932.8 | +0.0 |
| 1.00 | 3384.2 | +0.0 |
| 2.33 | 5792.7 | +0.0 |
| 4.00 | 7335.6 | +0.0 |
| 4.53 | 7653.3 | +0.0 |
| 4.63 | 7707.8 | +49.3 ← switched on |
| 5.00 | 7896.9 | +220.4 ← switched on |
| 10.00 | 9247.2 | +1433.0 ← switched on |

All six DESI DR2 bins (`z ≤ 2.33`) are insensitive to the correction, and the BAO `χ²/dof` does not change — the results of §6 are fully preserved.

### 11.5 The threshold's trace and the reach of the probes

The correction to `D_M` becomes appreciable only deep beyond the threshold:

| `z` | object | correction to `D_M` |
|---|---|---|
| 4.53 | detachment threshold | +0.0% |
| 5.00 | JWST galaxies | +2.8% |
| 6.00 | quasars | +7.1% |
| 8.00 | reionization | +12.4% |
| 11.00 | first galaxies | +16.6% |
| 1089.92 | CMB (calibration) | +26.9% |

The reach of the probes: DESI DR2 (`z ≤ 2.33`), DESI DR3 (~2027, `z ≤ 3.55`), and DESI-II/LBG (~2029, `2.5 < z < 3.5`) — all below the threshold. The `D_H/r_d` departure in the window `2.33 < z < 4.53` (up to ~5% at `z ≈ 3.5`, §7) is a trace of the core, not of `μ_d`. The threshold itself affects only `D_M` above `z ≈ 4.5`, where no precise standard rulers exist apart from the CMB calibration itself: no second independent anchor for `μ_d` exists in reachable data — until a precise distance scale appears at `z ≳ 5–6`.

### 11.6 What the second tier resolves and what it does not

Exactly one number is closed — the distance to the last-scattering surface, exactly and with no consequences for the late Universe. What remains open is the following, and this is a matter of principle.

**The borrowed sound horizon.** The closure target `d_LSS = r*/θ*` uses `r* = 144.43 Mpc` — a quantity not measured directly but computed within ΛCDM from the physics of the photon–baryon plasma, which DSIC does not have. The honest formulation: the directly measured quantity is the angle `θ*`, and it is against `θ*` that the calibration is made; the origin of the acoustic scale itself, `~147 Mpc` — showing up consistently in the late ruler as well, `r_d ≈ 148 Mpc` (§6.2) — is declared an open question of the core.

**The expansion rate is untouched — and this is the narrowest point, stated as a number.** The optical correction changes the path but not `H(z)`: at `z = 5–10` the rate remains +15…+47% above ΛCDM (§7), which means cosmic time at those epochs is compressed. Quantitatively (common scale `H₀ = 67.4`; `Ω_m = 0.303` from the joint fit of §6, whence the ΛCDM `t₀` here is 13.95 Gyr — both models are put on an equal footing; the DSIC `t₀` = 13.38 Gyr):

| `z` | `t_DSIC`, Gyr | `t_ΛCDM`, Gyr | ratio |
|---|---|---|---|
| 3.0 | 1.73 | 2.18 | 0.79 |
| 5.0 | 0.81 | 1.19 | 0.68 |
| 6.0 | 0.60 | 0.95 | 0.64 |
| 7.5 | 0.42 | 0.71 | 0.59 |
| 10.0 | 0.25 | 0.48 | 0.52 |
| 12.0 | 0.18 | 0.37 | 0.48 |

A Salpeter estimate sharpens the picture. At the Eddington limit, one e-folding of a black hole's mass takes ~45 Myr; growth from a stellar seed of `~10² M☉` to the masses `~10⁹ M☉` observed in quasars at `z ≈ 7.5` requires ~16 e-foldings ≈ 0.72 Gyr. ΛCDM has 0.71 Gyr at its disposal — just barely enough (which by itself already motivates heavy seeds); DSIC has 0.42 Gyr ≈ 9 e-foldings — either direct-collapse seeds of `~10⁵ M☉` with continuous accretion, or super-Eddington regimes, become mandatory. This is not a refutation — super-Eddington accretion is discussed in earnest, not least because of the "overmassive" black holes discovered by JWST — but it is the nearest real test of the core and probably the first point of serious resistance to the model; it is stated as a number, not as a caveat.

The core's picture also offers a natural, though not rigorously derivable, mechanism behind these observations. The detachment of the monolith (§10) is an event of the **average medium**, not of every point: when space seeps through and objectness as a whole parts, individual patches of limiting density may for various reasons **fail to part** and retain the density of the Mush — just as in any decompressing substance some clumps remain pressed together. Such regions that never admitted space are ready-made superheavy black hole seeds, existing **from the very start of the expansion** rather than grown by accretion from stellar remnants. In this picture, the "overmassive" black holes discovered by JWST at high `z` are not an anomaly demanding record-fast growth but an expected trace of inhomogeneous detachment: in some places the monolith parted, in others it remained Mush. A rigorous derivation of the distribution of such seeds is not provided by the core (it would require the statistics of binding, Appendix B), but their very existence does not contradict the picture — it follows from it, and the compressed time budget at `z = 5–12` ceases to be a pure vulnerability: the model expects a share of the massive early holes to be "innate".

**The CMB is not a single number.** Planck constrains the entire structure of the acoustic peaks, the damping tail, and the polarization — thousands of multipoles. The second tier closes one parameter (`θ*`); everything else requires a theory of perturbations and a composition (baryons, photons, neutrinos), which the model does not have.

### 11.7 Interpretation of the relic imprint

The imprint is a freeze-frame of the shape of objectness at the late stage of the Mush, when space had already begun to seep through and light could for the first time escape the monolith: the contrast on the map is the density contrast of the seeping space. The Mush at that point is not a strictly zero state but an offset from zero: `μ_cmb ≈ 0.9999995` — space has already seeped through slightly. The sign of the contrast in this picture: **cold/dark patches ↔ dense objectness (the monolith), hot ones ↔ pores of space that have seeped through**; the smallness of the contrast, `~10⁻⁵`, is qualitatively consistent with the pores of space only just being born near the threshold. For sobriety's sake, let me note: this sign coincides with the standard large-scale Sachs–Wolfe effect for adiabatic perturbations (overdense regions look colder) — a pleasant consistency, but not evidence, since the amplitude and the spectrum are not derived from the model. Finally, an important distinction: `μ_cmb` is the observed "floor" of the early zone; `μ_d` (`z ≈ 4.53`) is its late edge, the detachment threshold.

---

## 12. Structure growth: a check on the DSIC background

Earlier (§8), structure growth was declared outside the model. Here a minimal, honestly bounded step is taken: I check whether the DSIC **background** is compatible with the standard growth of linear perturbations — not deriving the growth equation from the core, but taking it in the GR limit on the DSIC geometry. This is a consistency check, not a derivation.

On a background with rate `H(z)` (§5.4), the matter density contrast `δ` in the sub-horizon GR limit obeys the standard equation

```
d²δ/d(ln a)² + [2 + dln H/dln a] · dδ/d(ln a) − (3/2) Ω_m(a) · δ = 0 ,
```

where `Ω_m(a) = Ω_m0 · a⁻³ / E(a)²`, `E = H/H₀`. Since the core introduces no densities, `Ω_m0` enters here as a growth parameter. From the solution one obtains the growth rate `f = dln δ/dln a` and the observable combination `fσ₈(z) = f(z)·σ₈,₀·δ(z)/δ(0)`, compared with measurements of redshift-space distortions (RSD).

The check is performed on the **"Gold-2017"** compilation (18 weakly correlated `fσ₈` points; for the three WiggleZ points [15] at `z = 0.44, 0.60, 0.73` their covariance matrix is taken into account). The key feature of `fσ₈` data is a degeneracy: they constrain not `Ω_m0` and `σ₈` separately but their combination `S₈ = σ₈·√(Ω_m0/0.3)`. The result is therefore reported through `S₈`:

```
DSIC (GR limit on the DSIC background), μ₀ = 0.76:
  S₈ = 0.761 ± 0.013         (stable along the entire degeneracy valley)
  at fixed Ω_m0 = 0.30:  σ₈ = 0.770,  χ²/dof ≈ 0.73
```

A direct comparison with ΛCDM on the same 18 points gives `S₈ = 0.760`, `χ²/dof ≈ 0.74`; the difference `Δχ² ≈ +0.1` is statistically indistinguishable. The DSIC `fσ₈(z)` curve is shallow, peaking around `z ≈ 0.5`, and on current data it describes the measurements no worse than the standard model. This does not prove that growth in DSIC is obliged to follow the GR-limit equation — the core does not dictate it — but it removes the objection of gross incompatibility: the DSIC background geometry is compatible with the observed growth of structure at reasonable `Ω_m0, σ₈`, and the recovered `S₈ ≈ 0.76` agrees with the standard value.

A caveat on systematics: the Alcock–Paczynski correction for the difference between the DSIC background distances and the fiducial cosmology of each RSD measurement was not applied; at `z < 1` this factor typically amounts to a few percent and is small against the data errors. A rigorous derivation of the growth equation from binding (Appendix B) remains a direction for development.

---

## 13. Conclusion

From four structural elements — the relativity of scale, the conservation of the two-pole norm, a single uniform flow, and the internal observer's modulus-only access — the scale factor `b(μ) = √(1−μ²)/μ²` follows uniquely, and from it the entire background kinematics: the redshift `1+z = b(μ₀)/b(μ_e)`, the rate `|H_τ| = (2−μ²)/(1−μ²)`, the sign-changing deceleration parameter `q(μ)` with the crossover `μ_cross ≈ 0.874`, closed-form distances, and a finite particle horizon. The single parameter `μ₀ ≈ 0.76` consistently passes six independent probes of the late Universe, remaining statistically indistinguishable from ΛCDM on its best data and sharply distinguishable from it at high `z`. In this picture, the acceleration of the expansion is not the action of a hidden component but the shadow of the drift of two scales, read by an observer who is itself made of one of them.

The second tier adds to this one entity — the binding of objects — and one measured constant. The detachment postulate P6 fixes the form of the correction to the light path by principle (light has nowhere to propagate until the average medium has parted; local clumps fall into themselves earlier), the CMB calibration fixes the threshold `μ_d ≈ 0.9806`; the late Universe is thereby untouched identically, the distance to the last-scattering surface is closed exactly, and distance duality and achromaticity hold by construction. The core parameter is measured: `μ₀ = 0.7600 ± 0.0091`. What remains open are the origin of the acoustic scale `~147 Mpc`, the compressed cosmic time at `z = 5–12` (at `z = 7.5` — 0.42 Gyr versus 0.71 in ΛCDM, §11.6), and structure growth — it is there that the model will meet its decisive tests.

---

## References

1. A. G. Riess et al., *Observational Evidence from Supernovae for an Accelerating Universe and a Cosmological Constant*, AJ **116**, 1009 (1998).

2. S. Perlmutter et al., *Measurements of Ω and Λ from 42 High-Redshift Supernovae*, ApJ **517**, 565 (1999).

3. D. M. Scolnic et al., *The Pantheon+ Analysis: The Full Dataset and Light-Curve Release*, ApJ **938**, 113 (2022); arXiv:2112.03863.

4. D. Brout et al., *The Pantheon+ Analysis: Cosmological Constraints*, ApJ **938**, 110 (2022); arXiv:2202.04077.

5. A. G. Riess et al., *A Comprehensive Measurement of the Local Value of the Hubble Constant* (SH0ES), ApJL **934**, L7 (2022); arXiv:2112.04510.

6. DESI Collaboration (M. Abdul-Karim et al.), *DESI DR2 Results II: Measurements of Baryon Acoustic Oscillations and Cosmological Constraints*, Phys. Rev. D **112**, 083515 (2025); arXiv:2503.14738.

7. Planck Collaboration (N. Aghanim et al.), *Planck 2018 results. VI. Cosmological parameters*, A&A **641**, A6 (2020); erratum A&A **652**, C4 (2021); arXiv:1807.06209.

8. E. Di Valentino et al., *In the Realm of the Hubble Tension — a Review of Solutions*, Class. Quantum Grav. **38**, 153001 (2021); arXiv:2103.01183.

9. M. Moresco et al., *Setting the Stage for Cosmic Chronometers* (cosmic chronometer covariance), ApJ **898**, 82 (2020); arXiv:2003.07362. (With the compilation of Moresco et al. 2012; Moresco 2015; Moresco et al. 2016. Original sources of the remaining `H(z)` points: Simon et al. 2005; Stern et al. 2010; Zhang et al. 2014; Ratsimbazafy et al. 2017; Borghi et al. 2022.)

10. O. Farooq & B. Ratra, *Hubble Parameter Measurement Constraints on the Cosmological Deceleration–Acceleration Transition Redshift*, ApJ **766**, L7 (2013); arXiv:1301.5243.

11. V. Sahni, A. Shafieloo & A. A. Starobinsky, *Two New Diagnostics of Dark Energy* (the Om diagnostic), Phys. Rev. D **78**, 103502 (2008); arXiv:0807.3548.

12. V. Sahni, T. D. Saini, A. A. Starobinsky & U. Alam, *Statefinder — a New Geometrical Diagnostic of Dark Energy*, JETP Lett. **77**, 201 (2003); arXiv:astro-ph/0201498.

13. S. Nesseris, G. Pantazis & L. Perivolaropoulos, *Tension and Constraints on Modified Gravity Parametrizations of G_eff(z) from Growth Rate and Planck Data*, Phys. Rev. D **96**, 023542 (2017); arXiv:1703.10538. (The "Gold-2017" compilation, 18 `fσ₈` points.)

14. B. Sagredo, S. Nesseris & D. Sapone, *Internal Robustness of Growth Rate Data*, Phys. Rev. D **98**, 083543 (2018); arXiv:1806.10822. (Table I — the `z, fσ₈, σ` values used in §12.)

15. C. Blake et al., *The WiggleZ Dark Energy Survey: Joint Measurements of the Expansion and Growth History at z < 1*, MNRAS **425**, 405 (2012); arXiv:1204.3674. (Covariance of the three `fσ₈` points at `z = 0.44, 0.60, 0.73`.)

---

## Reproducibility

All empirical claims of the paper are reproducible with a single self-contained script, `dsic_test.py`, from the repository https://github.com/az-zase/DSIC

The provenance of every hard-coded figure is documented in the file itself; all reference data pass self-checks.

Script blocks:

Part I: `[0]` line-by-line check of the hard-coded DESI DR2 against the official Table IV; `[0b]` machine verification of the core's internal identities (monotonicity of `b(μ)`, accuracy of `D_H = dχ/dz` to `~10⁻⁸`, agreement of the closed-form `χ(μ)` with the numerical `∫c dz/H`, the zero of `q` at `μ_cross = √(3−√5)`, the age formula); `[1]` the anchored SN fit to Pantheon+SH0ES (auto-download of the catalog and covariances, STAT+SYS, explicit `M`) with 1σ parameter uncertainties from the `χ²` Hessian; `[2]` the STAT-ONLY control; `[3]` the joint SN+BAO fit with the full 12×12 covariance and the recovery of the `r_d` ruler; `[4]` the `Om(z)` diagnostic with an `H₀` scan in two chronometer-covariance modes (diagonal and the full Moresco et al. 2020); `[5]` high-`z` predictions; `[6]` cosmic time `t(z)`, DSIC vs ΛCDM.

Part II: `[7]` self-check of the Planck 2018 reference data (`d_LSS = r*/θ*`); `[8]` the core's shortfall and the fixing of `μ_d` by the angle `θ*` (the borrowed `r*` is flagged explicitly); `[9]` verification of the identically zero correction below the threshold and the "amplitude ↔ threshold" degeneracy table; `[10]` the trace of `μ_d` in `z`, the corrected particle horizon, and the consistency flag of the P6 formulation (JWST galaxies above the threshold).

Part III: `[11]` (structure growth): integration of the standard linear growth equation `δ'' + (2 + dlnH/dlna)δ' − (3/2)Ω_m(a)δ = 0` on the DSIC background (GR limit), computation of `fσ₈(z)` and the check against the "Gold-2017" compilation (18 points, the WiggleZ covariance included); the result is reported through `S₈` owing to the `Ω_m0`–`σ₈` degeneracy, plus a comparison with ΛCDM on the same points.

---

## Appendix A. Summary of notation

| layer | symbol | meaning | sign / range |
|---|---|---|---|
| 0 | `φ = ωt` | flow phase (time) | monotonically increasing |
| 0 | `ω`, `φ̇` | flow speed | `> 0` |
| 1 | `R` | scale budget (norm) | `> 0` |
| 1 | `S = R·sin(φ/2)` | scale of space | signed (odd) |
| 1 | `O = R·cos(φ/2)` | scale of objectness | signed (even) |
| 1 | `σ_S, σ_O` | sign phase of the poles | `±1` |
| 1 | `O_full = O·σ_S` | objectness with inherited sign | signed |
| 1 | `W = −R²/2` | Wronskian (quadrature) | const |
| 2 | `P = \|S\| = R√(1−μ²)` | observable scale of space | `≥ 0` |
| 2 | `Q = \|O\| = R·μ` | observable scale of objectness | `≥ 0` |
| 2 | `μ = Q/R = \|cos(φ/2)\|` | observable phase (fraction of objectness) | `∈ [0,1]`, period `2π` |
| 3 | `L(μ) ∝ μ` | the observer's ruler | — |
| 3 | `τ`, `dτ = dt/μ` | the observer's proper time | — |
| 3 | `c` | speed of light | const |
| 3 | `b = √(1−μ²)/μ²` | observable scale factor | — |
| 3 | `z` | redshift, `1+z = b(μ₀)/b(μ_e)` | `≥ 0` |
| 3 | `s = ±1` | branch direction (−1 expansion, +1 contraction) | — |
| 3 | `k = (2/π)ω` | speed of `μ` along a branch | `> 0` |
| 3 | `H_τ`, `q` | expansion rate, deceleration parameter | — |
| 3 | `w_eff = (2q−1)/3` | effective equation of state (GR translation, §5.6) | — |
| 3 | `χ`, `d_L`, `DM` | comoving dist., luminosity dist., distance modulus | — |
| 3 | `μ₀, μ_e` | phase today / at emission | `∈ [0,1]`; physical reception branch `μ_e ∈ (μ₀, 1)` |
| 3 | `z_t` | redshift of the decel.→accel. transition | `≈ 0.77` |

**Notation of Part II (the second tier):**

| symbol | meaning | status |
|---|---|---|
| binding | the second entity of the second tier (absent from the core) | postulated (P6) |
| `μ_d = 0.980640` | detachment threshold (`z_detach ≈ 4.53`) | measured constant (calibrated to `d_LSS`) |
| `A ≡ 1` | amplitude of the `1/P` correction | fixed by postulate P6 |
| `D_M_corr` | comoving distance with the `arcsin` correction added | — |
| `d_LSS, z*, r*, θ*` | Planck 2018 reference quantities | external data |
| `r_ij`, `M_i`, `G` | scale relation of a pair, masses, Newton's constant | Appendix B module, not from the core |
| `g_ij`, `w_ij` | pair binding, space-distribution weight | Appendix B module |
| `κ` | space redistribution rate | free parameter of Appendix B |

---

## Appendix B. A development program: gravity and the distribution of space

> **Status of this section.** This is not a result but a sketch of a direction. The module below is phenomenological: it uses the Newtonian quantities `G, M_i`, which the `S ↔ O` flow does not contain; the parameter `κ` is free; no confrontation with data has been performed. It is placed in an appendix deliberately — to separate the tested core (Part I) and the calibrated second tier (Part II) from hypotheses not yet brought to quantitative testing. The module's value is architectural: it shows that binding (Part II) can also have a late, gravitational face, closing the early (the detachment threshold) and the late (the distribution of space) onto a single entity.

A late manifestation of the same binding. When objects fall into themselves (§3.4), gravitationally bound groups fall in concert — a **common center of infall** shows through in them: the more strongly objects are bound, the more distinct the common center and the less space is released between them; between weakly bound ones — more. The loop is self-reinforcing: dispersal weakens the bond → the region receives more space → dispersal accelerates. Dense structures stay compact, sparse ones receive more of the seeping space — qualitatively this agrees with what is observed: space accrues not where it is dense but where it is empty. Space is redistributed, the global balance is conserved, the anisotropy is local.

Module notation: `r_ij` — the scale relation of a pair; `M_i` — an object's mass; `G = 1`; `κ` — a free parameter (the redistribution rate).

```
binding:              g_ij = G·M_i·M_j / r_ij²
potential gradient:   |∇Φ_ij| = | Σ_k G·M_k·(mid_ij − pos_k)/|mid_ij − pos_k|³ |
space weight:         w_ij = (1/|∇Φ_ij|) / Z ,     Z = Σ_{i≠j} 1/|∇Φ_ij|
increment:            Δr_ij = κ·|ΔP|·w_ij / (M_i + M_j)
evolution:            r_ij(φ+Δφ) = r_ij(φ) + Δr_ij         (growth only)
conservation:         Σ_{i≠j} Δr_ij·(M_i + M_j)/κ = |ΔP(φ)|
```

Properties of the weights: `w_ij ≥ 0`, `Σ w_ij = 1`; a pair in a void (a flat field) receives a large weight, a pair in a cluster (a steep field) — a small one.

The module's status is phenomenological: it uses the Newtonian quantities `G, M_i`, which the `S ↔ O` flow does not contain; `κ` is free; no confrontation with data has been performed. But architecturally it now stands on the same footing as §11: binding is one entity with two faces. The early face is the detachment threshold (the condition under which the monolith's binding breaks and the observer's rulers acquire meaning); the late face is the distribution of the seeping space in the places of least gravitational binding.
The minimal program: to derive `μ_d` as a condition on binding from the same root as the weights `w_ij` — then both constructions would become consequences of one quantity rather than two add-ons.

Thank you for reading my thrilling paper to the end.

Adventurer — A. Zaseev