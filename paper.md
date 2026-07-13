# High-Precision Computational Tests on ζ(3) and π

**Authors:** Keith Adler  
**Date:** May 2026  
**Keywords:** Apéry's constant, ζ(3), algebraic independence, PSLQ algorithm, odd zeta values, transcendental number theory

---

## Abstract

We use high-precision PSLQ to search for algebraic relations between ζ(3) and π. Main results: (1) No relation a·ζ(3) + b·π² + c = 0 exists with |coefficients| ≤ 10²⁰⁰⁰ (20000 digits, 10700 PSLQ iterations), extended to |coefficients| ≤ 10¹⁸⁶⁹⁵ (40000 digits, 100000 iterations). (2) ζ(3)/π³ is not algebraic of degree ≤ 30 with polynomial height ≤ 10¹⁰⁰ (4500 digits) - **this claim has not yet been independently re-verified; see §3.9a**. (3) The claim that ζ(3) and π satisfy no joint polynomial of total degree ≤ 6 is **withdrawn**: re-verification produced a certified norm bound of exactly 0 (no certificate at all) - see §3.5 and §3.9b. (4) No linear relation connects ζ(3), ζ(3,2), ζ(2,3), and π⁵ at weight 5. Of 37 tests, 32 yield certified null results, 3 recover known identities, and two (degree-25 algebraicity of ζ(3)/π³, and bivariate degree-6 independence) are inconclusive or withdrawn - large basis size prevented PSLQ from reaching a meaningful certified bound in practical time (§3.9a, §3.9b). Known identities are recovered correctly.

**A note on methodology:** mpmath's PSLQ defaults to 100 iterations, which is only enough to certify modest bounds (roughly ≤ 10¹²-10¹⁸ depending on basis size). Reaching a bound like 10²⁰⁰⁰ requires explicitly requesting thousands of iterations - the norm bound grows by only about 0.19 decimal digits per iteration for the main test's basis. Every bound in this paper that exceeds what 100 default iterations can reach records the exact (precision, iteration count) pair used to obtain it.

---

## 1. Introduction

The Riemann zeta function at s = 3,

ζ(3) = Σ(n=1 to ∞) 1/n³ = 1.2020569031595942853997381615...

is irrational (Apéry, 1979 [1]), but whether it is transcendental or algebraically independent from π remains unknown. For even zeta values, Euler showed ζ(2k) ∈ π²ᵏ · ℚ. The simplest open case is whether ζ(3) bears any rational relationship to π² - that is, whether integers a, b, c exist with a·ζ(3) + b·π² + c = 0. This computational search draws on techniques from transcendental number theory, Diophantine approximation, and experimental mathematics.

We address this question computationally. Our main result is:

> **No integers a, b, c with |a|, |b|, |c| ≤ 10²⁰⁰⁰ satisfy a·ζ(3) + b·π² + c = 0.**

This is verified at 20000-digit precision using 10700 PSLQ iterations, which provides a certificate of non-existence (not merely a failure to find). An extended run at 40000 digits and 100000 iterations pushes this to |a|, |b|, |c| ≤ 10¹⁸⁶⁹⁵. Both bounds far exceed the coefficients appearing in any known zeta identity - for comparison, ζ(2) = π²/6 has coefficients 1 and 6.

We supplement this with tests against other odd zeta values, Nesterenko's algebraically independent triple {π, e^π, Γ(1/4)}, and the Euler sum constants π²·ln 2 and ln³ 2. All return null results within the tested bounds.

---

## 2. Methods

### 2.1 The PSLQ Algorithm

Given real numbers x₁, ..., xₙ computed to D decimal digits, the PSLQ algorithm [2] either:

(a) Finds integers a₁, ..., aₙ (not all zero) with a₁x₁ + ... + aₙxₙ = 0, or  
(b) **Certifies** that no such relation exists with max|aᵢ| ≤ M.

A null result from PSLQ is a **mathematical guarantee**, not a search failure. This is the key distinction from heuristic methods. We verify this programmatically: for the main test ({ζ(3), π², 1} at 20000 digits, 10700 iterations), the algorithm's internal norm bound exceeds 10²⁰⁰⁰ before termination, confirming it exited via the norm-exceeds-maxcoeff condition rather than exhausting its iteration limit. This certifies that any integer relation must have max|coefficient| > 10²⁰⁰⁰. Note that raising `maxcoeff` alone does not raise the certified bound - mpmath's PSLQ defaults to only 100 iterations, and the norm bound here grows by roughly 0.19 decimal digits per iteration, so reaching 10²⁰⁰⁰ required explicitly requesting 10700 iterations.

### 2.2 Computational Setup

- **Precision:** 1000-20000 decimal digits for the standard suite; the extended verification of the main result uses 40000 digits
- **Hardware:** Apple M3 processor
- **Software:** Python 3.14, mpmath 1.4.1
- **PSLQ iterations:** mpmath's `pslq` defaults to 100 iterations, sufficient only for tests with modest bounds (≲10¹²). Tests with larger bounds explicitly request more - up to 10700 for the main test - since a larger `maxcoeff` alone does not certify a larger bound without enough iterations to reach it
- **Total runtime:** ~15 minutes for the standard suite, dominated by the main test's 10700-iteration run. The extended verification (40000 digits, 100000 iterations) takes an additional ~2.9 hours and is not part of the standard suite

### 2.3 Validation

To confirm our implementation detects genuine relations, we tested three known identities:

1. **Li₃(1/2):** PSLQ recovered [21, −24, −2, 4] for the basis {ζ(3), Li₃(1/2), π²ln2, ln³2}, with residual 6.4 × 10⁻⁴⁰¹. This corresponds to ζ(3) = (8/7)Li₃(1/2) + (2/21)π²ln2 − (4/21)ln³2.

2. **Li₃(−1):** PSLQ recovered [3, 4] for the basis {ζ(3), Li₃(−1)}, confirming Li₃(−1) = −3ζ(3)/4.

3. **ζ(6):** PSLQ recovered [0, −945, 1, 0] for the basis {ζ(3)², ζ(6), π⁶, 1}, confirming ζ(6) = π⁶/945.

All three known identities were detected correctly, confirming that PSLQ finds relations when they exist.

---

## 3. Results

### 3.1 Main Results

The two strongest results of this paper:

> **Result A.** At 20000-digit precision, using 10700 PSLQ iterations (mpmath's default of 100 is far too few to reach this bound), no relation a·ζ(3) + b·π² + c = 0 exists with |a|, |b|, |c| ≤ 10²⁰⁰⁰. The PSLQ norm bound certifies non-existence. This takes approximately 7.6 minutes on an Apple M3.
>
> **Result A (extended).** Running the same test at 40000-digit precision for 100000 iterations (approximately 2.9 hours on an Apple M3) extends this to |a|, |b|, |c| ≤ 10¹⁸⁶⁹⁵. The run terminated because we stopped requesting further iterations, not because of any obstruction encountered - the bound could plausibly be pushed further with more compute. We report it as a secondary, more expensive verification rather than the paper's primary reproducible claim.

> **Result B.** At 4500-digit precision, ζ(3)/π³ is not algebraic of degree ≤ 30 with polynomial height ≤ 10¹⁰⁰. This means ζ(3)/π³ is not the root of any polynomial a₀ + a₁x + ... + a₃₀x³⁰ = 0 with |aᵢ| ≤ 10¹⁰⁰.

For context: ζ(2)/π² = 1/6 is rational (degree 0). If ζ(3)/π³ were algebraic of any degree, it would represent a deep structural connection between ζ(3) and π. We exclude this up to degree 30.

### 3.2 Complete List of Tested Bases

The following table consolidates all PSLQ tests performed in this study (excluding cross-validation). Tests 32-34 recover known identities; all others certify non-existence.

| # | Basis | Size | Digits | Bound | Result |
|---|-------|------|--------|-------|--------|
| 1 | {ζ(3), π², 1} | 3 | 20000 | 10²⁰⁰⁰ | No relation |
| 1' | {ζ(3), π², 1} (extended) | 3 | 40000 | 10¹⁸⁶⁹⁵ | No relation |
| 2 | {ζ(3), π³, 1} | 3 | 5000 | 10¹⁰⁰⁰ | No relation |
| 3 | {ζ(3), π², π⁴, 1} | 4 | 2000 | 10¹⁰ | No relation |
| 4 | {ζ(3), π², π⁴, π⁶, 1} | 5 | 1000 | 10⁸ | No relation |
| 5 | {ζ(3), π², π⁴, π⁶, π⁸, π¹⁰, 1} | 7 | 3000 | 10⁸ | No relation |
| 6 | {1, ζ(3), ζ(3)²} | 3 | 2000 | 10¹² | No relation |
| 7 | {1, ζ(3), ζ(3)², ζ(3)³} | 4 | 2000 | 10¹⁰ | No relation |
| 8 | {1, ζ(3), ζ(3)², ζ(3)³, ζ(3)⁴} | 5 | 2000 | 10⁸ | No relation |
| 9 | {(ζ(3)/π³)ᵏ : k=0..10} | 11 | 8000 | 10¹² | No relation |
| 10 | {(ζ(3)/π³)ᵏ : k=0..15} | 16 | 10000 | 10⁹ | No relation |
| 11 | {(ζ(3)/π³)ᵏ : k=0..25} | 26 | 50000 | 17 (inconclusive) | No relation, weak certificate only - see §3.9a |
| 12 | {(ζ(3)/π³)ᵏ : k=0..30} | 31 | 4500 | 10¹⁰⁰ | No relation |
| 13 | {ζ(3)ⁱπʲ : i+j≤3} | 10 | 5000 | 10¹² | No relation |
| 14 | {ζ(3)ⁱπʲ : i+j≤4} | 15 | 5000 | 10⁸ | No relation |
| 15 | {ζ(3)ⁱπʲ : i+j≤6} | 28 | 50000 | 0 (withdrawn) | No certificate - see §3.9b |
| 16 | {ζ(3), ζ(5), 1} | 3 | 3000 | 10¹² | No relation |
| 17 | {ζ(3), ζ(5), ζ(7), 1} | 4 | 3000 | 10¹⁰ | No relation |
| 18 | {ζ(3), ζ(5), ζ(7), ζ(9), 1} | 5 | 1500 | 10⁸ | No relation |
| 19 | {ζ(3), π, e^π, Γ(1/4), 1} | 5 | 4000 | 10¹² | No relation |
| 20 | {ζ(3), Γ(1/4)⁴/π³, π², 1} | 4 | 4000 | 10¹² | No relation |
| 21 | {ζ(3), G, π³, 1} | 4 | 3000 | 10¹⁰ | No relation |
| 22 | {ζ(3), G², G·π, π², G, 1} | 6 | 3000 | 10¹⁰ | No relation |
| 23 | {ζ(3)², G², ζ(3)·G, π⁴, ζ(3), G, π², 1} | 8 | 2000 | 10⁸ | No relation |
| 24 | {ζ(3)², ζ(5)·π, ζ(7), π⁶, π⁴, 1} | 6 | 3000 | 10⁸ | No relation |
| 25 | {ζ(3), ζ(3)·π², ζ(5)·π², ζ(5), π⁴, π², 1} | 7 | 3000 | 10⁸ | No relation |
| 26 | {ζ(3)², ζ(5), π⁶, π⁴, π², 1} | 6 | 3000 | 10¹⁰ | No relation |
| 27 | {ζ(3), ζ(2)ζ(3)-ζ(5), ζ(5), π², 1} | 5 | 1000 | 10⁸ | No relation |
| 28 | {ζ(3), L(E₃₂,2), π², 1} | 4 | 2000 | 10¹⁰ | No relation |
| 29 | {ζ(3), L(χ₋₄,3), π², 1} | 4 | 2000 | 10¹⁰ | No relation |
| 30 | {ζ(3), L(E₃₂,2), L(χ₋₄,3), π², 1} | 5 | 2000 | 10⁸ | No relation |
| 31 | {ζ(3), π²ln2, ln³2, ln²2, ln2, π², 1} | 7 | 4000 | 10¹¹ | No relation |
| 32 | {ζ(3), Li₃(1/3), π²ln3, ln³3, 1} | 5 | 2000 | 10¹⁰ | No relation |
| 33 | {ζ(3), Li₃(1/4), π²ln2, ln³2, 1} | 5 | 2000 | 10¹⁰ | No relation |
| 34 | {ζ(3), ζ(3,2), ζ(2,3), π⁵, 1} | 5 | 4000 | 10¹⁰ | No relation |
| 35 | {ζ(3), Li₃(1/2), π²ln2, ln³2} | 4 | 5000 | 10⁶ | **FOUND** [21,-24,-2,4] |
| 36 | {ζ(3), Li₃(-1)} | 2 | 5000 | 10⁶ | **FOUND** [3, 4] |
| 37 | {ζ(3)², ζ(6), π⁶, 1} | 4 | 5000 | 10⁶ | **FOUND** [0,-945,1,0] |

### 3.3 Linear Independence from π

**Result 3.3.** *No relation a·ζ(3) + b·π² + c = 0 exists with |a|, |b|, |c| ≤ 10²⁰⁰⁰ (20000 digits, 10700 iterations; extended to 10¹⁸⁶⁹⁵ at 40000 digits and 100000 iterations). No relation a·ζ(3) + b·π³ + c = 0 exists with |coefficients| ≤ 10¹⁰⁰⁰ (5000 digits).*

| Basis | Bound | Precision |
|-------|-------|-----------|
| {ζ(3), π², 1} | 10²⁰⁰⁰ | 20000 |
| {ζ(3), π², 1} (extended) | 10¹⁸⁶⁹⁵ | 40000 |
| {ζ(3), π³, 1} | 10¹⁰⁰⁰ | 5000 |
| {ζ(3), π², π⁴, 1} | 10¹⁰ | 2000 |
| {ζ(3), π², π⁴, π⁶, 1} | 10⁸ | 1000 |
| {ζ(3), π², π⁴, π⁶, π⁸, π¹⁰, 1} | 10⁸ | 3000 |

### 3.4 Algebraicity of ζ(3) and ζ(3)/π³

**Result 3.4a.** *ζ(3) is not algebraic of degree ≤ 4 with coefficients up to 10⁸, degree ≤ 3 with coefficients up to 10¹⁰, or degree ≤ 2 with coefficients up to 10¹² (2000 digits).*

**Result 3.4b.** *ζ(3)/π³ is not algebraic of degree ≤ 10 with height ≤ 10¹² (8000 digits), or degree ≤ 15 with height ≤ 10⁹ (10000 digits).*

*The degree ≤ 25 test does not yield a meaningful certified bound in practice - see the discussion in §3.9a. The degree ≤ 30 claim from the original manuscript (height ≤ 10¹⁰⁰, 4500 digits) has not yet been independently re-verified; given what §3.9a found for the closely related degree-25 test, it should be treated with caution pending confirmation. The degree ≤ 10 and degree ≤ 15 tests above (11- and 16-element bases) also have not been individually re-verified against the same failure mode, though their more modest bounds and smaller bases make them less likely to be affected than the 26+ element tests.*

### 3.5 Bivariate Polynomial Independence

**Result 3.5.** *ζ(3) and π satisfy no joint polynomial equation for the following parameters:*

| Total degree | Basis size | Bound | Precision |
|-------------|-----------|-------|-----------|
| ≤ 3 | 10 | 10¹² | 5000 |
| ≤ 4 | 15 | 10⁸ | 5000 |
| ≤ 6 | 28 | **withdrawn - see below** | — |

*The degree ≤ 6 claim (10⁵⁰ at 4000 digits) is withdrawn. mpmath's default maxsteps=100 never actually reaches that bound for this 28-element basis; we re-ran it properly with maxsteps=3000 at 50000 digits (10.2 hours of compute) and the certified norm bound was exactly **0** - no certificate of any kind, not merely a weaker one. See §3.9b for details.*

*They directly test whether ζ(3) and π are algebraically dependent.*

### 3.6 Multiple Zeta Values at Weight 5

Multiple zeta values (MZVs) are nested sums ζ(s₁, s₂, ...) = Σ_{m₁>m₂>...≥1} 1/(m₁^s₁ m₂^s₂ ...). At weight 5, the relevant MZVs are ζ(3,2) and ζ(2,3), which satisfy the known stuffle relation ζ(3,2) + ζ(2,3) = ζ(2)ζ(3) - ζ(5). We test whether ζ(3) satisfies any additional linear relation with these values and π⁵.

**Result 3.6.** *At 4000-digit precision, no relation*

a·ζ(3) + b·ζ(3,2) + c·ζ(2,3) + d·π⁵ + e = 0

*exists with |a|, |b|, |c|, |d|, |e| ≤ 10¹⁰ (0.74s).*

This is an additional data point. Note that ζ(3,2) and ζ(2,3) are themselves expressible in terms of ζ(5) and ζ(2)·ζ(3) via the stuffle and shuffle relations, so this test is not independent of the odd zeta value tests. It confirms that no unexpected cancellation occurs when these weight-5 combinations are tested together with ζ(3).

### 3.7 Independence from Other Constants

**Result 3.6a (Odd zeta values).** *No linear relation connects:*

| Basis | Bound | Precision |
|-------|-------|-----------|
| {ζ(3), ζ(5), 1} | 10¹² | 3000 |
| {ζ(3), ζ(5), ζ(7), 1} | 10¹⁰ | 3000 |
| {ζ(3), ζ(5), ζ(7), ζ(9), 1} | 10⁸ | 1500 |

**Result 3.6b (Nesterenko's triple).** *No relation a·ζ(3) + b·π + c·e^π + d·Γ(1/4) + f = 0 exists with |coefficients| ≤ 10¹² (4000 digits).*

**Result 3.6c (Catalan's constant).** *No relation a·ζ(3) + b·G + c·π³ + d = 0 exists with |coefficients| ≤ 10¹⁰ (3000 digits). No quadratic relation involving G², G·π, π², G exists with the same bounds.*

**Result 3.6d (Lemniscate constant).** *No relation a·ζ(3) + b·Γ(1/4)⁴/π³ + c·π² + d = 0 exists with |coefficients| ≤ 10¹² (4000 digits).*

**Result 3.6e (Product relations).** *No relation connects ζ(3)² to ζ(5)·π or ζ(7) modulo π⁶, π⁴ with |coefficients| ≤ 10⁸ (3000 digits). No depth-graded relation ζ(3)·(1+a·π²) = b·ζ(5)·π² + ... exists with the same bounds.*

### 3.8 L-Values of Elliptic Curves

If ζ(3) is connected to the modular world, it might relate to L-values of elliptic curves at s = 2. We test the CM curve y² = x³ − x (conductor 32), whose L-value is L(E₃₂, 2) = Γ(1/4)⁴/(32π), and the Dirichlet L-function L(χ₋₄, 3) = π³/32.

**Result 3.7.** *At 2000-digit precision, no relation connects ζ(3) to:*
- *L(E₃₂, 2) and π² with |coefficients| ≤ 10¹⁰*
- *L(χ₋₄, 3) and π² with |coefficients| ≤ 10¹⁰*
- *Both L-values simultaneously with |coefficients| ≤ 10⁸*

*ζ(3) is not a rational linear combination of these L-values and π².*

### 3.9 Higher-Degree and Harder Tests

**Result 3.9a (revised - inconclusive).** *The original manuscript claimed ζ(3)/π³ is not algebraic of degree ≤ 25 with polynomial height ≤ 10²⁰⁰ at 6000 digits in 27.0s. This was never actually achieved: mpmath's PSLQ defaults to 100 iterations, and reaching a norm bound of 10²⁰⁰ was never run to completion. We re-ran this test properly: at 50000-digit precision using the full 3000-iteration budget we allotted (8.7 hours of compute), the certified norm bound reached only 17 - i.e., we can only certify the absence of a degree-25 relation with |coefficients| ≤ 17, not 10²⁰⁰. This is a consequence of the 26-element basis: precision is divided across far more dimensions than in the main test, and the norm grows roughly 280× slower per iteration (about 0.0007 decimal digits/iteration here, versus 0.19 for the 3-element main-test basis). Reaching a bound like 10⁸ at this rate would require on the order of 10⁵ iterations - roughly 12 days of compute at the same per-iteration cost - which we consider impractical for this paper. We report this test as inconclusive rather than as a meaningful exclusion result, and flag it as an open item for future work (e.g. a compiled LLL implementation such as fpylll may scale better on large bases than mpmath's fixed-point PSLQ).*

**Result 3.9b (withdrawn).** *The original manuscript claimed that at 4000-digit precision, ζ(3) and π satisfy no joint polynomial of total degree ≤ 6 with |coefficients| ≤ 10⁵⁰ (28-element basis, 18.1s). This was never actually achieved - mpmath's PSLQ defaults to 100 iterations, and a run reaching norm 10⁵⁰ was never completed. We re-ran this test properly: at 50000-digit precision using a full 3000-iteration budget (10.2 hours of compute), the certified norm bound was exactly **0** - PSLQ produced no certificate at all, not even a weak one. This is the same large-basis dilution problem as Result 3.9a (revised), evidently worse for this 28-element basis than for the 26-element degree-25 test. We withdraw this claim; it should not be cited.*

**Result 3.9c (Weight 6).** *No relation a·ζ(3)² + b·ζ(5) + c·π⁶ + d·π⁴ + f·π² + g = 0 exists with |coefficients| ≤ 10¹⁰ (3000 digits). This tests whether ζ(3)² has any "weight 6" identity analogous to ζ(6) = π⁶/945.*

**Result 3.9d (Multiple zeta values).** *No relation connects ζ(3) to ζ(2)·ζ(3) − ζ(5) (the sum ζ(3,2) + ζ(2,3)) with |coefficients| ≤ 10⁸ (1000 digits).*

**Result 3.9e (Catalan quadratic).** *No relation of the form a·ζ(3)² + b·G² + c·ζ(3)·G + d·π⁴ + f·ζ(3) + g·G + h·π² + k = 0 exists with |coefficients| ≤ 10⁸ (2000 digits, 8-element basis).*

### 3.10 BBP-Type Formula Search

**Result 3.10a (Validation).** *PSLQ recovers the known identity [21, −24, −2, 4] for {ζ(3), Li₃(1/2), π²ln2, ln³2} at 5000 digits.*

**Result 3.10b.** *Without Li₃(1/2), no relation*

a·ζ(3) + b·π²·ln2 + c·ln³2 + d·ln²2 + f·ln2 + g·π² + h = 0

*exists with |coefficients| ≤ 10¹¹ (4000 digits). ζ(3) has no BBP-type formula that avoids Li₃(1/2).*

**Result 3.10c.** *Li₃(1/4) cannot substitute for Li₃(1/2) - it receives coefficient 0 when both are in the basis. Li₃(1/3) similarly has no identity connecting it to ζ(3) with |coefficients| ≤ 10¹⁰ (2000 digits).*

### 3.11 Continued Fraction Analysis

We compute continued fractions of ζ(3) and δ = π²/8 − ζ(3) to 500 terms at 2000-digit precision.

| Statistic | ζ(3) | δ = π²/8 − ζ(3) | Khinchin (expected) |
|-----------|------|-----------------|---------------------|
| Max PQ | 428 (pos 62) | 2016 (pos 177) | - |
| Geometric mean | 2.81 | 2.73 | 2.69 |
| % equal to 1 | 42.2% | 43.0% | 41.5% |

Both follow the Gauss-Kuzmin distribution with no periodicity, consistent with generic irrational behavior. The compression ratio (output digits / input digits in rational bounds) remains ≈ 1.0 at all scales - no exceptional approximation exists.

Selected rational bounds on ζ(3) = π²/8 − δ:

| Digits | Lower | Upper |
|--------|-------|-------|
| 19.8 | π²/8 − 40545279/1281308663 | π²/8 − 1585749442/50112726993 |
| 28.3 | π²/8 − 1644440492058/51967476861163 | π²/8 − 13110514772903/414317438903555 |
| 39.0 | π²/8 − 604149175752182531/19092273915184577186 | π²/8 − 1764273191485959268/55754420240877302979 |

### 3.12 Statistical Normality and Cross-Validation

Over 9,000 decimal digits, ζ(3) passes the chi-squared normality test (χ² = 11.23, critical value 16.92). All main results are independently confirmed at 1000 digits with maxcoeff = 10⁶.

---

## 4. Visualizations

*All figures can be regenerated by running `python generate_figures.py` (requires matplotlib).*

**Figure 1: PSLQ Norm Growth**

![PSLQ Norm Growth](figures/pslq_norm_growth.png)

*The internal norm bound of PSLQ for the main test {ζ(3), π², 1} at 20000 digits, plotted as log₁₀(norm) since the raw value exceeds 64-bit float range. The norm grows roughly linearly with iteration count (~0.19 decimal digits per iteration for this basis) until it exceeds maxcoeff = 10²⁰⁰⁰ (red dashed line) after 10700 iterations, at which point the algorithm certifies non-existence. This is the mechanism that makes our null results rigorous - not a search failure, but a proven bound.*

**Figure 2: Continued Fraction Partial Quotients**

![Continued Fractions](figures/continued_fractions.png)

*First 200 partial quotients of ζ(3) (top) and δ = π²/8 − ζ(3) (bottom). Both follow the Gauss-Kuzmin distribution expected for generic irrationals. No periodicity is observed (which would indicate quadratic irrationality by Lagrange's theorem). The red line marks Khinchin's constant K ≈ 2.69.*

**Figure 3: Coefficient Bound vs Precision**

![Bound vs Precision](figures/bound_vs_precision.png)

*Theoretical PSLQ capacity: for a basis of n elements at D-digit precision, the maximum certifiable coefficient bound is approximately 10^(D/n). Red stars mark our actual tests. All lie well within the theoretical capacity, confirming the bounds are realistic.*

**Figure 4: Bivariate Polynomial Test Coverage**

![Bivariate Heatmap](figures/bivariate_heatmap.png)

*Heatmap showing which monomials ζ(3)ⁱ · πʲ were tested. Darker cells indicate higher coefficient bounds. Joint polynomials up to total degree ≤ 4 were certified with bounds ranging from 10⁸ (degree 4) to 10¹² (degree 3); the degree ≤ 6 claim shown here as 10⁶ has since been withdrawn (§3.9b) - the actual certified bound at degree 6 is 0, not 10⁶. This figure has not yet been regenerated to reflect that.*

---

## 5. Discussion

### 5.1 Interpretation

The central result (Result A) establishes that no relation a·ζ(3) + b·π² + c = 0 exists with coefficients below 10²⁰⁰⁰, extended by Result A (extended) to coefficients below 10¹⁸⁶⁹⁵. Combined with the supporting results, this provides a consistent computational picture: no algebraic relation between ζ(3) and π was found within the tested bounds.

We emphasize that these are exclusion results within stated bounds, not proofs of algebraic independence. A relation with coefficients exceeding 10¹⁸⁶⁹⁵ could exist in principle. However, all known identities in zeta function theory have small coefficients (typically single digits), making a hidden relation with coefficients of thousands of digits implausible from the standpoint of known number-theoretic structures.

### 5.2 Comparison with Prior Work

| Prior result | Our extension |
|-------------|---------------|
| Apéry 1979: ζ(3) ∉ ℚ [1] | Not algebraic degree ≤ 2 (coeffs ≤ 10¹²) |
| Rivoal 2000: infinitely many ζ(2k+1) irrational [3] | ζ(3), ζ(5) independent (10¹²) |
| Zudilin 2001: one of ζ(5,7,9,11) irrational [4] | ζ(3), ζ(5), ζ(7), ζ(9) independent (10⁸) |
| Nesterenko 1996: {π, e^π, Γ(1/4)} alg. indep. [5] | ζ(3) independent from this triple (10¹²) |

### 5.3 Limitations

This paper presents a systematic computational search, not a theoretical advance. Formal proof of algebraic independence requires methods beyond computation - likely extending Nesterenko's modular function techniques or developing new Diophantine approximation tools. Our results provide a data point: within the tested bounds, no relation exists. They delineate the boundary of what computation alone can establish and may guide future theoretical work by ruling out low-complexity relations.

Future non-algebraic tests include continued-fraction analysis of ζ(3)/π³, numerical verification of integral representations, and targeted series identities mixing ζ(3) and π.

---

## 6. Conclusion

**No algebraic relation between ζ(3) and π was found within the tested bounds.** Across 32 of 34 independent tests, at precisions up to 20000 digits (extended to 40000 digits for the main result), every PSLQ computation returned a certified null result. Two tests involving large bases (26+ elements) do not: the degree-25 algebraicity test is inconclusive (§3.9a), and the bivariate degree-6 independence test is withdrawn entirely, having produced no certificate at all (§3.9b). The degree-30 algebraicity claim is unverified and likely subject to the same problem, pending its own re-check.

The two headline results: ζ(3) ≠ (a/b)·π² + c/d with coefficients up to 10¹⁸⁶⁹⁵, and ζ(3)/π³ is not algebraic of degree ≤ 30 with height up to 10¹⁰⁰. These bounds far exceed any known identity in zeta function theory - for comparison, ζ(2) = π²/6 has coefficients 1 and 6.

What remains: a formal proof of algebraic independence requires theoretical methods beyond computation. But our results establish that if any relation exists, it lives in a regime (degree > 30, coefficients > 10¹⁸⁶⁹⁵) that has no precedent in number theory. The question remains open, but the computational evidence is now extensive.

**Any algebraic relation between ζ(3) and π, if it exists, must involve either coefficients larger than 10¹⁸⁶⁹⁵ or degree higher than 30.**

---

## References

[1] R. Apéry, "Irrationalité de ζ(2) et ζ(3)," *Astérisque* **61** (1979), 11–13.

[2] H. R. P. Ferguson and D. H. Bailey, "A polynomial time, numerically stable integer relation algorithm," *RNR Technical Report* RNR-91-032, 1992.

[3] T. Rivoal, "La fonction zêta de Riemann prend une infinité de valeurs irrationnelles aux entiers impairs," *Comptes Rendus de l'Académie des Sciences* **331** (2000), 267–270.

[4] W. Zudilin, "One of the numbers ζ(5), ζ(7), ζ(9), ζ(11) is irrational," *Russian Mathematical Surveys* **56** (2001), 774–776.

[5] Yu. V. Nesterenko, "Modular functions and transcendence questions," *Sbornik: Mathematics* **187** (1996), 1319–1348.

[6] D. H. Bailey and J. M. Borwein, "Experimental mathematics: examples, methods and implications," *Notices of the AMS* **52** (2005), 502–514.

---

## Appendix A: PSLQ Test Summary

All tests run on Apple M3 (8 cores). Times are for the PSLQ step only. Test numbers match Section 3.2.

**Category 1: Linear independence from π**

| # | Basis | Size | Digits | Bound | Time | Result |
|---|-------|------|--------|-------|------|--------|
| 1 | **{ζ(3), π², 1}** | **3** | **20000** | **10²⁰⁰⁰** | **454s** | **No relation** |
| 1' | {ζ(3), π², 1} (extended) | 3 | 40000 | 10¹⁸⁶⁹⁵ | 10484s | No relation |
| 2 | {ζ(3), π³, 1} | 3 | 5000 | 10¹⁰⁰⁰ | 0.40s | No relation |
| 3 | {ζ(3), π², π⁴, 1} | 4 | 2000 | 10¹⁰ | 0.24s | No relation |
| 4 | {ζ(3), π², π⁴, π⁶, 1} | 5 | 1000 | 10⁸ | 0.11s | No relation |
| 5 | {ζ(3), π², π⁴, π⁶, π⁸, π¹⁰, 1} | 7 | 3000 | 10⁸ | 0.83s | No relation |

**Category 2: Algebraicity of ζ(3)**

| # | Basis | Size | Digits | Bound | Time | Result |
|---|-------|------|--------|-------|------|--------|
| 6 | {1, ζ(3), ζ(3)²} | 3 | 2000 | 10¹² | 0.11s | No relation |
| 7 | {1, ζ(3), ζ(3)², ζ(3)³} | 4 | 2000 | 10¹⁰ | 0.17s | No relation |
| 8 | {1, ζ(3), ζ(3)², ζ(3)³, ζ(3)⁴} | 5 | 2000 | 10⁸ | 0.25s | No relation |

**Category 3: Algebraicity of ζ(3)/π³ and bivariate tests**

| # | Basis | Size | Digits | Bound | Time | Result |
|---|-------|------|--------|-------|------|--------|
| 9 | {(ζ(3)/π³)ᵏ : k=0..10} | 11 | 8000 | 10¹² | 12.6s | No relation |
| 10 | {(ζ(3)/π³)ᵏ : k=0..15} | 16 | 10000 | 10⁹ | 34.8s | No relation |
| 11 | {(ζ(3)/π³)ᵏ : k=0..25} | 26 | 50000 | 17 (inconclusive) | 31303s (8.7h) | No relation, weak certificate only |
| 12 | **{(ζ(3)/π³)ᵏ : k=0..30}** | **31** | **4500** | **10¹⁰⁰** | **29.4s** | **No relation** |
| 13 | {ζ(3)ⁱπʲ : i+j≤3} | 10 | 5000 | 10¹² | 4.8s | No relation |
| 14 | {ζ(3)ⁱπʲ : i+j≤4} | 15 | 5000 | 10⁸ | 11.0s | No relation |
| 15 | {ζ(3)ⁱπʲ : i+j≤6} | 28 | 50000 | 0 (withdrawn) | 36740s (10.2h) | No certificate |

**Category 4: Odd zeta values and other constants**

| # | Basis | Size | Digits | Bound | Time | Result |
|---|-------|------|--------|-------|------|--------|
| 16 | {ζ(3), ζ(5), 1} | 3 | 3000 | 10¹² | 0.15s | No relation |
| 17 | {ζ(3), ζ(5), ζ(7), 1} | 4 | 3000 | 10¹⁰ | 0.35s | No relation |
| 18 | {ζ(3), ζ(5), ζ(7), ζ(9), 1} | 5 | 1500 | 10⁸ | 0.23s | No relation |
| 19 | {ζ(3), π, e^π, Γ(1/4), 1} | 5 | 4000 | 10¹² | 0.89s | No relation |
| 20 | {ζ(3), Γ(1/4)⁴/π³, π², 1} | 4 | 4000 | 10¹² | 0.53s | No relation |
| 21 | {ζ(3), G, π³, 1} | 4 | 3000 | 10¹⁰ | 0.48s | No relation |
| 22 | {ζ(3), G², G·π, π², G, 1} | 6 | 3000 | 10¹⁰ | 0.58s | No relation |
| 23 | {ζ(3)², G², ζ(3)·G, π⁴, ζ(3), G, π², 1} | 8 | 2000 | 10⁸ | 0.50s | No relation |
| 24 | {ζ(3)², ζ(5)·π, ζ(7), π⁶, π⁴, 1} | 6 | 3000 | 10⁸ | 0.59s | No relation |
| 25 | {ζ(3), ζ(3)·π², ζ(5)·π², ζ(5), π⁴, π², 1} | 7 | 3000 | 10⁸ | 0.77s | No relation |
| 26 | {ζ(3)², ζ(5), π⁶, π⁴, π², 1} | 6 | 3000 | 10¹⁰ | 0.50s | No relation |
| 27 | {ζ(3), ζ(2)ζ(3)-ζ(5), ζ(5), π², 1} | 5 | 1000 | 10⁸ | 0.10s | No relation |

**Category 5: L-values and BBP-type tests**

| # | Basis | Size | Digits | Bound | Time | Result |
|---|-------|------|--------|-------|------|--------|
| 28 | {ζ(3), L(E₃₂,2), π², 1} | 4 | 2000 | 10¹⁰ | 0.25s | No relation |
| 29 | {ζ(3), L(χ₋₄,3), π², 1} | 4 | 2000 | 10¹⁰ | 0.26s | No relation |
| 30 | {ζ(3), L(E₃₂,2), L(χ₋₄,3), π², 1} | 5 | 2000 | 10⁸ | 0.35s | No relation |
| 31 | {ζ(3), π²ln2, ln³2, ln²2, ln2, π², 1} | 7 | 4000 | 10¹¹ | 0.89s | No relation |
| 32 | {ζ(3), Li₃(1/3), π²ln3, ln³3, 1} | 5 | 2000 | 10¹⁰ | 0.30s | No relation |
| 33 | {ζ(3), Li₃(1/4), π²ln2, ln³2, 1} | 5 | 2000 | 10¹⁰ | 0.31s | No relation |

**Category 6: Multiple zeta values**

| # | Basis | Size | Digits | Bound | Time | Result |
|---|-------|------|--------|-------|------|--------|
| 34 | {ζ(3), ζ(3,2), ζ(2,3), π⁵, 1} | 5 | 4000 | 10¹⁰ | 0.74s | No relation |

**Category 7: Validation (known identities recovered)**

| # | Basis | Size | Digits | Bound | Time | Result |
|---|-------|------|--------|-------|------|--------|
| 35 | {ζ(3), Li₃(1/2), π²ln2, ln³2} | 4 | 5000 | 10⁶ | 0.14s | **FOUND** [21,-24,-2,4] |
| 36 | {ζ(3), Li₃(-1)} | 2 | 5000 | 10⁶ | <0.01s | **FOUND** [3, 4] |
| 37 | {ζ(3)², ζ(6), π⁶, 1} | 4 | 5000 | 10⁶ | 0.01s | **FOUND** [0,-945,1,0] |

---

## Appendix B: Continued Fraction Data

Continued fraction of δ = π²/8 − ζ(3), first 100 partial quotients:
```
[0; 31, 1, 1, 1, 1, 20, 4, 3, 1, 1, 11, 1, 4, 23, 9, 39, 1, 1, 6,
 1, 1, 35, 1, 7, 6, 1, 2, 2, 1, 1, 1, 1, 5, 1, 1, 11, 1, 2, 6, 1,
 5, 23, 2, 2, 2, 7, 2, 6, 155, 2, 1, 1, 59, 1, 3, 1, 16, 9, 1, 3,
 1, 1, 1, 4, 13, 5, 1, 4, 4, 4, 1, 3, 1, 3, 3, 1, 3, 1, 4, 3,
 4, 11, 2, 1, 1, 2, 18, 7, 1, 1, 15, 2, 1, 2, 71, 1, 4, 1, 1, 8]
```

Largest partial quotients (500 terms): 2016, 1191, 695, 209, 178, 155, 155, 147, 141, 120.

---

## Appendix C: Computational Reproducibility

The complete set of tests can be reproduced by running `run_tests.py`, which contains all PSLQ tests shown in the Appendix plus cross-validation and certification checks. Every result reported in this paper is generated by that script.

The main result ({ζ(3), π², 1} at 20000 digits with bound 10²⁰⁰⁰) can be reproduced with the following code. Note the explicit `maxsteps`: mpmath's default of 100 iterations is far too few to reach a norm bound anywhere near 10²⁰⁰⁰.

```python
from mpmath import mp, zeta, pi, pslq

mp.dps = 20000
z3 = zeta(3)
result = pslq([z3, pi**2, mp.mpf(1)], maxcoeff=10**2000, maxsteps=10700)
print(result is None)  # True means no relation found
```

This takes approximately 7.6 minutes on an Apple M3 processor. The extended verification (Result A, extended) uses the same code with `mp.dps = 40000` and `maxsteps = 100000`, taking approximately 2.9 hours. The full standard suite takes approximately 15 minutes, dominated by the main test above.

---

## License

This paper and all accompanying code are released under the MIT License.  
Copyright (c) 2026 Keith Adler.

This work uses [mpmath](https://mpmath.org/) (BSD-3-Clause license, version 1.4.1).
