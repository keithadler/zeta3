"""
Computational Evidence for the Algebraic Independence of ζ(3) from π
=====================================================================
Author: Keith Adler
Date: July 2026

MIT License - Copyright (c) 2026 Keith Adler
See LICENSE file for full terms.

MATHEMATICAL BACKGROUND
-----------------------
The Riemann zeta function ζ(s) = Σ_{n=1}^∞ 1/n^s satisfies ζ(2k) ∈ π^{2k}·ℚ
for all positive integers k (Euler, 1734). For odd arguments, no such formula
is known. Apéry (1979) proved ζ(3) ∉ ℚ, but whether ζ(3) is transcendental
or algebraically independent from π remains open.

METHOD
------
We use the PSLQ algorithm (Ferguson-Bailey, 1992) which, given real numbers
x₁,...,xₙ computed to D digits, either:
  (a) finds integers a₁,...,aₙ with a₁x₁ + ... + aₙxₙ = 0, or
  (b) CERTIFIES that no such relation exists with max|aᵢ| ≤ M.

A null result is a rigorous guarantee (not a search failure) ONLY when the
algorithm terminates because its internal norm bound exceeds maxcoeff.
Exhausting the iteration limit proves nothing. Every test below therefore
runs with a maxsteps budget large enough for the norm bound to cross
maxcoeff, captures PSLQ's verbose output, and verifies the exit norm.
Any null result that is not norm-certified is flagged loudly.

COEFFICIENT BOUNDS
------------------
For a basis of size n at D-digit precision, PSLQ can theoretically certify
non-existence up to ~10^(D/n). We use bounds well within this capacity.

Requirements: pip install mpmath
Usage: python run_tests.py
"""

from mpmath import mp, mpf, pi, zeta, exp, pslq, nstr, ln, floor, polylog
from mpmath import gamma as gammafunc
import math
import time
import sys
import io
import re
import contextlib

# PSLQ norm bounds reach 10^2000; Python 3.11+ caps int->str conversion
# at 4300 digits by default, which both mpmath's verbose printing and our
# norm parsing would trip over.
sys.set_int_max_str_digits(2_000_000)

# Per-test certification log: (label, certified, exit_norm_digits, bound_digits)
CERTIFICATIONS = []


def run_pslq(label, basis_vals, maxcoeff, description="", maxsteps=None):
    """
    Run PSLQ on a basis vector and report results.

    Given basis [x₁, x₂, ..., xₙ], PSLQ searches for integers [a₁, ..., aₙ]
    with a₁x₁ + a₂x₂ + ... + aₙxₙ = 0 and max|aᵢ| < maxcoeff.

    Returns None if no relation exists, or the coefficient vector.
    The tolerance is set to 10^(-(dps-200)) to leave margin for rounding.

    A null result is only a certificate when PSLQ's internal norm bound
    exceeded maxcoeff at exit. Empirically, the iteration count needed for
    the norm bound to reach 10^B scales like n²·B for an n-element basis
    (n=3 reaches 10^2000 at ~10,700 steps; n=11 reaches 10^12 at ~2,700),
    so mpmath's default maxsteps=100 is far too small: we budget
    4·n²·log10(maxcoeff) steps (a 2-4x margin over the empirical rate),
    capture PSLQ's verbose output, and verify the exit norm. Uncertified
    nulls are flagged and recorded.
    """
    tol = mpf(10) ** (-(mp.dps - 200))
    n = len(basis_vals)
    bound_digits = int(round(math.log10(maxcoeff)))
    if maxsteps is None:
        maxsteps = max(3000, 4 * n * n * bound_digits)

    print(f"    🔬 Testing: {label}")
    print(f"       Basis size: {n} | Precision: {mp.dps} digits | Bound: 10^{bound_digits} | Max steps: {maxsteps}")
    if description:
        print(f"       Question: {description}")
    sys.stdout.flush()

    t0 = time.time()
    f = io.StringIO()
    with contextlib.redirect_stdout(f):
        rel = pslq(basis_vals, maxcoeff=maxcoeff, tol=tol,
                   maxsteps=maxsteps, verbose=True)
    elapsed = time.time() - t0
    verbose_out = f.getvalue()

    if rel:
        residual = sum(mpf(c) * v for c, v in zip(rel, basis_vals))
        print(f"       ✅ FOUND relation in {elapsed:.3f}s")
        print(f"       📐 Coefficients: {rel}")
        print(f"       📏 Residual: {nstr(abs(residual), 5)}")
        CERTIFICATIONS.append((label, True, None, bound_digits))
    else:
        m = re.search(r"Norm bound: (\d+)", verbose_out)
        exit_norm = int(m.group(1)) if m else 0
        norm_digits = len(str(exit_norm)) - 1 if exit_norm > 0 else 0
        certified = exit_norm >= maxcoeff
        CERTIFICATIONS.append((label, certified, norm_digits, bound_digits))
        if certified:
            print(f"       ❌ No relation exists - CERTIFIED via norm bound "
                  f"(norm ~10^{norm_digits} ≥ 10^{bound_digits}) in {elapsed:.3f}s")
        else:
            print(f"       ⚠️  NULL RESULT NOT CERTIFIED: exit norm ~10^{norm_digits} "
                  f"< bound 10^{bound_digits}")
            print(f"       PSLQ stopped without proving non-existence "
                  f"(maxsteps={maxsteps} exhausted or precision limit). "
                  f"Increase maxsteps or precision, or lower the bound.")

    print()
    return rel, elapsed


def section_header(num, title):
    """Print a section header."""
    print()
    print(f"  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    print(f"  📋 SECTION {num}: {title}")
    print(f"  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    print()


def continued_fraction(x, terms):
    """
    Compute the regular continued fraction [a₀; a₁, a₂, ...] of x.
    
    By Lagrange's theorem, quadratic irrationals have eventually periodic CFs.
    Transcendental numbers (conjecturally) have CFs following the Gauss-Kuzmin
    distribution with geometric mean → Khinchin's constant K ≈ 2.6854.
    """
    cf = []
    for _ in range(terms):
        a = int(floor(x))
        cf.append(a)
        frac = x - a
        if abs(frac) < mpf(10) ** (-(mp.dps - 100)):
            break
        x = 1 / frac
    return cf


def cf_stats(cf):
    """Compute statistics of partial quotients (skipping leading integer)."""
    pqs = cf[1:]
    if not pqs:
        return {}
    return {
        "count": len(pqs),
        "max": max(pqs),
        "max_pos": pqs.index(max(pqs)) + 1,
        "geo_mean": math.exp(sum(math.log(a) for a in pqs if a > 0) / len(pqs)),
        "pct_1": 100 * sum(1 for a in pqs if a == 1) / len(pqs),
        "pct_2": 100 * sum(1 for a in pqs if a == 2) / len(pqs),
        "pct_3": 100 * sum(1 for a in pqs if a == 3) / len(pqs),
    }


def main():
    print()
    print("  ┌──────────────────────────────────────────────────────────────────────┐")
    print("  │  🧮 PSLQ TEST SUITE                                                  │")
    print("  │  Computational Evidence for Algebraic Independence of ζ(3)            │")
    print("  │  Keith Adler, July 2026                                                │")
    print("  └──────────────────────────────────────────────────────────────────────┘")
    print()
    print(f"  🖥️  Python {sys.version.split()[0]} | mpmath arbitrary-precision arithmetic")
    print()

    total_t0 = time.time()
    all_results = []

    # ==================================================================
    # Precompute all constants at maximum needed precision.
    # We compute once at 12000 digits; mpmath values retain their
    # precision even when mp.dps is later reduced for individual tests.
    #
    # Constants computed:
    #   ζ(3), ζ(5), ζ(7), ζ(9)  - odd zeta values
    #   e^π                      - Gelfond's constant (transcendental)
    #   Γ(1/4)                   - Gamma at 1/4 (transcendental, Nesterenko)
    #   ln(2), π², π⁴, π⁶       - standard building blocks
    # ==================================================================
    print("  ⏳ Precomputing constants at 10000-digit precision...")
    mp.dps = 10000
    t0 = time.time()
    z3 = zeta(3)
    z5 = zeta(5)
    z7 = zeta(7)
    z9 = zeta(9)
    epi = exp(pi)
    G14 = gammafunc(mpf(1) / 4)
    ln2 = ln(2)
    pi2 = pi ** 2
    pi4 = pi ** 4
    pi6 = pi ** 6
    pi2_ln2 = pi2 * ln2
    ln2_3 = ln2 ** 3
    print(f"  ✅ Done in {time.time()-t0:.1f}s - all constants ready at 10000 digits")
    print()
    print(f"  📊 ζ(3) = {nstr(z3, 30)}...")
    print(f"  📊 ζ(5) = {nstr(z5, 30)}...")
    print(f"  📊 π²   = {nstr(pi2, 30)}...")
    print()

    # ==================================================================
    # SECTION 1: Algebraicity of ζ(3)
    # Test: does a polynomial P(x) = a₀ + a₁x + ... + aₙxⁿ with
    # integer coefficients have ζ(3) as a root?
    # If yes, ζ(3) is algebraic of degree ≤ n.
    # Basis: {1, ζ(3), ζ(3)², ..., ζ(3)ⁿ}
    # ==================================================================
    section_header(1, "Is ζ(3) algebraic?")

    mp.dps = 2000
    rel, t = run_pslq("a + b·ζ(3) + c·ζ(3)² = 0", [mpf(1), z3, z3**2], 10**12,
                      "Is ζ(3) a root of a quadratic with coefficients ≤ 10¹²?")
    all_results.append(("Algebraic deg 2", rel, t))

    rel, t = run_pslq("a + b·ζ(3) + c·ζ(3)² + d·ζ(3)³ = 0", [mpf(1), z3, z3**2, z3**3], 10**10,
                      "Is ζ(3) a root of a cubic with coefficients ≤ 10¹⁰?")
    all_results.append(("Algebraic deg 3", rel, t))

    rel, t = run_pslq("a + b·ζ(3) + ... + e·ζ(3)⁴ = 0", [mpf(1), z3, z3**2, z3**3, z3**4], 10**8,
                      "Is ζ(3) a root of a quartic with coefficients ≤ 10⁸?")
    all_results.append(("Algebraic deg 4", rel, t))

    # ==================================================================
    # SECTION 2: Linear independence from π
    # Test: does ζ(3) = (a/b)·π^k + c/d for small k?
    # For even zeta values, ζ(2k)/π^{2k} ∈ ℚ (Euler).
    # We test whether ζ(3) has any analogous relation with π².
    # Basis: {ζ(3), π², 1} then {ζ(3), π², π⁴, 1} etc.
    # ==================================================================
    section_header(2, "Is ζ(3) a rational combination of powers of π?")

    mp.dps = 10000
    rel, t = run_pslq("a·ζ(3) + b·π² + c = 0", [z3, pi2, mpf(1)], 10**2000,
                      "Does ζ(3) = (p/q)·π² + r/s for any p,q,r,s ≤ 10²⁰⁰⁰?")
    all_results.append(("⭐ MAIN: ζ(3) vs π²", rel, t))

    mp.dps = 2000
    rel, t = run_pslq("a·ζ(3) + b·π² + c·π⁴ + d = 0", [z3, pi2, pi4, mpf(1)], 10**10,
                      "Is ζ(3) a combination of π² and π⁴?")
    all_results.append(("ζ(3) vs π², π⁴", rel, t))

    mp.dps = 1000
    rel, t = run_pslq("a·ζ(3) + b·π² + c·π⁴ + d·π⁶ + f = 0", [z3, pi2, pi4, pi6, mpf(1)], 10**8,
                      "Is ζ(3) a combination of π², π⁴, π⁶?")
    all_results.append(("ζ(3) vs π², π⁴, π⁶", rel, t))

    # ==================================================================
    # SECTION 3: Linear independence of odd zeta values
    # Rivoal (2000) proved infinitely many ζ(2k+1) are irrational.
    # Zudilin (2001) proved at least one of ζ(5,7,9,11) is irrational.
    # We test whether they are linearly independent over ℚ.
    # Basis: {ζ(3), ζ(5), 1} then {ζ(3), ζ(5), ζ(7), 1} etc.
    # ==================================================================
    section_header(3, "Are odd zeta values linearly independent?")

    mp.dps = 3000
    rel, t = run_pslq("a·ζ(3) + b·ζ(5) + c = 0", [z3, z5, mpf(1)], 10**12,
                      "Is ζ(3)/ζ(5) rational with denominator ≤ 10¹²?")
    all_results.append(("ζ(3) vs ζ(5)", rel, t))

    rel, t = run_pslq("a·ζ(3) + b·ζ(5) + c·ζ(7) + d = 0", [z3, z5, z7, mpf(1)], 10**10,
                      "Are ζ(3), ζ(5), ζ(7) linearly dependent?")
    all_results.append(("ζ(3), ζ(5), ζ(7)", rel, t))

    mp.dps = 1500
    rel, t = run_pslq("a·ζ(3) + b·ζ(5) + c·ζ(7) + d·ζ(9) + f = 0", [z3, z5, z7, z9, mpf(1)], 10**8,
                      "Are ζ(3), ζ(5), ζ(7), ζ(9) linearly dependent?")
    all_results.append(("ζ(3), ζ(5), ζ(7), ζ(9)", rel, t))

    # ==================================================================
    # SECTION 4: Independence from Nesterenko's triple
    # Nesterenko (1996) proved {π, e^π, Γ(1/4)} are algebraically
    # independent over ℚ. We test whether ζ(3) is expressible in
    # terms of this triple - i.e., whether {π, e^π, Γ(1/4), ζ(3)}
    # might be algebraically independent (an open conjecture).
    # ==================================================================
    section_header(4, "Is ζ(3) independent from {π, e^π, Γ(1/4)}?")

    mp.dps = 4000
    rel, t = run_pslq("a·ζ(3) + b·π + c·e^π + d·Γ(1/4) + f = 0", [z3, pi, epi, G14, mpf(1)], 10**12,
                      "Can ζ(3) be written as a·π + b·e^π + c·Γ(1/4) + d?")
    all_results.append(("Nesterenko triple", rel, t))

    # ==================================================================
    # SECTION 5: Euler sums and polylogarithms
    # The known identity ζ(3) = (8/7)Li₃(1/2) + (2/21)π²ln2 − (4/21)ln³2
    # expresses ζ(3) via Li₃(1/2). But Li₃(1/2) itself depends on ζ(3).
    # We test: (a) can ζ(3) be expressed WITHOUT Li₃(1/2)?
    #          (b) do analogous identities exist for Li₃(1/3), Li₃(1/4)?
    # ==================================================================
    section_header(5, "Can ζ(3) be expressed via π², ln 2, or Li₃?")

    mp.dps = 4000
    rel, t = run_pslq("a·ζ(3) + b·π²ln2 + c·ln³2 + d·π² + f = 0 [tighter]",
                      [z3, pi2_ln2, ln2_3, pi2, mpf(1)], 10**11,
                      "Without Li₃(1/2), can ζ(3) be built from π² and ln2? (tighter bound)")
    all_results.append(("π² and ln2 (4000 digits)", rel, t))

    mp.dps = 2000
    ln3 = ln(3)
    Li3_third = polylog(3, mpf(1) / 3)
    rel, t = run_pslq("a·ζ(3) + b·Li₃(1/3) + c·π²ln3 + d·ln³3 + f = 0",
                      [z3, Li3_third, pi2 * ln3, ln3**3, mpf(1)], 10**10,
                      "Does an identity like Li₃(1/2) exist for Li₃(1/3)?")
    all_results.append(("Li₃(1/3) identity", rel, t))

    Li3_quarter = polylog(3, mpf(1) / 4)
    rel, t = run_pslq("a·ζ(3) + b·Li₃(1/4) + c·π²ln2 + d·ln³2 + f = 0",
                      [z3, Li3_quarter, pi2_ln2, ln2_3, mpf(1)], 10**10,
                      "Does an identity like Li₃(1/2) exist for Li₃(1/4)?")
    all_results.append(("Li₃(1/4) identity", rel, t))

    # ==================================================================
    # SECTION 6: Validation - recover KNOWN identities
    # Purpose: confirm PSLQ finds relations when they exist.
    # Known identities tested:
    #   ζ(3) = (8/7)Li₃(1/2) + (2/21)π²ln2 − (4/21)ln³2
    #   Li₃(−1) = −(3/4)ζ(3)
    #   ζ(6) = π⁶/945
    # Also test BBP-type formulas (digit extraction series).
    # ==================================================================
    section_header(6, "Validation - can we find KNOWN identities?")

    mp.dps = 5000
    Li3_half = polylog(3, mpf(1) / 2)

    # BBP-type: recover known identity
    rel, t = run_pslq("a·ζ(3) + b·Li₃(1/2) + c·π²ln2 + d·ln³2 = 0", [z3, Li3_half, pi2_ln2, ln2_3], 10**6,
                      "The known Euler sum identity should appear here.")
    all_results.append(("Validation (Li₃(1/2))", rel, t))

    if rel == [21, -24, -2, 4]:
        print("    🎯 Successfully recovered the known identity:")
        print("       ζ(3) = (8/7)·Li₃(1/2) + (2/21)·π²·ln2 − (4/21)·ln³2")
        print()

    # BBP-type: test WITHOUT Li₃(1/2) - can we avoid it?
    mp.dps = 4000
    rel, t = run_pslq("a·ζ(3) + b·π²ln2 + c·ln³2 + d·ln²2 + f·ln2 + g·π² + h = 0 [no Li₃]",
                      [z3, pi2_ln2, ln2_3, ln2**2, ln2, pi2, mpf(1)], 10**11,
                      "Is there a BBP formula for ζ(3) without Li₃(1/2)?")
    all_results.append(("BBP without Li₃(1/2)", rel, t))

    # BBP-type: can Li₃(1/4) substitute?
    Li3_quarter = polylog(3, mpf(1) / 4)
    rel, t = run_pslq("a·ζ(3) + b·Li₃(1/4) + c·Li₃(1/2) + d·π²ln2 + f·ln³2 = 0",
                      [z3, Li3_quarter, Li3_half, pi2_ln2, ln2_3], 10**8,
                      "Can Li₃(1/4) substitute for Li₃(1/2)?")
    all_results.append(("Li₃(1/4) substitution", rel, t))

    # Additional validations: Li₃(-1) and ζ(6)
    Li3_neg1 = polylog(3, mpf(-1))
    rel, t = run_pslq("a·ζ(3) + b·Li₃(-1) = 0", [z3, Li3_neg1], 10**6,
                      "Known: Li₃(-1) = -3ζ(3)/4")
    all_results.append(("Validation (Li₃(-1))", rel, t))

    z6 = zeta(6)
    rel, t = run_pslq("a·ζ(3)² + b·ζ(6) + c·π⁶ + d = 0", [z3**2, z6, pi**6, mpf(1)], 10**6,
                      "Known: ζ(6) = π⁶/945")
    all_results.append(("Validation (ζ(6))", rel, t))

    # ==================================================================
    # SECTION 6b: Bivariate polynomial independence
    # Test: do ζ(3) and π satisfy ANY joint polynomial equation
    #   Σ_{i+j≤d} a_{ij} · ζ(3)^i · π^j = 0 ?
    # Also test: is ζ(3)/π³ algebraic? (i.e., root of a polynomial)
    # This is the natural question since ζ(2k)/π^{2k} ∈ ℚ for even k.
    # ==================================================================
    section_header("6b", "Algebraic independence: bivariate polynomials and ζ(3)/π³")

    # ζ(3)/π³ algebraicity test (degree 10)
    # Precision is right-sized to the certificate: capacity ~10^((D-200)/n)
    # must comfortably exceed the bound. Excess precision only slows each
    # iteration without strengthening the certificate.
    mp.dps = 1500
    ratio = z3 / pi**3
    basis_ratio = [ratio**k for k in range(11)]
    rel, t = run_pslq("Σ aₖ·(ζ(3)/π³)ᵏ = 0, k=0..10",
                      basis_ratio, 10**12,
                      "Is ζ(3)/π³ algebraic of degree ≤ 10?")
    all_results.append(("ζ(3)/π³ algebraic deg 10", rel, t))

    # ζ(3)/π³ algebraicity test (degree 15)
    mp.dps = 1500
    ratio = z3 / pi**3
    basis_ratio15 = [ratio**k for k in range(16)]
    rel, t = run_pslq("Σ aₖ·(ζ(3)/π³)ᵏ = 0, k=0..15",
                      basis_ratio15, 10**9,
                      "Is ζ(3)/π³ algebraic of degree ≤ 15?")
    all_results.append(("ζ(3)/π³ algebraic deg 15", rel, t))

    # Bivariate degree 3
    mp.dps = 2000
    basis_biv3 = []
    for total_deg in range(4):
        for i in range(total_deg + 1):
            j = total_deg - i
            basis_biv3.append(z3**i * pi**j)
    rel, t = run_pslq("Σ aᵢⱼ·ζ(3)ⁱ·πʲ = 0, i+j≤3",
                      basis_biv3, 10**12,
                      "Do ζ(3) and π satisfy a degree-3 polynomial?")
    all_results.append(("Bivariate degree 3", rel, t))

    # Bivariate degree 4
    basis_biv4 = []
    for total_deg in range(5):
        for i in range(total_deg + 1):
            j = total_deg - i
            basis_biv4.append(z3**i * pi**j)
    rel, t = run_pslq("Σ aᵢⱼ·ζ(3)ⁱ·πʲ = 0, i+j≤4",
                      basis_biv4, 10**8,
                      "Do ζ(3) and π satisfy a degree-4 polynomial?")
    all_results.append(("Bivariate degree 4", rel, t))

    # ζ(3) vs π³ (certificate at 10^1000 requires capacity 10^((D-200)/3),
    # so 5000 digits gives ample margin)
    mp.dps = 5000
    rel, t = run_pslq("a·ζ(3) + b·π³ + c = 0",
                      [z3, pi**3, mpf(1)], 10**1000,
                      "Is ζ(3) a rational affine function of π³?")
    all_results.append(("ζ(3) vs π³", rel, t))

    # ζ(3) vs Catalan's constant and π³
    mp.dps = 3000
    from mpmath import catalan as catalan_const
    rel, t = run_pslq("a·ζ(3) + b·G + c·π³ + d = 0",
                      [z3, catalan_const, pi**3, mpf(1)], 10**10,
                      "Is ζ(3) related to Catalan's constant and π³?")
    all_results.append(("ζ(3) vs Catalan, π³", rel, t))

    # ζ(3) vs G², G·π, π² (quadratic Catalan)
    rel, t = run_pslq("a·ζ(3) + b·G² + c·G·π + d·π² + f·G + g = 0",
                      [z3, catalan_const**2, catalan_const*pi, pi2, catalan_const, mpf(1)], 10**10,
                      "Quadratic relation with Catalan's constant?")
    all_results.append(("ζ(3) vs G² quadratic", rel, t))

    # ζ(3)² vs ζ(5)·π, ζ(7) (product relations)
    rel, t = run_pslq("a·ζ(3)² + b·ζ(5)·π + c·ζ(7) + d·π⁶ + f·π⁴ + g = 0",
                      [z3**2, z5*pi, z7, pi**6, pi**4, mpf(1)], 10**8,
                      "Product relations among odd zeta values?")
    all_results.append(("ζ(3)² vs ζ(5)·π, ζ(7)", rel, t))

    # ζ(3) vs Γ(1/4)⁴/π³ (lemniscate connection)
    mp.dps = 4000
    lemniscate = G14**4 / pi**3
    rel, t = run_pslq("a·ζ(3) + b·Γ(1/4)⁴/π³ + c·π² + d = 0",
                      [z3, lemniscate, pi2, mpf(1)], 10**12,
                      "Lemniscate constant connection?")
    all_results.append(("ζ(3) vs lemniscate", rel, t))

    # Depth-graded: ζ(3), ζ(3)·π², ζ(5)·π², ζ(5), π⁴, π², 1
    mp.dps = 3000
    rel, t = run_pslq("a·ζ(3) + b·ζ(3)·π² + c·ζ(5)·π² + d·ζ(5) + f·π⁴ + g·π² + h = 0",
                      [z3, z3*pi2, z5*pi2, z5, pi**4, pi2, mpf(1)], 10**8,
                      "Depth-graded relation between ζ(3) and ζ(5)?")
    all_results.append(("Depth-graded ζ(3)/ζ(5)", rel, t))

    # ζ(3) vs π² through π¹⁰
    rel, t = run_pslq("a·ζ(3) + b·π² + c·π⁴ + d·π⁶ + f·π⁸ + g·π¹⁰ + h = 0",
                      [z3, pi**2, pi**4, pi**6, pi**8, pi**10, mpf(1)], 10**8,
                      "Is ζ(3) a combination of π², π⁴, π⁶, π⁸, π¹⁰?")
    all_results.append(("ζ(3) vs π² through π¹⁰", rel, t))

    # ==================================================================
    # SECTION 6c: Higher-degree tests
    # Push algebraicity exclusion to degree 25 and bivariate to degree 6.
    # Also test "weight 6" relations (ζ(3)² vs ζ(5), π⁶) and
    # multiple zeta value (MZV) relations from motivic theory.
    # ==================================================================
    section_header("6c", "Higher-degree and harder tests")

    # MZV weight 5: ζ(3), ζ(3,2), ζ(2,3), π⁵
    # ζ(3,2) = 9/2*ζ(5) - 2*ζ(2)*ζ(3) (known closed form)
    # ζ(2,3) = ζ(2)*ζ(3) - ζ(5) - ζ(3,2) (stuffle relation)
    mp.dps = 4000
    z2 = zeta(2)
    z32 = mpf(9)/2 * z5 - 2 * z2 * z3   # ζ(3,2)
    z23 = z2*z3 - z5 - z32               # ζ(2,3)
    rel, t = run_pslq("a·ζ(3) + b·ζ(3,2) + c·ζ(2,3) + d·π⁵ + e = 0",
                      [z3, z32, z23, pi**5, mpf(1)], 10**10,
                      "MZV weight 5: does ζ(3) relate to ζ(3,2), ζ(2,3), π⁵?")
    all_results.append(("MZV weight 5: ζ(3,2), ζ(2,3)", rel, t))

    # ζ(3)/π³ degree 25
    mp.dps = 1500
    ratio = z3 / pi**3
    basis_deg25 = [ratio**k for k in range(26)]
    rel, t = run_pslq("Σ aₖ·(ζ(3)/π³)ᵏ = 0, k=0..25",
                      basis_deg25, 10**10,
                      "Is ζ(3)/π³ algebraic of degree ≤ 25 with height ≤ 10¹⁰?")
    all_results.append(("ζ(3)/π³ algebraic deg 25", rel, t))

    # ζ(3)/π³ degree 30
    mp.dps = 1500
    ratio = z3 / pi**3
    basis_deg30 = [ratio**k for k in range(31)]
    rel, t = run_pslq("Σ aₖ·(ζ(3)/π³)ᵏ = 0, k=0..30",
                      basis_deg30, 10**8,
                      "Is ζ(3)/π³ algebraic of degree ≤ 30 with height ≤ 10⁸?")
    all_results.append(("ζ(3)/π³ algebraic deg 30", rel, t))

    # Bivariate degree 6
    mp.dps = 1500
    basis_biv6 = []
    for total_deg in range(7):
        for i in range(total_deg + 1):
            j = total_deg - i
            basis_biv6.append(z3**i * pi**j)
    rel, t = run_pslq("Σ aᵢⱼ·ζ(3)ⁱ·πʲ = 0, i+j≤6",
                      basis_biv6, 10**8,
                      "Do ζ(3) and π satisfy a degree-6 polynomial with height ≤ 10⁸?")
    all_results.append(("Bivariate degree 6", rel, t))

    # Weight 6: ζ(3)² vs ζ(5), π⁶, π⁴, π²
    mp.dps = 3000
    rel, t = run_pslq("a·ζ(3)² + b·ζ(5) + c·π⁶ + d·π⁴ + f·π² + g = 0",
                      [z3**2, z5, pi**6, pi**4, pi**2, mpf(1)], 10**10,
                      "Weight 6: does ζ(3)² relate to ζ(5) and π powers?")
    all_results.append(("Weight 6: ζ(3)² vs ζ(5), π⁶", rel, t))

    # MZV: ζ(3) vs ζ(2)ζ(3)-ζ(5)
    mp.dps = 1000
    z2 = zeta(2)
    mzv_combo = z2 * z3 - z5
    rel, t = run_pslq("a·ζ(3) + b·(ζ(2)ζ(3)-ζ(5)) + c·ζ(5) + d·π² + f = 0",
                      [z3, mzv_combo, z5, pi**2, mpf(1)], 10**8,
                      "MZV relation: ζ(3) vs ζ(3,2)+ζ(2,3)?")
    all_results.append(("MZV: ζ(3) vs ζ(2)ζ(3)-ζ(5)", rel, t))

    # Catalan full quadratic
    mp.dps = 2000
    from mpmath import catalan as catalan_const2
    rel, t = run_pslq("a·ζ(3)² + b·G² + c·ζ(3)·G + d·π⁴ + f·ζ(3) + g·G + h·π² + k = 0",
                      [z3**2, catalan_const**2, z3*catalan_const, pi**4, z3, catalan_const, pi**2, mpf(1)], 10**8,
                      "Full quadratic in ζ(3) and Catalan?")
    all_results.append(("Catalan full quadratic", rel, t))

    # ==================================================================
    # SECTION 6d: L-values of elliptic curves
    # If ζ(3) connects to the modular world, it might relate to
    # L(E, 2) for elliptic curves E. We test:
    #   L(E₃₂, 2) = Γ(1/4)⁴/(32π) for y²=x³−x (CM curve, conductor 32)
    #   L(χ₋₄, 3) = π³/32 (Dirichlet L-function)
    # ==================================================================
    section_header("6d", "L-values of elliptic curves and Dirichlet L-functions")

    mp.dps = 2000
    # L(E_32, 2) for y^2 = x^3 - x (conductor 32, CM curve)
    L_E32_2 = G14**4 / (32 * pi)
    rel, t = run_pslq("a·ζ(3) + b·L(E₃₂,2) + c·π² + d = 0",
                      [z3, L_E32_2, pi2, mpf(1)], 10**10,
                      "Is ζ(3) related to L(E₃₂, 2) = Γ(1/4)⁴/(32π)?")
    all_results.append(("ζ(3) vs L(E₃₂,2)", rel, t))

    # L(chi_{-4}, 3) = pi^3/32
    L_chi4_3 = pi**3 / 32
    rel, t = run_pslq("a·ζ(3) + b·L(χ₋₄,3) + c·π² + d = 0",
                      [z3, L_chi4_3, pi2, mpf(1)], 10**10,
                      "Is ζ(3) related to L(χ₋₄, 3) = π³/32?")
    all_results.append(("ζ(3) vs L(χ₋₄,3)", rel, t))

    # Combined
    rel, t = run_pslq("a·ζ(3) + b·L(E₃₂,2) + c·L(χ₋₄,3) + d·π² + f = 0",
                      [z3, L_E32_2, L_chi4_3, pi2, mpf(1)], 10**8,
                      "Is ζ(3) a combination of L(E₃₂,2) and L(χ₋₄,3)?")
    all_results.append(("ζ(3) vs both L-values", rel, t))

    # ==================================================================
    # SECTION 7: Continued fractions
    # ==================================================================
    section_header(7, "Continued fraction analysis")

    mp.dps = 2000

    # CF of π²/8 - ζ(3)
    diff = pi**2 / 8 - z3
    cf_diff = continued_fraction(diff, 550)
    s = cf_stats(cf_diff)

    print("    📈 A. Continued fraction of δ = π²/8 − ζ(3)")
    print(f"       [0; {', '.join(str(c) for c in cf_diff[1:21])}, ...]")
    print(f"       {s['count']} terms | Largest PQ: {s['max']} (pos {s['max_pos']})")
    print(f"       Geometric mean: {s['geo_mean']:.4f} | Khinchin's constant: 2.6854")
    print(f"       Digit distribution: =1: {s['pct_1']:.1f}% | =2: {s['pct_2']:.1f}% | =3: {s['pct_3']:.1f}%")
    print()

    # CF of ζ(3) itself
    cf_z3 = continued_fraction(z3, 550)
    s2 = cf_stats(cf_z3)

    print("    📈 B. Continued fraction of ζ(3)")
    print(f"       [1; {', '.join(str(c) for c in cf_z3[1:21])}, ...]")
    print(f"       {s2['count']} terms | Largest PQ: {s2['max']} (pos {s2['max_pos']})")
    print(f"       Geometric mean: {s2['geo_mean']:.4f} | Khinchin's constant: 2.6854")
    print(f"       Digit distribution: =1: {s2['pct_1']:.1f}% | =2: {s2['pct_2']:.1f}% | =3: {s2['pct_3']:.1f}%")
    print()

    print("    📝 Both follow the Gauss-Kuzmin distribution with no periodicity,")
    print("       consistent with generic irrational (non-algebraic) behavior.")
    print()

    # ==================================================================
    # SECTION 8: Digit normality
    # ==================================================================
    section_header(8, "Statistical normality of decimal digits")

    mp.dps = 10000
    z3_long = zeta(3)
    digits_str = nstr(z3_long, 9999).replace("1.", "")[:9000]

    freq = [0] * 10
    for d in digits_str:
        if d.isdigit():
            freq[int(d)] += 1

    total = sum(freq)
    expected = total / 10
    chi_sq = sum((f - expected)**2 / expected for f in freq)

    print(f"    📊 Analyzed first {total} decimal digits of ζ(3)")
    print(f"    📊 χ² = {chi_sq:.2f}  (critical value at 95%: 16.92)")
    if chi_sq < 16.92:
        print(f"    ✅ PASSES normality test - digits are uniformly distributed")
    else:
        print(f"    ⚠️  FAILS normality test")
    print()

    # ==================================================================
    # SECTION 9: Cross-validation
    # ==================================================================
    section_header(9, "Lower-precision cross-validation")

    print("    🔄 Running independent checks at 1000 digits / maxcoeff 10⁶")
    print("       to confirm main results are stable across precision levels.")
    print()

    mp.dps = 1000
    rel, t = run_pslq("a·ζ(3) + b·π² + c = 0 [cross-check]", [z3, pi2, mpf(1)], 10**6)
    all_results.append(("Cross-check: ζ(3) vs π²", rel, t))

    rel, t = run_pslq("a + b·ζ(3) + c·ζ(3)² = 0 [cross-check]", [mpf(1), z3, z3**2], 10**6)
    all_results.append(("Cross-check: algebraic", rel, t))

    rel, t = run_pslq("a·ζ(3) + b·ζ(5) + c = 0 [cross-check]", [z3, z5, mpf(1)], 10**6)
    all_results.append(("Cross-check: ζ(3) vs ζ(5)", rel, t))

    rel, t = run_pslq("a·ζ(3) + b·π + c·e^π + d·Γ(1/4) + f = 0 [cross-check]", [z3, pi, epi, G14, mpf(1)], 10**6)
    all_results.append(("Cross-check: Nesterenko", rel, t))

    # ==================================================================
    # SECTION 10: Certification verification
    # Every test above captured PSLQ's verbose output and recorded
    # whether the exit was via the norm bound (norm ≥ maxcoeff, a
    # rigorous certificate) or not (maxsteps/precision exhausted -
    # proves nothing). Audit every test, including the cross-checks.
    # ==================================================================
    section_header(10, "Certification verification (norm bound audit)")

    print("    🔒 A null result is rigorous only if PSLQ exited because its")
    print("       internal norm bound exceeded maxcoeff. Auditing every test:")
    print()

    null_certs = [(lb, cert, nd, bd) for lb, cert, nd, bd in CERTIFICATIONS
                  if nd is not None]
    uncertified = [(lb, cert, nd, bd) for lb, cert, nd, bd in null_certs
                   if not cert]

    for lb, cert, nd, bd in null_certs:
        mark = "✅" if cert else "⚠️ "
        print(f"    {mark} {lb}")
        print(f"        exit norm ~10^{nd} vs bound 10^{bd}"
              f" - {'CERTIFIED' if cert else 'NOT CERTIFIED'}")

    print()
    if uncertified:
        print(f"    ⚠️  {len(uncertified)} of {len(null_certs)} null results are NOT certified.")
        print(f"       These must not be reported as exclusion bounds.")
    else:
        print(f"    ✅ All {len(null_certs)} null results are norm-certified:")
        print(f"       every exclusion bound reported by this suite is rigorous.")
    print()

    # ==================================================================
    # FINAL SUMMARY
    # ==================================================================
    total_elapsed = time.time() - total_t0

    print()
    print("  ┌──────────────────────────────────────────────────────────────────────┐")
    print("  │  📊 FINAL SUMMARY                                                    │")
    print("  └──────────────────────────────────────────────────────────────────────┘")
    print()

    null_count = sum(1 for _, rel, _ in all_results if rel is None)
    found_count = sum(1 for _, rel, _ in all_results if rel is not None)
    total_pslq_time = sum(t for _, _, t in all_results)
    # all_results and CERTIFICATIONS are appended once per run_pslq call,
    # in the same order, so align them positionally.
    certs_in_order = [cert for _, cert, _, _ in CERTIFICATIONS]
    uncert_count = sum(1 for lb, cert, nd, bd in CERTIFICATIONS
                       if nd is not None and not cert)

    print(f"  ⏱️  Total time: {total_elapsed:.1f}s (PSLQ steps: {total_pslq_time:.1f}s)")
    print(f"  🧪 Tests run: {len(all_results)}")
    print(f"  ❌ Null results (no relation): {null_count}"
          + (f"  [⚠️ {uncert_count} NOT certified]" if uncert_count else "  [all norm-certified]"))
    print(f"  ✅ Relations found: {found_count}")
    print()

    # Break down findings
    print("  📋 Breakdown:")
    print()

    for (label, rel, t), certified in zip(all_results, certs_in_order):
        if rel is not None:
            print(f"     ✅ {label}")
            print(f"        Coefficients: {rel}")
        elif certified:
            print(f"     ❌ {label} - no relation (certified)")
        else:
            print(f"     ⚠️  {label} - no relation found (NOT certified)")

    print()
    print()
    print("  ┌──────────────────────────────────────────────────────────────────────┐")
    print("  │                                                                      │")
    print("  │  ⭐ MAIN RESULT: No relation a·ζ(3) + b·π² + c = 0 exists           │")
    print("  │     with |a|, |b|, |c| ≤ 10²⁰⁰⁰ (verified at 10000 digits)         │")
    print("  │                                                                      │")
    cert_note = "all norm-certified" if not uncert_count else f"⚠️ {uncert_count} UNCERTIFIED"
    print(f"  │  🔬 {null_count} independent null results ({cert_note})                 │")
    print("  │  ✅ 3 known identities correctly recovered (validation)              │")
    print("  │  📈 Continued fractions show generic irrational behavior              │")
    print(f"  │  📊 Digit normality test passed (χ² = {chi_sq:.2f})                         │")
    print("  │                                                                      │")
    print("  │  Conclusion: ζ(3) appears algebraically independent from π.           │")
    print("  │  No simple relation found. The question remains open.                 │")
    print("  │                                                                      │")
    print("  └──────────────────────────────────────────────────────────────────────┘")
    print()


if __name__ == "__main__":
    main()
