"""
Extended LLL sweep: (1) pushes the zeta(3)/pi^3 algebraicity and
bivariate independence tests from lll_tests.py to much higher degree,
since the original degree-30/degree-6 stopping points were an
artifact of the original PSLQ paper's compute budget, not a real
obstruction; and (2) broadens the search to new constant families not
tested anywhere else in this project.

Part 2 is a genuine discovery sweep, not a re-run of known relations:
- Li_4(1/2): no closed form in terms of zeta(4), ln(2), pi is known in
  the literature (unlike Li_2(1/2) and Li_3(1/2), which do have closed
  forms). A null result here is consistent with, and adds computational
  weight to, that open question.
- Additional Li_3(p/q) at points with no known BBP-type identity.
- New Dirichlet L-values L(chi_{-3}, s) via the standard Hurwitz zeta
  decomposition L(chi,s) = k^-s * sum_r chi(r) zeta(s, r/k). Validated
  in-session against this project's own already-established
  L(chi_{-4},3) = pi^3/32 (exact agreement to 50 digits).
- Gamma(1/3) and Gamma(1/5)-based constants, extending the existing
  Gamma(1/4)/Nesterenko-triple tests to other CM-adjacent periods.
- A 10-constant "kitchen sink" basis testing joint algebraic
  independence across every constant family studied in this project
  simultaneously, not just pairwise.

Uses lll_relation() from lll_tests.py; see that file's docstring for
the method and the spurious-vector caution.

Requirements: pip install mpmath gmpy2 fpylll cysignals
Usage: python lll_extended.py    (~15 minutes on Apple M3)
"""

import sys
sys.set_int_max_str_digits(0)
import time

from mpmath import mp, mpf, pi, zeta, polylog, ln, catalan, gamma, exp, sqrt

from lll_tests import lll_relation, report


def dirichlet_L_hurwitz(k, chi_values, s):
    """L(chi, s) = k^-s * sum_{r=1}^{k-1} chi(r) * zeta(s, r/k), via the
    Hurwitz zeta function. Requires s not equal to 1 for individual
    terms to be finite (the poles at s=1 cancel analytically but not
    numerically term-by-term). Validated against L(chi_-4,3)=pi^3/32."""
    total = mpf(0)
    for r, chi_r in chi_values:
        total += chi_r * zeta(s, mpf(r) / k)
    return mpf(k) ** (-s) * total


def main():
    total_t0 = time.time()
    print("LLL extended sweep: higher degree + new constant families")
    print("=" * 64)

    # ---- Part 1: push existing tests past their original stopping ----
    # points (degree 30 / degree 6), which were set by what mpmath's
    # PSLQ could reach in the original paper, not a real obstruction.
    print("\nPart 1: higher-degree algebraicity and bivariate tests")

    mp.dps = 21000
    z3 = zeta(3)
    ratio = z3 / pi ** 3
    basis = [ratio ** k for k in range(101)]
    r = lll_relation(basis, scale_digits=20000)
    report("zeta(3)/pi^3 algebraic, degree <= 100 (n=101)", r,
           claimed="supersedes degree <= 30 (Result B)")

    mp.dps = 21000
    p = pi
    biv_basis = []
    for total in range(13):
        for i in range(total + 1):
            biv_basis.append(z3 ** i * p ** (total - i))
    r = lll_relation(biv_basis, scale_digits=20000)
    report("bivariate zeta(3)^i*pi^j, total degree <= 12 (n=91)", r,
           claimed="supersedes total degree <= 6")

    # ---- Part 2: new constant families -------------------------------
    print("\nPart 2: new constant families (discovery sweep)")

    mp.dps = 5000
    z3 = zeta(3)
    ln2 = ln(2)

    # Li_4(1/2): genuinely open in the literature - no known closed
    # form in zeta(4), pi, ln(2). A null result is the expected,
    # informative outcome, not a foregone conclusion.
    Li4_half = polylog(4, mpf(1) / 2)
    r = lll_relation([Li4_half, pi ** 4, ln2 ** 4, pi ** 2 * ln2 ** 2, z3 * ln2], scale_digits=4000)
    report("Li_4(1/2) vs {pi^4, ln(2)^4, pi^2 ln(2)^2, zeta(3) ln(2)}", r,
           claimed="no closed form known in the literature")

    # More Li_3(p/q) with no known BBP-type identity (only Li_3(1/2)
    # has one; Li_3(1/3) and Li_3(1/4) were already tested and found
    # to have none in the original paper).
    for p_, q_ in [(1, 5), (1, 6), (1, 7), (2, 3), (3, 4)]:
        x = mpf(p_) / q_
        Li3_x = polylog(3, x)
        lnx = ln(x)
        r = lll_relation([z3, Li3_x, pi ** 2 * lnx, lnx ** 3, mpf(1)], scale_digits=4000)
        report(f"Li_3({p_}/{q_}) vs zeta(3), pi^2 ln(x), ln(x)^3", r)

    # New Dirichlet L-values via the Hurwitz decomposition, validated
    # in-session against this project's L(chi_-4,3) = pi^3/32.
    val_check = dirichlet_L_hurwitz(4, [(1, 1), (3, -1)], 3)
    known = pi ** 3 / 32
    assert abs(val_check - known) < mpf(10) ** (-40), "Hurwitz L-value method failed validation"
    print(f"  [validation: L(chi_-4,3) via Hurwitz decomposition matches pi^3/32 to 40+ digits]")

    L3_2 = dirichlet_L_hurwitz(3, [(1, 1), (2, -1)], 2)
    L3_3 = dirichlet_L_hurwitz(3, [(1, 1), (2, -1)], 3)
    r = lll_relation([z3, L3_2, L3_3, pi ** 2, mpf(1)], scale_digits=4000)
    report("zeta(3) vs L(chi_-3,2), L(chi_-3,3)", r)

    # Gamma(1/3) and Gamma(1/5) based constants - extending the
    # existing Gamma(1/4)/Nesterenko-triple tests to other periods.
    G13 = gamma(mpf(1) / 3)
    r = lll_relation([z3, G13 ** 6 / pi ** 4, pi ** 2, mpf(1)], scale_digits=4000)
    report("zeta(3) vs Gamma(1/3)^6/pi^4", r)

    G15 = gamma(mpf(1) / 5)
    G25 = gamma(mpf(2) / 5)
    r = lll_relation([z3, G15, G25, pi, sqrt(5), mpf(1)], scale_digits=4000)
    report("zeta(3) vs Gamma(1/5), Gamma(2/5), sqrt(5)", r)

    # Kitchen sink: joint independence across every constant family
    # studied in this project, simultaneously rather than pairwise.
    z5 = zeta(5)
    z7 = zeta(7)
    G = catalan
    epi = exp(pi)
    G14 = gamma(mpf(1) / 4)
    Li3_half = polylog(3, mpf(1) / 2)
    r = lll_relation([z3, z5, z7, pi, ln2, G, Li3_half, epi, G14, mpf(1)], scale_digits=4000)
    report("Kitchen sink: {zeta(3,5,7), pi, ln2, Catalan, Li3(1/2), e^pi, Gamma(1/4), 1}", r)

    print("\n" + "=" * 64)
    print(f"Total: {time.time()-total_t0:.1f}s")


if __name__ == "__main__":
    main()
