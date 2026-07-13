"""
Extends the LLL search to constant families this paper had not yet
given the same algebraicity/bivariate treatment given to zeta(3)/pi^3:

- Catalan's constant G: is G/pi^2 algebraic? does G satisfy any joint
  polynomial with pi (mirroring the zeta(3),pi bivariate test)?
- zeta(5) and zeta(7): the paper tests LINEAR independence of the odd
  zetas from each other and from pi, but never tests whether zeta(5)
  or zeta(7) individually satisfy an algebraicity-with-pi relation
  the way zeta(3) does (i.e. is zeta(2k+1)/pi^(2k+1) algebraic?), nor
  whether zeta(5) and zeta(7) satisfy a joint (non-linear) polynomial.
- The Euler-Mascheroni constant gamma: does not appear anywhere else
  in this paper. Unlike zeta(3), gamma's irrationality is not even
  proven, so a null result here is a genuinely different flavor of
  open question - we test gamma's algebraicity alone, its linear
  independence from pi and e, and whether gamma, pi, e satisfy any
  joint polynomial relation.

Uses lll_relation() from lll_tests.py; see that file's docstring for
the method and the spurious-vector caution.

Requirements: pip install mpmath gmpy2 fpylll cysignals
Usage: python lll_new_constants.py    (~1 minute on Apple M3)
"""

import sys
sys.set_int_max_str_digits(0)
import time

from mpmath import mp, mpf, pi, catalan, euler, e, zeta

from lll_tests import lll_relation, report


def bivariate_basis(x, y, deg):
    basis = []
    for total in range(deg + 1):
        for i in range(total + 1):
            basis.append(x ** i * y ** (total - i))
    return basis


def trivariate_basis(x, y, z, deg):
    basis = []
    for total in range(deg + 1):
        for i in range(total + 1):
            for j in range(total + 1 - i):
                k = total - i - j
                basis.append(x ** i * y ** j * z ** k)
    return basis


def main():
    total_t0 = time.time()
    print("LLL sweep: new constant families (Catalan, zeta(5)/zeta(7), gamma)")
    print("=" * 68)

    mp.dps = 21000
    G = catalan
    z5 = zeta(5)
    z7 = zeta(7)
    p = pi

    print("\nCatalan's constant")
    ratio = G / p ** 2
    r = lll_relation([ratio ** k for k in range(31)], scale_digits=20000)
    report("G/pi^2 algebraic, degree <= 30 (n=31)", r,
           claimed="mirrors the zeta(3)/pi^3 algebraicity test - not previously tested for G")

    r = lll_relation(bivariate_basis(G, p, 8), scale_digits=20000)
    report("bivariate G^i * pi^j, total degree <= 8 (n=45)", r,
           claimed="mirrors the zeta(3),pi bivariate test - not previously tested for G")

    print("\nOdd zeta values beyond zeta(3): algebraicity and joint polynomials")
    ratio5 = z5 / p ** 5
    r = lll_relation([ratio5 ** k for k in range(31)], scale_digits=20000)
    report("zeta(5)/pi^5 algebraic, degree <= 30 (n=31)", r,
           claimed="only linear independence of zeta(5) from pi was tested before")

    ratio7 = z7 / p ** 7
    r = lll_relation([ratio7 ** k for k in range(31)], scale_digits=20000)
    report("zeta(7)/pi^7 algebraic, degree <= 30 (n=31)", r,
           claimed="only linear independence of zeta(7) from pi was tested before")

    mp.dps = 5000
    z5_ = zeta(5)
    z7_ = zeta(7)
    r = lll_relation(bivariate_basis(z5_, z7_, 5), scale_digits=4000)
    report("bivariate zeta(5)^i * zeta(7)^j, total degree <= 5 (n=21)", r,
           claimed="only linear independence was tested before; this is a genuine joint polynomial search")

    print("\nEuler-Mascheroni constant (does not appear elsewhere in this paper;"
          " unlike zeta(3), gamma's irrationality itself is unproven)")
    gam = euler
    r = lll_relation([gam, p, e, mpf(1)], scale_digits=4000)
    report("gamma vs pi, e (linear)", r)

    r = lll_relation([gam ** k for k in range(6)], scale_digits=4000)
    report("gamma algebraic, degree <= 5 (n=6)", r)

    r = lll_relation(trivariate_basis(gam, p, e, 6), scale_digits=4000)
    report("trivariate gamma, pi, e joint polynomial, total degree <= 6 (n=84)", r)

    print("\n" + "=" * 68)
    print(f"Total: {time.time()-total_t0:.1f}s")


if __name__ == "__main__":
    main()
