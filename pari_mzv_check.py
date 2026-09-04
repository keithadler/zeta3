"""
Independent cross-check of this project's multiple-zeta-value evaluators
against PARI/GP's zetamult (via cypari2).

lll_zagier_check.py and lll_weight11.py compute zeta(6,2), zeta(8,2)
(trigamma / Hurwitz-tail reductions) and zeta(5,3,3), zeta(7,2,2)
(iterated Hurwitz-tail reduction) with mpmath, summing via nsum, and
validate them against closed forms and theorem-guaranteed relations.
Those controls are internal to mpmath. This script compares the same
evaluators against a second, unrelated implementation - PARI's zetamult,
which uses a different algorithm (Akhilesh's) - digit for digit.

Convention: decreasing, zeta(s1,...,sk) = sum_{n1>...>nk>=1} prod n_i^-s_i,
which is also PARI's convention for zetamult([s1,...,sk]).

Requirements: pip install mpmath gmpy2 cypari2
Usage: python pari_mzv_check.py [digits]     (default 300; ~8 min at 300
       without gmpy2, ~1100 digits reproduces the paper's working precision
       and takes correspondingly longer)
"""

import sys
import time

from mpmath import mp, mpf, zeta, nsum, inf, nstr
import cypari2

pari = cypari2.Pari()


def Z(s, m):
    return zeta(s, m)


def mzv2(a, b):
    """zeta(a,b): identical to lll_zagier_check.py / lll_weight11.py."""
    return zeta(a) * zeta(b) - nsum(lambda m: Z(b, m) / m ** a, [1, inf])


def mzv3(a, b, c):
    """zeta(a,b,c): identical to lll_weight11.py."""
    A = nsum(lambda m: Z(b, m) / m ** a, [1, inf])
    B = nsum(lambda n: Z(c, n) / n ** b, [1, inf])
    C = nsum(lambda n: Z(c, n) * Z(a, n + 1) / n ** b, [1, inf])
    return zeta(a) * mzv2(b, c) - zeta(c) * A + zeta(a) * B - C


def pari_mzv(idx, digits):
    """PARI zetamult, extracted exactly as floor(value * 10^digits)."""
    p = pari.zetamult(list(idx), precision=int((digits + 40) * 3.33))
    n = int(pari.floor(p * pari(10) ** digits))
    return mpf(n) / mpf(10) ** digits


def main():
    digits = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    mp.dps = digits + 10
    print(f"MZV cross-check: mpmath evaluators vs PARI zetamult at {digits} digits")
    print("=" * 68)

    # PARI must itself reproduce a closed form before we trust it as referee.
    closed_42 = zeta(3) ** 2 - mpf(4) / 3 * zeta(6)
    d = abs(pari_mzv((4, 2), digits) - closed_42)
    print(f"  PARI zeta(4,2) vs closed form zeta(3)^2 - 4/3 zeta(6): {nstr(d, 3)}")
    assert d < mpf(10) ** (-(digits - 2)), "PARI referee failed its own control"

    worst = mpf(0)
    for idx, f in [((6, 2), lambda: mzv2(6, 2)),
                   ((8, 2), lambda: mzv2(8, 2)),
                   ((5, 3, 3), lambda: mzv3(5, 3, 3)),
                   ((7, 2, 2), lambda: mzv3(7, 2, 2))]:
        t0 = time.time()
        v = f()
        d = abs(v - pari_mzv(idx, digits))
        worst = max(worst, d)
        print(f"  zeta{idx}: |mpmath - PARI| = {nstr(d, 3)}   ({time.time() - t0:.0f}s)")
        sys.stdout.flush()
        assert d < mpf(10) ** (-(digits - 2)), f"zeta{idx} disagrees with PARI"

    print(f"\nAll four MZVs agree with PARI to {digits} digits (worst |diff| = {nstr(worst, 3)}).")


if __name__ == "__main__":
    main()
