"""
LLL-based integer relation tests via fpylll - an independent, far
faster engine than mpmath's PSLQ, used to cross-validate this paper's
results and to re-establish the large-basis exclusions that mpmath's
pure-Python PSLQ could not certify in practical time.

METHOD
------
Given reals x_1..x_n computed to D digits, build the n x (n+1) lattice

    B = [ I_n | round(N * x_i) ],   N = 10^S,  S < D (with margin)

Any integer relation sum(a_i x_i) = 0 makes the combination
sum(a_i row_i) short: its last coordinate is the accumulated rounding
error, at most 0.5*||a||_1, so the whole vector has norm about ||a||
rather than the generic det^(1/n) ~ 10^(S/n). LLL finds short vectors.

EXCLUSION CERTIFICATE
---------------------
fplll guarantees its output basis is LLL-reduced. For an LLL-reduced
basis, ||b_1|| <= 2^((n-1)/2) * lambda_1 (proven for delta = 3/4;
conservative for the delta = 0.99 used here). A relation with
coefficient vector a maps to a lattice vector of norm at most
||a|| * sqrt(1 + n/4). Therefore any relation satisfies

    ||a||_2  >=  ||b_1|| / (2^((n-1)/2) * sqrt(1 + n/4))

and we report the right-hand side (log10) as the certified exclusion
bound. Note this bounds the Euclidean norm ||a||_2; the corresponding
bound on max|a_i| is smaller by at most sqrt(n) (under one order of
magnitude for every basis here).

SECOND, GUARANTEE-FREE CERTIFICATE (exact Gram-Schmidt)
-------------------------------------------------------
The bound above rests on fplll's promise that its output is LLL-reduced
(documented for the default wrapper: (2*delta-1, 2*eta-1/2)-reduced).
A second certificate needs no such promise. For ANY basis b_1..b_n of a
lattice, lambda_1 >= min_i ||b_i*|| where b_i* are the Gram-Schmidt
vectors, and ||b_i*||^2 = det(G_i) / det(G_{i-1}) for the leading
principal minors G_i of the Gram matrix B B^T. Those determinants are
exact integers (computed with FLINT via python-flint), so

    ||a||_2  >=  min_i ||b_i*|| / sqrt(1 + n/4)

holds regardless of how the basis was produced. It is also tighter (no
2^((n-1)/2) loss): degree 100 improves from 10^183 to ~10^196. Both
bounds are printed; the exact-GSO one is the headline.

SPURIOUS-VECTOR CAUTION
-----------------------
When no relation exists, LLL still returns balanced vectors with all
coordinates ~10^(S/n) whose residual sum(a_i x_i) ~ |last|/N is small
automatically. The distinguishing test for a genuine relation: its
residual is zero to the full precision of the inputs (bounded by
||a||_1 * 10^-D), while a spurious vector's residual is pinned near
|last|/N >> 10^-D. We therefore verify every candidate at input
precision D, which must exceed S by a margin (we use 1000 digits).

Requirements: pip install mpmath gmpy2 fpylll cysignals python-flint
Usage: python lll_tests.py           (~2 minutes on Apple M3)
"""

import sys
sys.set_int_max_str_digits(0)
import math
import time

from fpylll import IntegerMatrix, LLL
from flint import fmpz_mat
from mpmath import mp, mpf, pi, zeta, polylog, ln, nstr, log10 as mp_log10


def _log10_int(q):
    """log10 of a (possibly huge) positive Python int, as a float."""
    return float(mp_log10(mpf(q)))


def gso_lower_bound_log10(rows):
    """log10 of min_i ||b_i*|| for the integer basis `rows` (list of equal-
    length integer lists), via exact Gram minors. Lower-bounds lambda_1 of
    the lattice spanned by the rows, with NO reducedness assumption.
    Returns None if the rows are linearly dependent (a zero GS vector)."""
    n = len(rows)
    B = fmpz_mat(rows)
    G = B * B.transpose()
    prev = 1
    best = None
    for i in range(1, n + 1):
        d = int(fmpz_mat([[G[r, c] for c in range(i)] for r in range(i)]).det())
        if d <= 0:
            return None
        val = 0.5 * (_log10_int(d) - _log10_int(prev))
        best = val if best is None else min(best, val)
        prev = d
    return best


def lll_relation(basis_vals, scale_digits):
    """Run one LLL integer-relation search. Inputs must already be
    computed at mp.dps >= scale_digits + margin.

    The margin is not optional. The spurious-vector filter below
    rejects a balanced no-relation vector because its residual
    (~10^(c - S) for coefficient digits c and scale S) exceeds the
    genuine-relation threshold (~10^(c - dps + 50)); that requires
    mp.dps - scale_digits > 50 + log10(n). We enforce a hard minimum
    gap of 100 digits because running with a thin margin silently
    reports spurious relations as genuine (this bit us in testing:
    a 40-digit gap returned plausible-looking 43-digit 'relations'
    on a basis with no relation)."""
    n = len(basis_vals)
    if mp.dps - scale_digits < 100:
        raise ValueError(
            f"mp.dps ({mp.dps}) must exceed scale_digits ({scale_digits}) "
            f"by at least 100; a thin precision margin makes the "
            f"spurious-relation filter unsound (see docstring)")
    N = mpf(10) ** scale_digits

    last_col = []
    for x in basis_vals:
        v = x * N
        last_col.append(int(v + mpf('0.5')) if v >= 0 else -int(-v + mpf('0.5')))

    B = IntegerMatrix(n, n + 1)
    for i in range(n):
        B[i, i] = 1
        B[i, n] = last_col[i]

    t0 = time.time()
    LLL.reduction(B)
    elapsed = time.time() - t0

    def row_norm_sq(i):
        return sum(B[i, j] ** 2 for j in range(n + 1))

    norms = sorted((row_norm_sq(i), i) for i in range(n))
    min_norm_sq = norms[0][0]

    slack_log10 = ((n - 1) / 2) * math.log10(2) + 0.5 * math.log10(1 + n / 4)
    exclusion_log10 = (0.5 * _log10_int(min_norm_sq) if min_norm_sq > 0 else 0.0) - slack_log10

    # Second certificate: exact Gram-Schmidt lower bound on lambda_1,
    # independent of any LLL-reducedness guarantee (see module docstring).
    t1 = time.time()
    rows = [[int(B[i, j]) for j in range(n + 1)] for i in range(n)]
    gso_log10 = gso_lower_bound_log10(rows)
    gso_exclusion_log10 = None
    if gso_log10 is not None:
        gso_exclusion_log10 = gso_log10 - 0.5 * math.log10(1 + n / 4)
        # Sanity: the GS bound can never exceed the shortest reduced row.
        assert gso_log10 <= 0.5 * _log10_int(min_norm_sq) + 1e-9
    gso_elapsed = time.time() - t1

    relation = None
    residual = None
    input_dps = mp.dps
    for norm_sq, i in norms:
        coeffs = [B[i, j] for j in range(n)]
        if all(c == 0 for c in coeffs):
            continue
        coeff_norm = math.isqrt(sum(c * c for c in coeffs)) or 1
        r = sum(mpf(c) * v for c, v in zip(coeffs, basis_vals))
        threshold = mpf(coeff_norm) * n * mpf(10) ** (-(input_dps - 50))
        if abs(r) < threshold:
            relation = coeffs
            residual = r
            break

    return {
        "relation": relation,
        "residual": residual,
        "exclusion_log10": exclusion_log10,
        "gso_exclusion_log10": gso_exclusion_log10,
        "gso_elapsed": gso_elapsed,
        "elapsed": elapsed,
        "n": n,
        "scale_digits": scale_digits,
    }


# Minimum exclusion bounds (log10 of ||a||_2, LLL-reducedness certificate)
# that this suite has established and the README/paper state. The run
# fails if any test comes in below its stated bound, so a regression in
# the arithmetic, the lattice construction, or the LLL engine cannot pass
# silently. A different fplll version may reduce to a slightly different
# basis, hence the 2-digit tolerance.
EXPECTED_MIN_LOG10 = {}
BOUND_TOLERANCE = 2.0
_shortfalls = []


def report(label, r, claimed=None, expect=None):
    if expect is not None and r["relation"] is None:
        EXPECTED_MIN_LOG10[label] = expect
        if r["exclusion_log10"] < expect - BOUND_TOLERANCE:
            _shortfalls.append((label, r["exclusion_log10"], expect))
    if r["relation"] is not None:
        print(f"  {label}: FOUND relation {r['relation']}")
        print(f"      residual: {nstr(abs(r['residual']), 5)}  ({r['elapsed']:.2f}s)")
    else:
        extra = f"  [prior claim: {claimed}]" if claimed else ""
        g = r.get("gso_exclusion_log10")
        gso = (f"; exact-GSO certificate ||a|| >= 10^{g:.1f} (+{r['gso_elapsed']:.1f}s)"
               if g is not None else "; exact-GSO certificate unavailable (dependent rows)")
        print(f"  {label}: no relation; LLL-reducedness certificate ||a|| >= 10^{r['exclusion_log10']:.1f}"
              f"  ({r['elapsed']:.2f}s){gso}{extra}")


def main():
    total_t0 = time.time()
    import platform
    import fpylll
    import flint
    import mpmath
    print("LLL (fpylll) integer relation test suite")
    print("=" * 64)
    print(f"Python {platform.python_version()} | mpmath {mpmath.__version__} | "
          f"fpylll {fpylll.__version__} | python-flint {flint.__version__} | "
          f"{platform.machine()} {platform.system()}")

    # ---- Validation: known identities must be found -------------------
    print("\nValidation (known identities):")
    mp.dps = 5000
    z3 = zeta(3)
    ln2 = ln(2)
    r = lll_relation([z3, polylog(3, mpf(1) / 2), pi ** 2 * ln2, ln2 ** 3], 4000)
    report("zeta(3)/Li3(1/2) identity", r)
    assert r["relation"] is not None and sorted(map(abs, r["relation"])) == [2, 4, 21, 24]

    r = lll_relation([z3, polylog(3, mpf(-1))], 4000)
    report("Li3(-1) identity", r)
    assert r["relation"] is not None and sorted(map(abs, r["relation"])) == [3, 4]

    r = lll_relation([z3 ** 2, zeta(6), pi ** 6, mpf(1)], 4000)
    report("zeta(6) identity", r)
    assert r["relation"] is not None and sorted(map(abs, r["relation"])) == [0, 0, 1, 945]

    # ---- Main result: cross-validation of PSLQ ------------------------
    print("\nMain result {zeta(3), pi^2, 1} (cross-validates PSLQ's 10^18695):")
    t0 = time.time()
    mp.dps = 61000
    z3 = zeta(3)
    pi2 = pi ** 2
    print(f"  (constants at 61000 digits: {time.time()-t0:.1f}s)")
    r = lll_relation([z3, pi2, mpf(1)], 60000)
    report("a*zeta(3) + b*pi^2 + c = 0", r, claimed="PSLQ certified 10^18695 in 2.9h", expect=19999)

    # ---- Algebraicity of zeta(3)/pi^3 ---------------------------------
    print("\nAlgebraicity of zeta(3)/pi^3 (degrees 10/15 verify PSLQ; 25/30 supersede withdrawn claims):")
    mp.dps = 21000
    z3 = zeta(3)
    ratio = z3 / pi ** 3
    for deg, claimed, expect in [(10, "PSLQ: 10^12", 1816), (15, "PSLQ: 10^9", 1247),
                                 (25, "withdrawn (was 10^200; real PSLQ bound 17)", 765),
                                 (30, "withdrawn (was 10^100; real PSLQ bound 2)", 640)]:
        basis = [ratio ** k for k in range(deg + 1)]
        r = lll_relation(basis, 20000)
        report(f"degree <= {deg} (n={deg+1})", r, claimed=claimed, expect=expect)

    # ---- Bivariate polynomial independence ----------------------------
    print("\nBivariate zeta(3)^i * pi^j independence:")
    for d, claimed, expect in [(3, "PSLQ: 10^12", 1998), (4, "PSLQ: 10^8", 1331),
                               (6, "withdrawn (was 10^50; real PSLQ bound 0)", 710)]:
        basis = []
        for total in range(d + 1):
            for i in range(total + 1):
                basis.append(z3 ** i * pi ** (total - i))
        r = lll_relation(basis, 20000)
        report(f"total degree <= {d} (n={len(basis)})", r, claimed=claimed, expect=expect)

    print("\n" + "=" * 64)
    print(f"Total: {time.time()-total_t0:.1f}s")
    if _shortfalls:
        print("\nFAIL: certified bounds below the stated minimums:")
        for label, got, exp in _shortfalls:
            print(f"  {label}: got 10^{got:.1f}, expected >= 10^{exp}")
        sys.exit(1)
    print(f"All {len(EXPECTED_MIN_LOG10)} stated exclusion bounds re-established.")


if __name__ == "__main__":
    main()
