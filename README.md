# High-Precision Computational Tests on ζ(3) and π

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

**Authors:** Keith Adler  
**Date:** May 2026

---

## What is this?

An open-source research project testing whether ζ(3) (Apéry's constant) is algebraically independent from π, using **two independent integer-relation engines**: mpmath's PSLQ and LLL lattice reduction via fpylll. Every meaningful exclusion is certified by at least one engine, and cross-validated by both wherever PSLQ can reach a useful bound.

**No relation was found within the tested bounds.** Strongest results: no relation a·ζ(3) + b·π² + c = 0 exists with coefficient norm up to **10¹⁹⁹⁹⁹** (LLL; independently confirmed to 10¹⁸⁶⁹⁵ by PSLQ), and ζ(3)/π³ is not algebraic of degree ≤ 100 with coefficient norm up to **10¹⁸³**.

## Main Results

| Result | Engine(s) | Bound |
|--------|-----------|-------|
| **No relation a·ζ(3) + b·π² + c = 0** | LLL (0.2s) + PSLQ (2.9h) | **10¹⁹⁹⁹⁹** / 10¹⁸⁶⁹⁵ |
| **ζ(3)/π³ not algebraic degree ≤ 100** | LLL (355s) | **10¹⁸³** |
| ζ(3)/π³ not algebraic degree ≤ 30 | LLL (7.6s) | 10⁶⁴⁰ |
| ζ(3)/π³ not algebraic degree ≤ 25 | LLL (4.4s) | 10⁷⁶⁵ |
| ζ(3)/π³ not algebraic degree ≤ 10 / ≤ 15 | LLL + PSLQ | 10¹⁸¹⁶ / 10¹²⁴⁷ |
| **ζ(3) and π: no joint polynomial degree ≤ 12** | LLL (250s) | **10²⁰⁶** |
| ζ(3) and π: no joint polynomial degree ≤ 6 | LLL (5.5s) | 10⁷¹⁰ |
| ζ(3), ζ(5), ζ(7), ζ(9) linearly independent | PSLQ | 10⁸-10¹² |
| ζ(3) independent from {π, e^π, Γ(1/4)} | PSLQ | 10¹² |
| No MZV relation: ζ(3), ζ(3,2), ζ(2,3), π⁵ | PSLQ | 10¹⁰ |
| Known Li₃(1/2) identity recovered | LLL + PSLQ | validation |

LLL bounds are on the Euclidean norm ‖a‖₂ of the coefficient vector (the max-coefficient bound is smaller by at most √n, under one order of magnitude). LLL certificates rest on fplll's guarantee of LLL-reducedness plus the proven approximation factor; PSLQ certificates on mpmath's norm-bound termination. See [`paper.md`](paper.md) §2 for both. Note the exclusion bound shrinks as degree grows (fixed precision spread across more basis dimensions) - the degree-100 bound is smaller than the degree-30 bound, but the degree reach is over 3× larger.

**🔬 Extended sweep (§3.9f, [`lll_extended.py`](lll_extended.py), ~11 min).** Beyond the original test families, we searched five previously-untested constant families - Li₄(1/2) (whose closed form, if any, is a known open question in the literature), five more Li₃(p/q) points, new Dirichlet L-values L(χ₋₃,2/3), Γ(1/3) and Γ(1/5)-based constants, and a 10-constant "kitchen sink" joint test - and pushed the degree-30/degree-6 tests to degree 100/degree 12 to check whether the original stopping points reflected a real obstruction or just the original 2005-era PSLQ paper's practical limits (the latter: they don't). All 12 additional tests found no relation, at LLL-certified bounds from 10³⁹⁸ to 10⁹⁹⁹.

**📜 Correction history (kept for transparency).** An earlier version of this repository claimed the degree-25 (10²⁰⁰), degree-30 (10¹⁰⁰), and bivariate degree-6 (10⁵⁰) bounds from PSLQ runs that never actually reached them (mpmath's PSLQ defaults to 100 iterations - far too few). Honest PSLQ re-verification (8.7-12.3 hours per test) produced only trivial bounds (17, 2, and 0 respectively), and the claims were withdrawn - mpmath's fixed-point PSLQ simply cannot make practical progress on 26+ element bases. The LLL engine then re-established all three exclusions in seconds, at bounds hundreds of orders of magnitude beyond the original claims. Full history: [`paper.md`](paper.md) §3.9.

## Quick Start

```bash
git clone https://github.com/keithadler/zeta3.git
cd zeta3
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python lll_tests.py     # LLL suite: all headline results, <30 seconds
python run_tests.py     # PSLQ suite: ~15 minutes
python lll_extended.py  # degree-100/degree-12 push + 5 new constant families, ~11 min
```

The LLL suite ([`lll_tests.py`](lll_tests.py)) reproduces every headline bound in under 30 seconds. The PSLQ suite (~15 minutes, dominated by the main test's 10700 iterations) provides the independent cross-validation. The extended sweep ([`lll_extended.py`](lll_extended.py), ~11 min) pushes the degree-30/degree-6 tests further and searches new constant families. The extended PSLQ verification (40000 digits, 100000 iterations, ~2.9 hours) is documented separately in [`paper.md`](paper.md) and is not part of the default suite.

## Repository Structure

```
├── README.md           This file
├── LICENSE             MIT License
├── requirements.txt    Python dependencies
├── CONTRIBUTING.md     How to contribute
├── paper.md            Full research paper
├── lll_tests.py        LLL (fpylll) test suite - all headline results, <30s
├── lll_extended.py     Degree-100/degree-12 push + 5 new constant families, ~11 min
├── run_tests.py        PSLQ test suite (34 of 37 tests; 3 large-basis tests are LLL-only - see correction history above)
├── generate_figures.py Optional figure generation (requires matplotlib)
└── figures/            Generated figures for the paper
```

## How It Works

Two independent methods:

**PSLQ** ([Ferguson-Bailey](https://www.cecm.sfu.ca/organics/papers/bailey/paper/html/node3.html)) takes real numbers computed to D digits and either finds an integer relation or **certifies** none exists with coefficients below a bound M, via its internal norm bound. Caveat discovered during this project: mpmath's implementation defaults to 100 iterations, and the norm bound only grows ~0.19 digits/iteration (3-element basis) - far less for large bases - so large `maxcoeff` values require explicitly requesting enough iterations, and 26+ element bases are impractical for it entirely.

**LLL** (via [fpylll](https://github.com/fplll/fpylll)) reduces the lattice [I | 10^S·x] whose short vectors correspond to integer relations. fplll's guarantee of LLL-reducedness plus the proven approximation factor yields exclusion certificates directly from the shortest reduced vector - and it handles 30+ element bases in seconds.

## The Paper

[`paper.md`](paper.md) contains the full write-up: plain-English summary, formal results, continued fraction analysis (500 terms), digit normality test, and complete appendix with all 37 test parameters.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). We especially welcome:
- Reproductions on different hardware
- Extensions to higher degree or new constant families
- Cross-validation with other software (Sage, ARPREC, Julia)

## Citation

```bibtex
@misc{adler2026zeta3,
  title={High-Precision Computational Tests on $\zeta(3)$ and $\pi$},
  author={Adler, Keith},
  year={2026},
  note={Available at https://github.com/keithadler/zeta3}
}
```

## License

[MIT](LICENSE)
