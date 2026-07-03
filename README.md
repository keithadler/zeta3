# High-Precision Computational Tests on ζ(3) and π

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

**Author:** Keith Adler  
**Date:** July 2026

---

## What is this?

An open-source research project testing whether ζ(3) (Apéry's constant) is algebraically independent from π. We use the PSLQ integer relation algorithm across 37 independent tests — 34 certified null results plus 3 recovered known identities. Every null result is verified against PSLQ's norm-bound exit condition, so each bound is a rigorous certificate rather than a failed search.

**No relation was found within the tested bounds.** If any algebraic relation between ζ(3) and π exists, it must have coefficients exceeding 10²⁰⁰⁰ (linear case) or degree greater than 30 with height above 10⁸.

## Main Results

| Result | Precision | Bound |
|--------|-----------|-------|
| **No relation a·ζ(3) + b·π² + c = 0** | 10000 digits | **10²⁰⁰⁰** |
| **ζ(3)/π³ not algebraic degree ≤ 25** | 1500 digits | **10¹⁰** |
| **ζ(3)/π³ not algebraic degree ≤ 30** | 1500 digits | **10⁸** |
| ζ(3) and π: no joint polynomial degree ≤ 6 | 1500 digits | 10⁸ |
| ζ(3), ζ(5), ζ(7), ζ(9) linearly independent | 1500-3000 digits | 10⁸-10¹² |
| ζ(3) independent from {π, e^π, Γ(1/4)} | 4000 digits | 10¹² |
| No MZV relation: ζ(3), ζ(3,2), ζ(2,3), π⁵ | 4000 digits | 10¹⁰ |
| Known Li₃(1/2) identity recovered | 5000 digits | 10⁶ |

All null results are certified by PSLQ's internal norm bound - not search failures.

## Quick Start

```bash
git clone https://github.com/keithadler/zeta3.git
cd zeta3
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python run_tests.py
```

Expected runtime: under an hour on Apple M3, dominated by the degree-30 algebraicity certificate (~18 minutes).

## Repository Structure

```
├── README.md           This file
├── LICENSE             MIT License
├── requirements.txt    Python dependencies
├── CONTRIBUTING.md     How to contribute
├── paper.tex           Research paper (LaTeX, Experimental Mathematics submission)
├── references.bib      BibTeX references
├── paper.md            Markdown version of the paper
├── paper.pdf           Rendered paper
├── make_pdf.py         Build paper.pdf from paper.md
├── run_tests.py        Reproducible test suite (37 tests)
├── generate_figures.py Optional figure generation (requires matplotlib)
└── figures/            Generated figures for the paper
```

## How It Works

The [PSLQ algorithm](https://www.cecm.sfu.ca/organics/papers/bailey/paper/html/node3.html) takes real numbers computed to D digits and either finds an integer relation or **certifies** none exists with coefficients below a bound M. When PSLQ's internal norm bound exceeds `maxcoeff`, non-existence is proven - not merely undetected.

## The Paper

[`paper.tex`](paper.tex) contains the full write-up: formal results, continued fraction analysis (500 terms), digit normality test, and complete appendix with all 37 test parameters.

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
