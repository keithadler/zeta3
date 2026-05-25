# High-Precision Computational Tests on ζ(3) and π

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

**Authors:** Keith Adler, William R. Adler  
**Date:** May 2026

---

## What is this?

An open-source research project testing whether ζ(3) (Apéry's constant) is algebraically independent from π. We use the PSLQ integer relation algorithm at up to 14000-digit precision across 35 independent tests.

**No relation was found within the tested bounds.** If any algebraic relation between ζ(3) and π exists, it must have degree > 30 or coefficients > 10¹⁸.

## Main Results

| Result | Precision | Bound |
|--------|-----------|-------|
| **No relation a·ζ(3) + b·π² + c = 0** | 10000 digits | **10¹⁸** |
| **ζ(3)/π³ not algebraic degree ≤ 30** | 14000 digits | **10⁸** |
| ζ(3) and π: no joint polynomial degree ≤ 6 | 4000 digits | 10⁶ |
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

Expected runtime: ~30 minutes on Apple M3 (dominated by degree-30 algebraicity test).

## Repository Structure

```
├── README.md           This file
├── LICENSE             MIT License
├── requirements.txt    Python dependencies
├── CONTRIBUTING.md     How to contribute
├── paper.md            Full research paper
├── run_tests.py        Reproducible test suite (35 tests)
├── generate_figures.py Optional figure generation (requires matplotlib)
└── figures/            Generated figures for the paper
```

## How It Works

The [PSLQ algorithm](https://www.cecm.sfu.ca/organics/papers/bailey/paper/html/node3.html) takes real numbers computed to D digits and either finds an integer relation or **certifies** none exists with coefficients below a bound M. When PSLQ's internal norm bound exceeds `maxcoeff`, non-existence is proven - not merely undetected.

## The Paper

[`paper.md`](paper.md) contains the full write-up: plain-English summary, formal results, continued fraction analysis (500 terms), digit normality test, and complete appendix with all 35 test parameters.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). We especially welcome:
- Reproductions on different hardware
- Extensions to higher degree or new constant families
- Cross-validation with other software (Sage, ARPREC, Julia)

## Citation

```bibtex
@misc{adler2026zeta3,
  title={High-Precision Computational Tests on $\zeta(3)$ and $\pi$},
  author={Adler, Keith and Adler, William R.},
  year={2026},
  note={Available at https://github.com/keithadler/zeta3}
}
```

## License

[MIT](LICENSE)
