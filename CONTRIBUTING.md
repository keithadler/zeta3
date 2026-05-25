# Contributing

We welcome contributions to this project.

## Running the Tests

The simplest contribution is running `python run_tests.py` on your hardware and reporting results. We're interested in:
- Timings on different processors (Intel, AMD, Apple Silicon, ARM)
- Confirmation that all null results reproduce
- Cross-validation with other arbitrary-precision libraries (ARPREC, Sage, Julia)

## Extending the Test Suite

To add a new PSLQ test, follow the existing pattern in `run_tests.py`:

```python
mp.dps = YOUR_PRECISION
rel, t = run_pslq("description of relation tested",
                  [basis_elements], maxcoeff,
                  "What question does this answer?")
all_results.append(("Short label", rel, t))
```

## Future Computational Directions (Non-Algebraic)

These directions go beyond PSLQ and may reveal structure that integer relation searches cannot detect:

1. **Continued fraction of ζ(3)/π³** — Compute it to a few thousand digits. If ζ(3)/π³ were rational, the continued fraction would terminate. Comparing its statistics (geometric mean, largest partial quotient) to those of δ = π²/8 − ζ(3) would add another data point on the "genericness" of ζ(3) relative to π.

2. **High-precision numerical search for integral identities** — Test whether ζ(3) numerically equals known or conjectured integral expressions involving π (there are many integral representations for ζ(3), some involving logs or other functions). Compare at 1000+ digits to rule out simple coincidences.

3. **Test specific known series mixing ζ(3) and π** — For example, there are formulas of the form ζ(3) = rational × π³ + a fast-converging series involving e^{2πk}. Compute the series part to very high precision and examine how clean the connection looks numerically.

4. **Randomized coefficient search** — Generate random small rational coefficients and test whether any combination of ζ(3), π, and related constants gets unusually close to zero. This is weaker than PSLQ but can probe regions of coefficient space that structured searches miss.

## Improving the Paper

The paper is in `paper.md`. We welcome:
- Corrections to mathematical statements
- Additional references to related work
- Clearer exposition
- Translations

## Reporting Issues

Please open an issue if you find:
- A test that gives different results on your machine
- A claim in the paper not backed by the code
- A numerical instability or precision issue

## If You Find a Relation

If PSLQ returns a non-None result on any null test, that would be extraordinary. Before reporting:
1. Verify at higher precision (double the digits)
2. Check the residual is genuinely small (< 10⁻¹⁰⁰)
3. Confirm the coefficients are within the stated bound
4. Open an issue immediately with full output

## Code Style

- Keep it simple - this is a research tool, not production software
- Document what each test checks and why it matters
- Include precision, bounds, and timing in all output
