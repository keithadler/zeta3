"""
Generate figures for the paper.
Optional - requires matplotlib: pip install matplotlib

Usage: python generate_figures.py
Output: figures/*.png
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpmath import mp, mpf, pi, zeta, pslq, floor, ln
import os

os.makedirs('figures', exist_ok=True)

# ============================================================
# Figure 1: PSLQ Norm Growth
# ============================================================
print("Generating Figure 1: PSLQ Norm Growth...")

mp.dps = 10000
z3 = zeta(3)
pi2 = pi**2

# Run with verbose to capture norm at each step
import io, contextlib
f = io.StringIO()
with contextlib.redirect_stdout(f):
    pslq([z3, pi2, mpf(1)], maxcoeff=10**18, verbose=True)
output = f.getvalue()

# Parse norms from output
norms = []
for line in output.strip().split('\n'):
    if 'Norm:' in line and '/' in line:
        parts = line.split('Norm:')
        if len(parts) > 1:
            norm_str = parts[1].strip()
            try:
                norm_val = int(norm_str)
                if norm_val > 0:
                    norms.append(norm_val)
            except ValueError:
                pass

fig, ax = plt.subplots(figsize=(10, 6))
ax.semilogy(range(len(norms)), norms, 'b-', linewidth=2, label='PSLQ norm bound')
ax.axhline(y=1e18, color='r', linestyle='--', linewidth=2, label='maxcoeff = 10¹⁸')
ax.set_xlabel('PSLQ Iteration', fontsize=12)
ax.set_ylabel('Norm Bound', fontsize=12)
ax.set_title('PSLQ Norm Growth for {ζ(3), π², 1} at 10000 digits', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
ax.set_ylim(bottom=1)
plt.tight_layout()
plt.savefig('figures/pslq_norm_growth.png', dpi=150)
plt.close()
print("  Done: figures/pslq_norm_growth.png")

# ============================================================
# Figure 2: Continued Fraction Partial Quotients
# ============================================================
print("Generating Figure 2: Continued Fraction Partial Quotients...")

mp.dps = 1000
z3 = zeta(3)
diff = pi**2/8 - z3

# Compute CFs
def get_cf(x, n):
    cf = []
    for _ in range(n):
        a = int(floor(x))
        cf.append(a)
        frac = x - a
        if abs(frac) < mpf(10)**(-900):
            break
        x = 1/frac
    return cf

cf_z3 = get_cf(z3, 210)[1:201]  # first 200 PQs
cf_delta = get_cf(diff, 210)[1:201]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

ax1.bar(range(len(cf_z3)), cf_z3, width=1, color='steelblue', alpha=0.7)
ax1.set_ylabel('Partial Quotient', fontsize=11)
ax1.set_title('Continued Fraction of ζ(3)', fontsize=13)
ax1.set_ylim(0, min(max(cf_z3)*1.1, 500))
ax1.axhline(y=2.69, color='red', linestyle='--', alpha=0.5, label='Khinchin K≈2.69')
ax1.legend()

ax2.bar(range(len(cf_delta)), cf_delta, width=1, color='darkorange', alpha=0.7)
ax2.set_xlabel('Position', fontsize=11)
ax2.set_ylabel('Partial Quotient', fontsize=11)
ax2.set_title('Continued Fraction of δ = π²/8 − ζ(3)', fontsize=13)
ax2.set_ylim(0, min(max(cf_delta)*1.1, 500))
ax2.axhline(y=2.69, color='red', linestyle='--', alpha=0.5, label='Khinchin K≈2.69')
ax2.legend()

plt.tight_layout()
plt.savefig('figures/continued_fractions.png', dpi=150)
plt.close()
print("  Done: figures/continued_fractions.png")

# ============================================================
# Figure 3: Coefficient Bound vs Precision
# ============================================================
print("Generating Figure 3: Coefficient Bound vs Precision...")

# Theoretical: for basis size n, max certifiable bound ≈ 10^(D/n)
precisions = np.array([1000, 2000, 3000, 4000, 5000, 8000, 10000, 12000])
basis_sizes = [3, 5, 10, 15, 26]

fig, ax = plt.subplots(figsize=(10, 6))
for n in basis_sizes:
    bounds = precisions / n
    ax.plot(precisions, bounds, 'o-', linewidth=2, markersize=6, label=f'n={n} elements')

# Mark our actual tests
actual_tests = [
    (5000, 15, 3, '●'),   # {z3, pi^2, 1}
    (8000, 12, 11, '●'),  # z3/pi^3 deg 10
    (12000, 9, 26, '●'),  # z3/pi^3 deg 25
    (5000, 12, 10, '●'),  # bivariate deg 3
]
for prec, log_bound, n, marker in actual_tests:
    ax.plot(prec, log_bound, 'r*', markersize=15, zorder=5)

ax.set_xlabel('Working Precision (digits)', fontsize=12)
ax.set_ylabel('log₁₀(max certifiable coefficient)', fontsize=12)
ax.set_title('Theoretical PSLQ Capacity: Bound ≈ 10^(D/n)', fontsize=14)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 13000)
plt.tight_layout()
plt.savefig('figures/bound_vs_precision.png', dpi=150)
plt.close()
print("  Done: figures/bound_vs_precision.png")

# ============================================================
# Figure 4: Bivariate Polynomial Test Heatmap
# ============================================================
print("Generating Figure 4: Bivariate Test Heatmap...")

# Show which (degree in ζ(3), degree in π) pairs were tested
max_deg = 7
tested = np.zeros((max_deg, max_deg))

# Mark tested regions
# Bivariate deg ≤ 3: all i+j ≤ 3
for i in range(4):
    for j in range(4):
        if i+j <= 3:
            tested[i, j] = 3  # strong (10^12)
# Bivariate deg ≤ 4: all i+j ≤ 4
for i in range(5):
    for j in range(5):
        if i+j <= 4:
            tested[i, j] = max(tested[i, j], 2)  # medium (10^8)
# Bivariate deg ≤ 6: all i+j ≤ 6
for i in range(7):
    for j in range(7):
        if i+j <= 6:
            tested[i, j] = max(tested[i, j], 1)  # weaker (10^6)

fig, ax = plt.subplots(figsize=(8, 7))
cmap = plt.cm.YlOrRd
im = ax.imshow(tested, cmap=cmap, origin='lower', aspect='equal', vmin=0, vmax=3)

ax.set_xlabel('Degree in π', fontsize=12)
ax.set_ylabel('Degree in ζ(3)', fontsize=12)
ax.set_title('Bivariate Polynomial Tests: ζ(3)ⁱ · πʲ', fontsize=14)
ax.set_xticks(range(max_deg))
ax.set_yticks(range(max_deg))

# Add text annotations
for i in range(max_deg):
    for j in range(max_deg):
        if tested[i, j] > 0:
            bounds = {1: '10⁶', 2: '10⁸', 3: '10¹²'}
            ax.text(j, i, bounds[int(tested[i, j])], ha='center', va='center', fontsize=9)

plt.colorbar(im, ax=ax, label='Strength (darker = higher bound)')
plt.tight_layout()
plt.savefig('figures/bivariate_heatmap.png', dpi=150)
plt.close()
print("  Done: figures/bivariate_heatmap.png")

print("\nAll figures generated in figures/")
