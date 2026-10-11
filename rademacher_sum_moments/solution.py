"""Complete working code for the task."""

import math
from fractions import Fraction

import numpy as np


def rademacher_law(c):
    """Exact law of S = sum_i c_i eps_i, as a dict {value: probability}."""
    c = [abs(int(x)) for x in np.asarray(c).ravel()]
    # eps_i is symmetric, so only |c_i| matters; equal magnitudes are grouped into a binomial step
    counts = {}
    for a in c:
        counts[a] = counts.get(a, 0) + 1
    law = {0: 1}                                      # value -> number of sign patterns
    for a, k in counts.items():
        if a == 0:
            law = {v: w * 2 ** k for v, w in law.items()}
            continue
        step = [(a * (2 * j - k), math.comb(k, j)) for j in range(k + 1)]
        new = {}
        for v, w in law.items():
            for s, b in step:
                new[v + s] = new.get(v + s, 0) + w * b
        law = new
    total = 2 ** len(c)
    return {v: Fraction(w, total) for v, w in sorted(law.items())}


def even_moment(c, m):
    """Exact E[S^(2m)] as a Python int, for any number of coefficients."""
    c = np.abs(np.asarray(c, dtype=np.int64).ravel())
    vals, mult = np.unique(c[c != 0], return_counts=True)
    vals = [int(v) for v in vals]
    mult = [int(k) for k in mult]
    # E exp(x S) = prod_i cosh(c_i x); with y = x^2, log cosh(sqrt(y)) = sum_j lam_j y^j
    f = [Fraction(1, math.factorial(2 * k)) for k in range(m + 1)]
    lam = [Fraction(0)] * (m + 1)                     # log f by f * (log f)' = f'
    for k in range(1, m + 1):
        s = k * f[k] - sum(j * lam[j] * f[k - j] for j in range(1, k))
        lam[k] = s / k
    h = [Fraction(0)] * (m + 1)                       # sum_i log cosh(c_i x), coefficient of x^(2j)
    for j in range(1, m + 1):
        power_sum = sum(k * v ** (2 * j) for v, k in zip(vals, mult))
        h[j] = lam[j] * power_sum
    g = [Fraction(1)] + [Fraction(0)] * m             # g = exp(h) by k g_k = sum_j j h_j g_(k-j)
    for k in range(1, m + 1):
        g[k] = sum(j * h[j] * g[k - j] for j in range(1, k + 1)) / k
    moment = g[m] * math.factorial(2 * m)
    assert moment.denominator == 1
    return int(moment)


def tail_probability(c, t):
    """Exact P(|S| > t)."""
    return sum((p for v, p in rademacher_law(c).items() if abs(v) > t), Fraction(0))


def khintchine_ratio(c, p):
    """||S||_p / ||S||_2 for real p >= 1."""
    c = np.asarray(c, dtype=np.int64).ravel()
    sigma2 = sum(int(x) ** 2 for x in c)
    if sigma2 == 0:
        raise ValueError("all coefficients are zero")
    if float(p) == round(p) and round(p) % 2 == 0:
        q = round(p)
        # exact ratio^q, then a single float conversion and root
        return float(Fraction(even_moment(c, q // 2), sigma2 ** (q // 2))) ** (1.0 / q)
    sigma = math.sqrt(sigma2)
    total = sum(float(w) * (abs(v) / sigma) ** p for v, w in rademacher_law(c).items())
    return total ** (1.0 / p)
