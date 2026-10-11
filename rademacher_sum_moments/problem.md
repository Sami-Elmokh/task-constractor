# Moments and tails of Rademacher sums

## Scientific objective

Let $\varepsilon_1, \dots, \varepsilon_n$ be independent random signs, $P(\varepsilon_i = 1) = P(\varepsilon_i = -1) = \tfrac12$, and let $c_1, \dots, c_n$ be integers. Study the random sum

$$S = \sum_{i=1}^{n} c_i \varepsilon_i$$

exactly: its law, its even moments, its tail, and the ratio of its $L^p$ and $L^2$ norms that appears in the Khintchine inequalities

$$\alpha_p \lVert S\rVert_2 \le \lVert S\rVert_p \le \beta_p \lVert S\rVert_2, \qquad \lVert S\rVert_p = \bigl(E\lvert S\rvert^p\bigr)^{1/p}.$$

## Formulas and assumptions

**Coefficients.** The $c_i$ are integers of any sign; some may be zero, and values may repeat.

**Law.** $P(S = s) = 2^{-n} N(s)$, where $N(s)$ is the number of sign vectors $(\varepsilon_1, \dots, \varepsilon_n) \in \{-1, 1\}^n$ such that $\sum_i c_i \varepsilon_i = s$.

**Tail.** $P(\lvert S\rvert > t)$, with a strict inequality, for any rational $t$.

**Ratio.** $R_p(c) = \lVert S\rVert_p / \lVert S\rVert_2$ for real $p \ge 1$. Note that $\lVert S\rVert_2^2 = \sum_i c_i^2$.

**Requirements.** Laws, moments and tail probabilities must be exact. Ratios must have a relative accuracy of $10^{-12}$. Each call must finish within 10 seconds for inputs within the limits stated below.

| Symbol or notation | Meaning | Units, if any |
| --- | --- | --- |
| $n$ | number of coefficients | — |
| $c_i$ | integer coefficients | — |
| $\varepsilon_i$ | independent random signs | — |
| $S$ | Rademacher sum $\sum_i c_i \varepsilon_i$ | — |
| $m$ | order of the even moment $E[S^{2m}]$ | — |
| $p$ | exponent of the norm $\lVert S\rVert_p$ | — |

## Required functions

### `rademacher_law(c)`

- Input: `c`, 1D `numpy` integer array of length $n \le 300$, $\lvert c_i\rvert \le 10^{12}$, such that $S$ takes at most $2 \times 10^4$ distinct values.
- Output: `dict` mapping each value $s$ with $P(S = s) > 0$ (Python `int`) to $P(S = s)$ (`fractions.Fraction`), keys in increasing order.
- Purpose: Exact law of $S$.

### `even_moment(c, m)`

- Input: `c`, 1D `numpy` integer array of length $n \le 10^5$, $\lvert c_i\rvert \le 10^6$; `m`, integer, $1 \le m \le 20$.
- Output: Python `int`, the exact value of $E[S^{2m}]$.
- Purpose: Exact even moments of $S$.

### `tail_probability(c, t)`

- Input: `c` as in `rademacher_law`; `t`, an `int`, a `fractions.Fraction` or a `float`.
- Output: `fractions.Fraction`, the exact value of $P(\lvert S\rvert > t)$.
- Purpose: Exact tail of $S$.

### `khintchine_ratio(c, p)`

- Input: `c`, 1D `numpy` integer array, not identically zero; `p`, float, $p \ge 1$. When $p$ is an even integer, `c` satisfies the limits of `even_moment`; otherwise those of `rademacher_law`.
- Output: float, $R_p(c)$.
- Purpose: Use `rademacher_law` and `even_moment` to produce the Khintchine ratio.

## Allowed libraries

`numpy` and the Python standard library.
