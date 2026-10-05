# Period–energy curve of a one-dimensional potential well

## Scientific objective

A particle of unit mass moves on a line in a smooth potential $f(x)$:

$$x''(t) = -f'(x(t)).$$

Its energy $h = \tfrac12 x'^2 + f(x)$ is conserved, and a particle trapped between two turning points oscillates periodically. Compute the **period as a function of energy** across a potential well, from the bottom of the well up to energies just below the escape threshold.

## Formulas and assumptions

**Turning points.** For an energy $h$ and a position $x_0$ with $f(x_0) < h$, the turning points $a < x_0 < b$ are the ends of the maximal open interval containing $x_0$ on which $f(x) < h$. Assume that $a$ and $b$ are finite, that $f(a) = f(b) = h$, and that $f'(a) \neq 0$ and $f'(b) \neq 0$.

**Period.**

$$T(h) = 2\int_a^b \frac{du}{\sqrt{2\,\bigl(h - f(u)\bigr)}}.$$

The required relative accuracy is $10^{-8}$ for every energy considered, including $h$ within $10^{-6}$ of a barrier top and $h - f(x_0)$ as small as $10^{-8}$.

**Barriers and escape energy.** Let $x_0$ be a strict local minimum of $f$ ($f'(x_0) = 0$, $f''(x_0) > 0$). The barriers of the well are the nearest local maxima $e_{\text{left}} < x_0 < e_{\text{right}}$ of $f$, and the escape energy is

$$h^* = \min\bigl(f(e_{\text{left}}),\, f(e_{\text{right}})\bigr).$$

Both barriers are assumed to exist, and all equilibria involved are non-degenerate ($f'' \neq 0$).

**Period–energy curve.** For fractions $q \in (0, 1)$ of the well depth, the energies are

$$h(q) = f(x_0) + q\,\bigl(h^* - f(x_0)\bigr).$$

| Symbol or notation | Meaning | Units, if any |
| --- | --- | --- |
| $x$ | position | dimensionless |
| $f(x)$, $f'(x)$ | potential energy per unit mass and its derivative | dimensionless |
| $h$ | energy of the orbit | same as $f$ |
| $x_0$ | a point of the orbit (turning points, period) or the bottom of the well (barriers, curve) | same as $x$ |
| $a, b$ | turning points | same as $x$ |
| $T$ | period | dimensionless time |
| $e_{\text{left}}, e_{\text{right}}$ | barriers of the well | same as $x$ |
| $h^*$ | escape energy | same as $f$ |
| $q$ | fraction of the well depth | — |

The potential is given as two Python callables, `f(x)` and its derivative `fp(x)`, each taking and returning a float.

## Required functions

### `turning_points(f, fp, h, x0)`

- Input: `f`, `fp` callables; `h` float; `x0` float with `f(x0) < h`.
- Output: 1D `numpy` array of 2 floats, `[a, b]`, accurate to about $10^{-12}$.
- Purpose: Find the turning points of the orbit of energy $h$ through $x_0$.

### `barrier_energy(f, fp, x0)`

- Input: `f`, `fp` callables; `x0` float, a strict local minimum of $f$.
- Output: 1D `numpy` array of 3 floats, `[e_left, e_right, h_star]`.
- Purpose: Find the barriers of the well and its escape energy.

### `oscillation_period(f, fp, h, x0)`

- Input: same as `turning_points`.
- Output: float, the period $T(h)$.
- Purpose: Use `turning_points` to compute the period.

### `period_energy_curve(f, fp, x0, q)`

- Input: `f`, `fp` callables; `x0` float, a strict local minimum of $f$; `q`, 1D array-like of floats in $(0, 1)$.
- Output: 1D `numpy` array of floats of length `len(q)`, the periods $T(h(q))$.
- Purpose: Use `barrier_energy` and `oscillation_period` to produce the period–energy curve of the well.

## Allowed libraries

`numpy`, `scipy`.
