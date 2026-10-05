# Period–energy curve of a one-dimensional potential well

## Scientific objective

A particle of unit mass moves on a line in a smooth potential $f(x)$:

$$x''(t) = -f'(x(t)).$$

Its energy $h = \tfrac12 x'^2 + f(x)$ is conserved. When the particle is trapped between two turning points, it oscillates periodically. The goal is to compute the **period as a function of energy** across a whole potential well, from small oscillations at the bottom of the well up to energies just below the escape threshold.

This curve is a basic characteristic of any anharmonic oscillator, for example a pendulum, a molecular bond or an atom in an optical lattice. Close to small oscillations the period tends to $2\pi/\sqrt{f''(x_0)}$. As the energy approaches the top of the lowest barrier of the well (the separatrix), the period diverges logarithmically. An accurate computation must handle the singular integrand at the turning points, very narrow barriers, orbits that pass over inner barriers, and asymmetric wells.

## Formulas and assumptions

**Turning points.** For an energy $h$ and a position $x_0$ with $f(x_0) < h$, the turning points $a < x_0 < b$ are the ends of the **maximal** open interval containing $x_0$ on which $f(x) < h$. You may assume that $a$ and $b$ are finite, that $f(a) = f(b) = h$, and that they are simple, i.e. $f'(a) \neq 0$ and $f'(b) \neq 0$.

Two situations must be handled correctly:
- the region where $f \ge h$ that separates the well from a neighbouring one can be extremely narrow, of order $10^{-3}$, when $h$ is just below the top of a barrier;
- if $h$ is above an inner barrier, the orbit passes over it and the interval extends to the next points where $f$ reaches $h$.

**Period.** The motion goes from $a$ to $b$ and back, so

$$T(h) = 2\int_a^b \frac{du}{\sqrt{2\,\bigl(h - f(u)\bigr)}}.$$

The integrand is infinite at both ends, an integrable $1/\sqrt{\cdot}$ singularity. The required relative accuracy is $10^{-8}$ in two extreme regimes as well:
- for $h$ within $10^{-6}$ of a barrier top;
- for tiny oscillations, with $h - f(x_0)$ as small as $10^{-8}$. Here $h - f(u)$ is the difference of two nearly equal numbers and must be evaluated without catastrophic cancellation.

**Barriers and escape energy.** Let $x_0$ be a strict local minimum of $f$, i.e. $f'(x_0) = 0$ and $f''(x_0) > 0$. The barriers of the well are the nearest local maxima $e_{\text{left}} < x_0 < e_{\text{right}}$ of $f$, which are unstable equilibria where $f' = 0$ and $f'' < 0$. The escape energy is

$$h^* = \min\bigl(f(e_{\text{left}}),\, f(e_{\text{right}})\bigr).$$

The potential need not be symmetric. Both barriers are assumed to exist, and all equilibria involved are non-degenerate ($f'' \neq 0$).

**Period–energy curve.** For a list of fractions $q \in (0, 1)$, the energies across the well are

$$h(q) = f(x_0) + q\,\bigl(h^* - f(x_0)\bigr),$$

so $q \to 0$ is the small-oscillation limit and $q \to 1$ is the escape threshold.

| Symbol or notation | Meaning | Units, if any |
| --- | --- | --- |
| $x$ | position | dimensionless |
| $f(x)$, $f'(x)$ | potential energy per unit mass and its derivative (minus the force) | dimensionless |
| $h$ | energy of the orbit | same as $f$ |
| $x_0$ | a point of the orbit (turning points, period) or the bottom of the well (barriers, curve) | same as $x$ |
| $a, b$ | turning points, $a < x_0 < b$ | same as $x$ |
| $T$ | period of the oscillation | dimensionless time |
| $e_{\text{left}}, e_{\text{right}}$ | nearest unstable equilibria around $x_0$ | same as $x$ |
| $h^*$ | escape energy of the well | same as $f$ |
| $q$ | fraction of the well depth, $0 < q < 1$ | — |

The potential is given as two Python callables: `f(x)` and its derivative `fp(x)`. Both accept a float and return a float. They are smooth, and vectorized evaluation is not guaranteed.

## Required functions

### `turning_points(f, fp, h, x0)`

- Input: `f`, `fp` callables; `h` float; `x0` float with `f(x0) < h`.
- Output: 1D `numpy` array of 2 floats, `[a, b]`, accurate to about $10^{-12}$.
- Purpose: Find the turning points of the orbit of energy $h$ through $x_0$.

### `barrier_energy(f, fp, x0)`

- Input: `f`, `fp` callables; `x0` float, a strict local minimum of $f$.
- Output: 1D `numpy` array of 3 floats, `[e_left, e_right, h_star]`.
- Purpose: Find the barriers that bound the well and its escape energy.

### `oscillation_period(f, fp, h, x0)`

- Input: same as `turning_points`.
- Output: float, the period $T(h)$ with relative accuracy $10^{-8}$.
- Purpose: Use `turning_points` and compute the period integral.

### `period_energy_curve(f, fp, x0, q)`

- Input: `f`, `fp` callables; `x0` float, a strict local minimum of $f$; `q`, 1D array-like of floats in $(0, 1)$.
- Output: 1D `numpy` array of floats of length `len(q)`, the periods $T(h(q))$.
- Purpose: Use `barrier_energy` and `oscillation_period` to produce the period across the well, up to the escape threshold.

## Allowed libraries

`numpy`, `scipy`.
