# Verification
**Scientific method checked using:**

- **Approach.** Energy conservation gives the period as

  $$T = 2\int_a^b \frac{du}{\sqrt{2\,(h - f(u))}}.$$

- **Numerical treatment.**
  - The substitution $u = \text{end} - L s^2$ removes the $1/\sqrt{\cdot}$ singularity at each turning point.
  - Simpson's rule on $f'$ evaluates $h - f(u)$ near the turning points without cancellation.
  - Geometric breakpoints resolve the logarithmic peak near a separatrix.
- **Narrow barriers.** The turning-point search detects barriers narrower than its step through a sign change of $f'$. Without this check, the pendulum barrier at $h = 1 - 10^{-6}$ (only $2.8 \times 10^{-3}$ wide) is stepped over, and the period comes out as $241$ instead of $34.56$.
- **Agreement with references.**
  - Maximum relative difference: $7 \times 10^{-10}$, just below escape.
  - Non-extreme tests: at most $4 \times 10^{-14}$.
- **Deliberate-error check.** In a copy of the task, one realistic mistake at a time was introduced into `solution.py` and the tests were run. Each of these nine mistakes made at least one test fail:
  1. fixed-step march without a narrow-barrier check;
  2. stopping at the first barrier even when $h$ is above it;
  3. a symmetric-well assumption;
  4. symmetric barriers;
  5. $\max$ instead of $\min$ for $h^*$ (the barrier test fails, and the curve test does not terminate because the orbit escapes);
  6. half period;
  7. the small-oscillation formula;
  8. a midpoint rule ignoring the endpoint singularities;
  9. energy fractions measured from $0$ instead of $f(x_0)$.

**Frontier-model check:**

- ChatGPT, given only `problem.md`, returned a solution that fails 5 of the 14 tests:
  - its turning-point search uses a geometrically growing step that jumps over the pendulum barrier at $h = 1 - 10^{-6}$, so the orbit crosses an unstable equilibrium: turning points $\pm 2.3 \times 10^{17}$ instead of $\pm 3.1402$, and a period of $1.3 \times 10^{18}$ instead of $34.56$;
  - the same search jumps over the left barrier of the tilted potential near escape and raises an error;
  - it evaluates $h - f(u)$ directly, which loses accuracy for tiny oscillations: $T = 6.28289$ instead of $6.28319$ at $h - f(x_0) = 10^{-8}$ (relative error $5 \times 10^{-5}$), and a relative error of $10^{-6}$ on the period–energy curve at $q = 10^{-8}$.

**Expected answers checked using:**

An independent script that never calls `solution.py`.

- **Turning points.**
  - Harmonic: closed form $\pm\sqrt{h/2}$.
  - Pendulum: closed form $\pm\arccos(-h)$.
  - Trapped lattice (orbit over three wells) and tilted potential (asymmetric well): a dense grid (spacing $10^{-5}$) plus bisection, confirmed by an ODE trajectory started at rest at $a$ that reaches $b$ within $3 \times 10^{-14}$.
- **Barriers.**
  - Exact values $[-\pi,\ \pi,\ 1]$ and $[-1,\ 1,\ 1/4]$.
  - Tilted potential: $[-\pi + \arcsin 0.2,\ \pi + \arcsin 0.2]$, with heights $0.3917$ and $1.6484$.
- **Periods.**
  - $2x^2$: $T = \pi$.
  - Pendulum: $T = 4K(m)$ with $m = (1+h)/2$, using `scipy.special.ellipk`.
  - Quartic well: $T = 4K(m)/\sqrt{1 - A^2/2}$ with $m = \dfrac{A^2/2}{1 - A^2/2}$, cross-checked against ODE integration to $6 \times 10^{-15}$.
  - Tilted potential and three-well orbit of the trapped lattice (no closed form): direct integration of $x'' = -f'(x)$ with `scipy.integrate.solve_ivp` (DOP853, `rtol = 1e-13`, `atol = 1e-14`). On the pendulum, this method agrees with $4K(m)$ to $4 \times 10^{-9}$ at $10^{-6}$ below the separatrix.

**Numerical tolerances chosen because:**

- **Reference accuracy.** Turning points and barriers are exact to about $10^{-13}$. Periods are exact, or accurate to about $10^{-9}$ near escape.
- **Chosen values.** `RTOL = 1e-7` and `ATOL = 1e-9`. These:
  - leave a margin of 100 over the least accurate reference;
  - are 10 times looser than the accuracy required in `problem.md`;
  - are far tighter than any modelling error, each of which changes the result by at least $10^{-3}$.
- **Deliberately sensitive test.** A plain `scipy.integrate.quad` of $1/\sqrt{2(h - f(u))}$ in $u$ passes every other test. For the tiny pendulum oscillation $h - f(x_0) = 10^{-8}$, however, it returns $6.2831881$ instead of $6.2831853$, because $h - f(u)$ loses about 8 significant digits. `problem.md` states this regime ($h - f(x_0)$ as small as $10^{-8}$) and the required accuracy $10^{-8}$ explicitly, so the failure reflects a numerical-stability error, not an unstated requirement.

**Important assumptions:**

- **Simple turning points.** $f'(a) \neq 0$ and $f'(b) \neq 0$, so every orbit tested is a closed periodic orbit, not a separatrix.
- **Non-degenerate equilibria.** $f'' \neq 0$ at all equilibria, and stability requires $f''(x_0) > 0$ strictly.
- **Bounded orbit.** Both barriers of the well exist, and $q \in (0, 1)$, so the orbit never reaches the escape energy.
- **Scalar evaluation.** `f` and `fp` are evaluated on scalar floats only.
