# Verification

- **Scientific method checked using:**
  - Energy conservation gives $dt = dx/\sqrt{2(h - f(x))}$, so the period is twice the travel time from $a$ to $b$. Near a simple turning point, $h - f(u) \approx |f'(\text{end})|\,|u - \text{end}|$. The substitution $u = \text{end} - L s^2$ (with $L = \text{end} - x_0$) removes the $1/\sqrt{\cdot}$ singularity. Close to a turning point, $h - f(u)$ is evaluated with Simpson's rule on $f'$ to avoid cancellation. Geometric breakpoints resolve the logarithmic peak that appears near a separatrix.
  - The turning-point search marches away from $x_0$ and detects a barrier narrower than the step through a sign change of $f'$. Without this check, at $h = 1 - 10^{-6}$ the pendulum's barrier around $\pi$ (only $2.8 \times 10^{-3}$ wide) is stepped over, and the computed period becomes 241 instead of 34.56.
  - The solution was compared with the references below. The maximum relative difference over all 14 tests is $7 \times 10^{-10}$, reached at 1e−6 below the escape energy. It is $7 \times 10^{-12}$ for the pendulum near the separatrix and $\le 4 \times 10^{-14}$ in every other test.
  - **Deliberate errors.** The task folder was copied, one realistic error was injected into the copied `solution.py`, and `pytest` was run. All 9 errors were caught:

    | Injected error | Result |
    | --- | --- |
    | Fixed-step march without a check for narrow barriers | 2 tests fail |
    | Stops at the first barrier even when $h$ is above it | 4 tests fail |
    | Assumes a symmetric well, $a = 2x_0 - b$ | 4 tests fail |
    | Assumes symmetric barriers, $h^* = f(e_{\text{right}})$ | `test_barrier_energy_difficult_case` fails; the curve test never terminates (the orbit escapes over the lower barrier) |
    | Higher barrier instead of the lower one ($\max$ instead of $\min$) | Same as above |
    | Half period only | 7 tests fail |
    | Small-oscillation formula $2\pi/\sqrt{f''(x_0)}$ | 4 tests fail |
    | Midpoint rule in $u$, ignoring the endpoint singularity | 7 tests fail |
    | Energy fraction measured from 0 instead of $f(x_0)$ | 2 tests fail |

- **Expected answers checked using:** an independent script that never calls `solution.py`.
  - **Turning points:** closed forms $\pm\sqrt{h/2}$ for the harmonic well $2x^2$ and $\pm\arccos(-h)$ for the pendulum. For the trapped lattice $0.01x^4 - \cos 2x$ and the tilted potential $-\cos x + 0.2x$: a dense grid (spacing $10^{-5}$) locates the connected component of $\{f < h\}$ containing $x_0$, and its ends are refined by bisection. An ODE integration started at rest at $a$ reaches $b$ after half a period to within $3 \times 10^{-14}$, which confirms that no barrier lies in between. For the trapped lattice, the inner barriers are at $x = \pm 1.613$ with height $1.064 < 1.5$, so the orbit at $h = 1.5$ spans three wells.
  - **Barriers:** exact values $[-\pi, \pi, 1]$ (pendulum) and $[-1, 1, 1/4]$ for $x^2/2 - x^4/4$ (zeros of $x - x^3$). For the tilted potential, $\sin e = -0.2$ with $\cos e < 0$ gives $e = -\pi + \arcsin 0.2$ and $\pi + \arcsin 0.2$, with heights 0.3917 and 1.6484.
  - **Periods:**
    - harmonic well $2x^2$: $T = 2\pi/\omega = \pi$ at every energy;
    - pendulum: $T = 4K(m)$ with $m = (1+h)/2$, the complete elliptic integral of the first kind (`scipy.special.ellipk`);
    - quartic well $x^2/2 - x^4/4$: amplitude $A$ from $h = A^2/2 - A^4/4$, then $T = 4K(m)/\sqrt{1 - A^2/2}$ with $m = (A^2/2)/(1 - A^2/2)$. Cross-checked against ODE integration to $6 \times 10^{-15}$;
    - trapped lattice and tilted potential (no closed form): direct integration of $x'' = -f'(x)$ (scipy `solve_ivp`, DOP853, `rtol = 1e-13`, `atol = 1e-14`), started at rest at $a$, with the period given by the next time $x' = 0$ while increasing. On the pendulum this ODE period agrees with $4K(m)$ to $4 \times 10^{-9}$ at $10^{-6}$ below the separatrix, which bounds the accuracy of the ODE references near escape.

- **Numerical tolerances chosen because:** turning points and barriers are exact to about $10^{-13}$. Periods are exact (elliptic integrals) or accurate to about $10^{-9}$ (ODE near the escape threshold). `RTOL = 1e-7` and `ATOL = 1e-9` leave a margin of about 100 over the least accurate reference and are 10 times looser than the accuracy required in `problem.md` ($10^{-8}$). Every modelling error above changes the result by at least $10^{-3}$, far above the tolerance.

  **A numerically sensitive case on purpose.** A plain adaptive `scipy.integrate.quad` of $1/\sqrt{2(h - f(u))}$ directly in $u$ passes 13 of the 14 tests. It is accurate to $10^{-9}$ even $10^{-6}$ below the separatrix. It fails only `test_oscillation_period_boundary_case`, the tiny pendulum oscillation with $h - f(x_0) = 10^{-8}$, where it returns 6.2831881 instead of 6.2831853 (relative error $4 \times 10^{-7}$). There $h - f(u)$ is computed as the difference of two numbers close to $-1$ and loses about 8 significant digits. `problem.md` states this regime and the need to avoid cancellation explicitly, so the failure reflects a numerical-stability error and not an ambiguity. Any method that evaluates $h - f(u)$ stably near the bottom passes, for example the approach used here ($f(\text{end}) - f(u)$ from Simpson's rule on $f'$) or a Taylor expansion of $f$.

- **Important assumptions:**
  - The turning points are simple ($f'(a), f'(b) \neq 0$), so the orbit is a closed periodic orbit and not a separatrix of infinite period.
  - All equilibria involved are non-degenerate ($f'' \neq 0$); stability requires $f''(x_0) > 0$ strictly.
  - Both barriers of the well exist. For an energy above the lower barrier of a potential that is unbounded below, such as the tilted one, the particle escapes and no period exists.
  - `f` and `fp` are evaluated on scalar floats only.
