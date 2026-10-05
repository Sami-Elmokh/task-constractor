# Verification

- **Scientific method checked using:**
  - In Fourier space, $y'' + y = \mu$ decouples into $(1 - k^2)\,y_k = c_k$. The non-resonant modes give the periodic particular solution $y_p = \sum_{k \ne \pm 1} c_k/(1-k^2)\, e^{ikt}$. The resonant forcing $a_1 \cos t + b_1 \sin t$ has the secular solution $\tfrac{t}{2}(a_1 \sin t - b_1 \cos t)$, whose slope at $t = 0$ is $-b_1/2$. A free oscillation $A \cos t + B \sin t$ then enforces $y(0) = y'(0) = 0$.
  - The solution was compared with a direct numerical integration of the ODE (scipy `solve_ivp`, DOP853, `rtol = atol = 1e-13`) for all three forcing types. The maximum difference is $6 \times 10^{-14}$ (normal case), $4 \times 10^{-13}$ (pure resonance) and $6 \times 10^{-12}$ (hard case, up to $t = 60$).
  - **Deliberate errors.** The task folder was copied, one realistic error was injected into the copied `solution.py`, and `pytest` was run. All 8 errors were caught:

    | Injected error | Failing tests |
    | --- | --- |
    | $1/N$ instead of $2/N$ normalization | 3 |
    | Sign of $b_1$ flipped (FFT sign convention) | 4 |
    | Wave numbers from `np.fft.fftfreq(N, 1/N)` (returns $1.0000000000000002$ for $N = 49$, so the resonant mode is divided by about $10^{-16}$) | 2 |
    | Wrong sign of $y_p'$ | 4 |
    | No secular term (resonance ignored) | 2 |
    | Slope $-b_1/2$ of the secular term forgotten | 2 |
    | Condition $y'(0) = 0$ not enforced | 3 |
    | Particular solution returned without the free oscillation | 3 |

- **Expected answers checked using:** an independent script that never calls `solution.py`.
  - `resonant_modes`: exact values $[2, -0.5]$ and $[0, 0]$ for trigonometric polynomials. For $\mu = e^{\cos(t - 0.4)}$, the generating function $e^{\cos s} = I_0(1) + 2\sum_{k \ge 1} I_k(1) \cos ks$ gives $a_1 = 2 I_1(1)\cos 0.4$ and $b_1 = 2 I_1(1) \sin 0.4$ (modified Bessel functions from `scipy.special.iv`). Aliasing from $k \ge 48$ is below $10^{-60}$.
  - `periodic_particular_solution`: hand-derived $y_p = 0.25 - \tfrac{0.7}{3}\cos 2t + 0.05 \sin 3t$ for the normal case, $y_p \equiv 0$ for a purely resonant force, and the Bessel series $I_0(1) + \sum_{k \ge 2} \tfrac{2 I_k(1)}{1-k^2}\cos k(t-0.4)$ for the hard case.
  - `driven_response`: closed forms $y_p(t) - y_p(0)\cos t - y_p'(0) \sin t$ (normal case) and $0.75(\sin t - t\cos t)$ for $\mu = 1.5 \sin t$ (boundary case). For the hard case, the Bessel series plus the secular term. All three were cross-checked against the direct ODE integration above.

- **Numerical tolerances chosen because:** the method is exact for the interpolant, so the only error is floating-point rounding, about $10^{-13}$. The independent references agree to $6 \times 10^{-12}$. `RTOL = 1e-8` and `ATOL = 1e-9` leave a margin of three orders of magnitude for any correct method, such as explicit DFT sums, variation of constants or exact quadrature. Every modelling error listed above changes the result by at least $10^{-3}$.

- **Important assumptions:**
  - $N$ is odd, so the interpolant has no ambiguous Nyquist term.
  - $\mu$ is the trigonometric interpolant of the samples, so the expected answers are exact for it. For a smooth force such as $e^{\cos(t-0.4)}$ with $N = 49$, the interpolant equals the continuous force to machine precision.
  - The test cases cover three regimes: no resonance (bounded periodic response), pure resonance (linear growth from $t = 0$), and a mixed force with all harmonics up to long times ($t = 60$).
