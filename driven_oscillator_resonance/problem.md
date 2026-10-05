# Driven oscillator: periodic response and resonance

## Scientific objective

An undamped harmonic oscillator with natural angular frequency $\omega_0 = 1$ starts from rest and is driven by a $2\pi$-periodic force $\mu(t)$:

$$y''(t) + y(t) = \mu(t), \qquad y(0) = 0, \qquad y'(0) = 0.$$

The force is known only through $N$ equally spaced samples over one period. Compute the exact motion $y(t)$ at arbitrary times $t \ge 0$. The force may or may not contain a component at the natural frequency, and the required functions must be valid in both cases.

## Formulas and assumptions

**Samples and interpolant.** The samples are $\mu_j = \mu(t_j)$ with

$$t_j = \frac{2\pi j}{N}, \qquad j = 0, 1, \dots, N-1, \qquad N \text{ odd}, \; N \ge 3.$$

Throughout the task, $\mu$ **is defined to be** the unique trigonometric polynomial of degree $M = (N-1)/2$ that interpolates the samples:

$$\mu(t) = \sum_{k=-M}^{M} c_k \, e^{ikt}, \qquad c_k = \frac{1}{N} \sum_{j=0}^{N-1} \mu_j \, e^{-ik t_j}.$$

All outputs are exact for this interpolant; no other approximation of $\mu$ is allowed.

**Resonant coefficients.** $a_1$ and $b_1$ are the real coefficients of the $k = \pm 1$ part of the interpolant, written $a_1 \cos t + b_1 \sin t$.

**Periodic particular solution.** With $\tilde\mu(t) = \mu(t) - a_1 \cos t - b_1 \sin t$, $y_p$ is the unique $2\pi$-periodic solution of $y_p'' + y_p = \tilde\mu$ that contains no $\cos t$ or $\sin t$ term.

**Full response.** $y(t)$ is the solution of the initial value problem above with the full interpolant $\mu$. Times up to $t \approx 100$ are required.

| Symbol or notation | Meaning | Units, if any |
| --- | --- | --- |
| $t$ | time | dimensionless (natural period $2\pi$) |
| $N$ | number of samples | — |
| $\mu_j$ | sampled force at $t_j$ | force per unit mass (dimensionless) |
| $c_k$ | complex Fourier coefficient of the interpolant | same as $\mu$ |
| $a_1, b_1$ | coefficients of $\cos t$ and $\sin t$ in the interpolant | same as $\mu$ |
| $y_p$ | periodic particular solution | displacement (dimensionless) |
| $y$ | response starting from rest | displacement (dimensionless) |

## Required functions

### `resonant_modes(mu)`

- Input: `mu`, 1D `numpy` array of floats of odd length $N \ge 3$, the samples $\mu_j$.
- Output: 1D `numpy` array of 2 floats, `[a1, b1]`.
- Purpose: Extract the component of the force at the natural frequency.

### `periodic_particular_solution(mu, t)`

- Input: `mu` as above; `t`, 1D array-like of floats (times).
- Output: `numpy` array of floats of shape `(2, len(t))`. Row 0 is $y_p(t)$ and row 1 is $y_p'(t)$.
- Purpose: Compute the periodic response to the non-resonant part of the force.

### `driven_response(mu, t)`

- Input: `mu` as above; `t`, 1D array-like of floats (times $t \ge 0$).
- Output: 1D `numpy` array of floats of length `len(t)`, the values $y(t)$.
- Purpose: Use `resonant_modes` and `periodic_particular_solution` to produce the response from rest.

## Allowed libraries

`numpy` only.
