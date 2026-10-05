# Driven oscillator: periodic response and resonance

## Scientific objective

An undamped harmonic oscillator with natural angular frequency $\omega_0 = 1$ starts from rest and is driven by a $2\pi$-periodic force $\mu(t)$:

$$y''(t) + y(t) = \mu(t), \qquad y(0) = 0, \qquad y'(0) = 0.$$

The force is not known as a formula. It is only known through $N$ equally spaced samples over one period. The goal is to compute the exact motion $y(t)$ at arbitrary times $t \ge 0$.

The physically important question is resonance. If the force has a component at the natural frequency, i.e. a $\cos t$ or $\sin t$ term, the amplitude grows without bound. If it has no such component, the motion stays bounded and is $2\pi$-periodic. The required functions must remain valid in both situations, including when the resonant component is the dominant part of the force.

## Formulas and assumptions

**Samples and interpolant.** The samples are $\mu_j = \mu(t_j)$ with

$$t_j = \frac{2\pi j}{N}, \qquad j = 0, 1, \dots, N-1, \qquad N \text{ odd}, \; N \ge 3.$$

Throughout the task, $\mu$ **is defined to be** the unique trigonometric polynomial of degree $M = (N-1)/2$ that interpolates the samples:

$$\mu(t) = \sum_{k=-M}^{M} c_k \, e^{ikt}, \qquad c_k = \frac{1}{N} \sum_{j=0}^{N-1} \mu_j \, e^{-ik t_j}.$$

Every required output is exact for this interpolant. No other discretization or approximation of $\mu$ is allowed.

**Resonant coefficients.** Writing the $k = \pm 1$ part of the interpolant as $a_1 \cos t + b_1 \sin t$ defines

$$a_1 = \frac{2}{N} \sum_{j=0}^{N-1} \mu_j \cos t_j, \qquad b_1 = \frac{2}{N} \sum_{j=0}^{N-1} \mu_j \sin t_j.$$

**Periodic particular solution.** Let $\tilde\mu(t) = \mu(t) - a_1 \cos t - b_1 \sin t$ be the forcing with its resonant part removed. The periodic particular solution $y_p$ is the unique $2\pi$-periodic solution of

$$y_p'' + y_p = \tilde\mu$$

that contains no $\cos t$ or $\sin t$ term. It generally does **not** satisfy the initial conditions.

**Full response.** $y(t)$ is the unique solution of the initial value problem above, with the full interpolant $\mu$, including its resonant part. It must be evaluated exactly at the requested times. Long times, $t$ up to about 100, must not lose accuracy.

| Symbol or notation | Meaning | Units, if any |
| --- | --- | --- |
| $t$ | time | dimensionless (natural period $2\pi$) |
| $N$ | number of samples, odd | — |
| $\mu_j$ | sampled force at $t_j$ | force per unit mass (dimensionless) |
| $c_k$ | complex Fourier coefficient of the interpolant | same as $\mu$ |
| $a_1, b_1$ | coefficients of $\cos t$ and $\sin t$ in the interpolant | same as $\mu$ |
| $y_p$ | periodic particular solution defined above | displacement (dimensionless) |
| $y$ | response starting from rest | displacement (dimensionless) |

## Required functions

### `resonant_modes(mu)`

- Input: `mu`, 1D `numpy` array of floats of odd length $N \ge 3$, the samples $\mu_j$.
- Output: 1D `numpy` array of 2 floats, `[a1, b1]`.
- Purpose: Extract the component of the force at the natural frequency.

### `periodic_particular_solution(mu, t)`

- Input: `mu` as above; `t`, 1D array-like of floats (times, any real values).
- Output: `numpy` array of floats of shape `(2, len(t))`. Row 0 is $y_p(t)$ and row 1 is $y_p'(t)$ at each requested time.
- Purpose: Compute the bounded periodic response to the non-resonant part of the force.

### `driven_response(mu, t)`

- Input: `mu` as above; `t`, 1D array-like of floats (times $t \ge 0$).
- Output: 1D `numpy` array of floats of length `len(t)`, the values $y(t)$.
- Purpose: Use `resonant_modes` and `periodic_particular_solution` to produce the exact response from rest, resonant or not.

## Allowed libraries

`numpy` only.
