# Driven oscillator near resonance

## Scientific objective

An undamped harmonic oscillator with natural angular frequency $\omega > 0$ starts from rest and is driven by a $2\pi$-periodic force $\mu(t)$:

$$y''(t) + \omega^2 y(t) = \mu(t), \qquad y(0) = 0, \qquad y'(0) = 0.$$

The force is known only through $N$ equally spaced samples over one period. Compute the motion $y(t)$ for any $\omega > 0$, including $\omega$ equal to, or extremely close to, the frequency of one of the harmonics of the force.

## Formulas and assumptions

**Samples and interpolant.** The samples are $\mu_j = \mu(t_j)$ with $t_j = 2\pi j / N$, $j = 0, \dots, N-1$, $N$ odd, $N \ge 3$. Throughout the task, $\mu$ **is defined to be** the unique trigonometric polynomial of degree $M = (N-1)/2$ that interpolates the samples, written in real form as

$$\mu(t) = a_0 + \sum_{k=1}^{M} \bigl(a_k \cos kt + b_k \sin kt\bigr), \qquad b_0 = 0.$$

**Mode responses.** For an integer $k \ge 0$, $C_k(t)$ and $S_k(t)$ are the solutions of

$$C_k'' + \omega^2 C_k = \cos kt, \qquad S_k'' + \omega^2 S_k = \sin kt,$$

with $C_k(0) = C_k'(0) = S_k(0) = S_k'(0) = 0$.

**Requirements.** $\omega$ may be exactly an integer or differ from an integer by as little as $10^{-12}$. Times range over $0 \le t \le 200$. All outputs must have a relative accuracy of $10^{-9}$, or an absolute accuracy of $10^{-11}$ for values close to zero.

| Symbol or notation | Meaning | Units, if any |
| --- | --- | --- |
| $t$ | time | dimensionless (force period $2\pi$) |
| $\omega$ | natural angular frequency of the oscillator | 1 / time |
| $N$, $M$ | number of samples, $M = (N-1)/2$ | — |
| $\mu_j$ | sampled force per unit mass | dimensionless |
| $a_k$, $b_k$ | real Fourier coefficients of the interpolant | same as $\mu$ |
| $C_k$, $S_k$ | responses from rest to $\cos kt$ and $\sin kt$ | displacement |
| $y$ | response from rest to $\mu$ | displacement |

## Required functions

### `forcing_coefficients(mu)`

- Input: `mu`, 1D `numpy` array of floats of odd length $N \ge 3$, the samples $\mu_j$.
- Output: `numpy` array of floats of shape `(2, M + 1)`; row 0 is $a_0, \dots, a_M$ and row 1 is $b_0, \dots, b_M$.
- Purpose: Decompose the force into harmonics.

### `mode_response(k, omega, t)`

- Input: `k`, integer $\ge 0$; `omega`, float $> 0$; `t`, 1D array-like of floats in $[0, 200]$.
- Output: `numpy` array of floats of shape `(2, len(t))`; row 0 is $C_k(t)$ and row 1 is $S_k(t)$.
- Purpose: Response of the oscillator to a single harmonic.

### `driven_response(mu, omega, t)`

- Input: `mu` as in `forcing_coefficients`; `omega`, float $> 0$; `t` as in `mode_response`.
- Output: 1D `numpy` array of floats of length `len(t)`, the values $y(t)$.
- Purpose: Use `forcing_coefficients` and `mode_response` to produce the response to the full force.

## Allowed libraries

`numpy` only.
