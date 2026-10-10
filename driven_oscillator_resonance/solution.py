"""Complete working code for the task."""

import numpy as np


def forcing_coefficients(mu):
    """Real Fourier coefficients of the trigonometric interpolant of the samples.

    Returns an array of shape (2, M + 1): row 0 is a_0, ..., a_M and row 1 is b_0, ..., b_M,
    with mu(t) = a_0 + sum_k (a_k cos kt + b_k sin kt) and b_0 = 0.
    """
    mu = np.asarray(mu, dtype=float)
    N = mu.size
    M = (N - 1) // 2
    c = np.fft.rfft(mu)[: M + 1] / N          # c_k for k = 0..M
    a = 2 * c.real
    b = -2 * c.imag
    a[0] = c[0].real
    b[0] = 0.0
    return np.array([a, b])


def mode_response(k, omega, t):
    """Responses from rest of y'' + omega^2 y = cos(kt) and = sin(kt), shape (2, len(t)).

    Written so that no division by omega - k occurs: exact at resonance and free of
    cancellation when omega is arbitrarily close to k.
    """
    t = np.atleast_1d(np.asarray(t, dtype=float))
    d = omega - k
    s = omega + k
    # sin(d t / 2) / (d / 2) = t * sinc(d t / (2 pi)), smooth through d = 0
    half_sin_ratio = t * np.sinc(d * t / (2 * np.pi))
    # (cos kt - cos wt) / (w^2 - k^2) = 2 sin(s t/2) sin(d t/2) / (s d)
    y_cos = np.sin(s * t / 2) * half_sin_ratio / s
    # (sin kt - (k/w) sin wt) / (w^2 - k^2) = (-2 cos(s t/2) sin(d t/2) / d + sin(wt)/w) / s
    y_sin = (-np.cos(s * t / 2) * half_sin_ratio + np.sin(omega * t) / omega) / s
    return np.array([y_cos, y_sin])


def driven_response(mu, omega, t):
    """Use the subproblem functions to produce the response from rest of y'' + omega^2 y = mu."""
    t = np.atleast_1d(np.asarray(t, dtype=float))
    a, b = forcing_coefficients(mu)
    y = np.zeros_like(t)
    for k in range(a.size):
        y_cos, y_sin = mode_response(k, omega, t)
        y += a[k] * y_cos + b[k] * y_sin
    return y
