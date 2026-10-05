"""Complete working code for the task."""

import numpy as np


def resonant_modes(mu):
    """Coefficients [a1, b1] of cos t and sin t in the trigonometric interpolant of the samples."""
    mu = np.asarray(mu, dtype=float)
    N = mu.size
    t = 2 * np.pi * np.arange(N) / N
    a1 = 2.0 / N * np.dot(mu, np.cos(t))
    b1 = 2.0 / N * np.dot(mu, np.sin(t))
    return np.array([a1, b1])


def periodic_particular_solution(mu, t):
    """Periodic solution y_p of y'' + y = mu with the k = +-1 modes of mu removed.

    Returns an array of shape (2, len(t)): y_p(t) and y_p'(t).
    """
    mu = np.asarray(mu, dtype=float)
    t = np.atleast_1d(np.asarray(t, dtype=float))
    N = mu.size
    c = np.fft.fft(mu) / N
    # integer wave numbers in FFT order (0, 1, ..., M, -M, ..., -1); np.fft.fftfreq returns
    # floats that are not always exactly 1 for the first harmonic
    k = np.fft.ifftshift(np.arange(-(N // 2), N // 2 + 1))
    nonres = np.abs(k) != 1
    y_hat = np.zeros_like(c)
    y_hat[nonres] = c[nonres] / (1.0 - k[nonres] ** 2)
    phase = np.exp(1j * np.outer(t, k))
    y_p = np.real(phase @ y_hat)
    dy_p = np.real(phase @ (1j * k * y_hat))
    return np.array([y_p, dy_p])


def driven_response(mu, t):
    """Use the subproblem functions to produce the exact solution of y'' + y = mu, y(0) = y'(0) = 0."""
    t = np.atleast_1d(np.asarray(t, dtype=float))
    a1, b1 = resonant_modes(mu)
    y_p, _ = periodic_particular_solution(mu, t)
    y0, dy0 = periodic_particular_solution(mu, [0.0])[:, 0]
    # secular response to the resonant forcing a1 cos t + b1 sin t; its slope at t = 0 is -b1/2
    y_s = 0.5 * t * (a1 * np.sin(t) - b1 * np.cos(t))
    # free oscillation that enforces y(0) = 0 and y'(0) = 0
    return y_p + y_s - y0 * np.cos(t) - (dy0 - 0.5 * b1) * np.sin(t)
