"""Complete working code for the task."""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


def turning_points(f, fp, h, x0):
    """Ends [a, b] of the maximal open interval containing x0 on which f < h."""
    dx = 1e-2  # marching step
    eps = np.finfo(float).eps

    def edge(sign):
        x = x0
        while True:
            y = x + sign * dx
            if f(y) >= h:
                return brentq(lambda u: f(u) - h, x, y, xtol=1e-15, rtol=4 * eps)
            # a barrier narrower than the step can hide between x and y: look for a maximum of f
            if sign * fp(x) > 0 and sign * fp(y) < 0:
                xm = brentq(fp, x, y, xtol=1e-15, rtol=4 * eps)
                if f(xm) >= h:
                    return brentq(lambda u: f(u) - h, x, xm, xtol=1e-15, rtol=4 * eps)
            x = y

    return np.array([edge(-1), edge(1)])


def barrier_energy(f, fp, x0):
    """Nearest maxima e_left < x0 < e_right of f around the minimum x0, and h* = min(f(e_left), f(e_right))."""
    dx = 1e-2  # marching step
    eps = np.finfo(float).eps

    def nearest_max(sign):
        x = x0 + sign * dx
        while True:
            y = x + sign * dx
            if sign * fp(y) <= 0:
                return brentq(fp, x, y, xtol=1e-15, rtol=4 * eps)
            x = y

    e_left, e_right = nearest_max(-1), nearest_max(1)
    return np.array([e_left, e_right, min(f(e_left), f(e_right))])


def oscillation_period(f, fp, h, x0):
    """Period T = 2 * int_a^b du / sqrt(2 (h - f(u))) of the closed orbit of energy h through x0."""
    a, b = turning_points(f, fp, h, x0)
    # geometric breakpoints resolve the logarithmic peak that appears near a separatrix
    breaks = np.concatenate(([0.0], np.logspace(-8, 0, 17)))

    def half(end):
        # u = end - L s^2 removes the 1/sqrt singularity at the turning point u = end
        L = end - x0

        def g(s):
            if s == 0.0:
                return 2 * abs(L) / np.sqrt(2 * fp(end) * L)
            delta = L * s * s
            # h - f(u) = f(end) - f(end - delta); close to the end, Simpson's rule on f'
            # avoids the cancellation of the direct difference
            if abs(delta) < 1e-3:
                d = delta * (fp(end) + 4 * fp(end - delta / 2) + fp(end - delta)) / 6
            else:
                d = f(end) - f(end - delta)
            return 2 * abs(L) * s / np.sqrt(2 * d)

        return sum(quad(g, lo, hi, epsabs=1e-13, epsrel=1e-11, limit=200)[0]
                   for lo, hi in zip(breaks[:-1], breaks[1:]))

    return 2 * (half(a) + half(b))


def period_energy_curve(f, fp, x0, q):
    """Use the subproblem functions to produce the periods at h = f(x0) + q (h* - f(x0)) across the well."""
    h_star = barrier_energy(f, fp, x0)[2]
    f0 = f(x0)
    return np.array([oscillation_period(f, fp, f0 + qi * (h_star - f0), x0) for qi in np.atleast_1d(q)])
