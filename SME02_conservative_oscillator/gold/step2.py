def oscillation_period(f, fp, h, x0):
    '''Period of the closed orbit of energy h of x'' = -f'(x) that passes through x0.
    Input
    f:   callable, potential energy f(x)
    fp:  callable, its derivative f'(x)
    h:   float, energy of the orbit, with f(x0) < h
    x0:  float, a position visited by the orbit
    Output
    T:   float, period T = 2 * integral_a^b du / sqrt(2 (h - f(u)))
    '''
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

    T = 2 * (half(a) + half(b))
    return T
