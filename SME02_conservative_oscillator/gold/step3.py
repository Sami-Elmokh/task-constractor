def barrier_energy(f, fp, x0):
    '''Escape energy of the potential well whose bottom is the stable equilibrium x0.
    Input
    f:   callable, potential energy f(x)
    fp:  callable, its derivative f'(x)
    x0:  float, strict local minimum of f (f'(x0) = 0, f''(x0) > 0)
    Output
    out: 1D float array [e_left, e_right, h_star]; e_left < x0 < e_right are the nearest
         local maxima of f (unstable equilibria) and h_star = min(f(e_left), f(e_right))
    '''
    dx = 1e-2
    eps = np.finfo(float).eps

    def nearest_max(sign):
        x = x0 + sign * dx
        while True:
            y = x + sign * dx
            if sign * fp(y) <= 0:
                return brentq(fp, x, y, xtol=1e-15, rtol=4 * eps)
            x = y

    e_left, e_right = nearest_max(-1), nearest_max(1)
    out = np.array([e_left, e_right, min(f(e_left), f(e_right))])
    return out
