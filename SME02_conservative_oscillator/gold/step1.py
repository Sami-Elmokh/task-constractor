def turning_points(f, fp, h, x0):
    '''Turning points of the orbit of energy h of x'' = -f'(x) that passes through x0.
    Input
    f:   callable, potential energy f(x)
    fp:  callable, its derivative f'(x)
    h:   float, energy of the orbit, with f(x0) < h
    x0:  float, a position visited by the orbit
    Output
    ab:  1D float array [a, b], ends of the maximal open interval containing x0 on which f < h
    '''
    dx = 1e-2
    eps = np.finfo(float).eps

    def edge(sign):
        x = x0
        while True:
            y = x + sign * dx
            if f(y) >= h:
                return brentq(lambda u: f(u) - h, x, y, xtol=1e-15, rtol=4 * eps)
            # a barrier narrower than dx can hide between x and y: look for a local maximum of f
            if sign * fp(x) > 0 and sign * fp(y) < 0:
                xm = brentq(fp, x, y, xtol=1e-15, rtol=4 * eps)
                if f(xm) >= h:
                    return brentq(lambda u: f(u) - h, x, xm, xtol=1e-15, rtol=4 * eps)
            x = y

    ab = np.array([edge(-1), edge(1)])
    return ab
