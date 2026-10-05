def forcing_resonant_modes(mu):
    '''Component of a sampled 2*pi-periodic forcing at the natural frequency omega0 = 1.
    Input
    mu:    1D float array of odd length N >= 3, samples mu_j = mu(t_j) at t_j = 2*pi*j/N
    Output
    modes: 1D float array [a1, b1], the cos(t) and sin(t) coefficients of the
           trigonometric interpolant of the samples
    '''
    mu = np.asarray(mu, dtype=float)
    N = mu.size
    t = 2 * np.pi * np.arange(N) / N
    a1 = 2.0 / N * np.dot(mu, np.cos(t))
    b1 = 2.0 / N * np.dot(mu, np.sin(t))
    modes = np.array([a1, b1])
    return modes
