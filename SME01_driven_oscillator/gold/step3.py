def driven_response(mu, t):
    '''Exact response of y'' + y = mu, y(0) = y'(0) = 0, at arbitrary times, resonant or not.
    Input
    mu:  1D float array of odd length N >= 3, samples of the forcing at t_j = 2*pi*j/N;
         the forcing is the trigonometric interpolant of these samples
    t:   1D float array of evaluation times
    Output
    y:   1D float array, y(t) at each requested time
    '''
    mu = np.asarray(mu, dtype=float)
    t = np.asarray(t, dtype=float)
    N = mu.size
    a1, b1 = forcing_resonant_modes(mu)
    c = np.fft.fft(mu) / N
    # integer wave numbers in FFT order (fftfreq returns floats that may miss k == 1 exactly)
    k = np.fft.ifftshift(np.arange(-(N // 2), N // 2 + 1))
    nonres = np.abs(k) != 1
    y_hat = np.zeros_like(c)
    y_hat[nonres] = c[nonres] / (1.0 - k[nonres] ** 2)
    # bounded part from the non-resonant modes
    y_p = np.real(np.exp(1j * np.outer(t, k)) @ y_hat)
    # secular part from the resonant mode a1*cos(t) + b1*sin(t)
    y_s = 0.5 * t * (a1 * np.sin(t) - b1 * np.cos(t))
    y0 = np.real(np.sum(y_hat))
    dy0 = np.real(np.sum(1j * k * y_hat)) - 0.5 * b1
    y = y_p + y_s - y0 * np.cos(t) - dy0 * np.sin(t)
    return y
