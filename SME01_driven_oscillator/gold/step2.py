def periodic_response(mu, tol):
    '''2*pi-periodic response of y'' + y = mu, y(0) = y'(0) = 0, on the sampling grid.
    Input
    mu:   1D float array of odd length N >= 3, samples of the forcing at t_j = 2*pi*j/N
    tol:  float, resonance threshold on sqrt(a1**2 + b1**2)
    Output
    phi:  None if sqrt(a1**2 + b1**2) > tol (the response is not periodic); otherwise a
          1D float array of length N with the response at t_j, computed for the forcing
          whose k = +-1 Fourier modes have been removed
    '''
    mu = np.asarray(mu, dtype=float)
    N = mu.size
    a1, b1 = forcing_resonant_modes(mu)
    if np.hypot(a1, b1) > tol:
        return None
    t = 2 * np.pi * np.arange(N) / N
    c = np.fft.fft(mu) / N
    # integer wave numbers in FFT order (fftfreq returns floats that may miss k == 1 exactly)
    k = np.fft.ifftshift(np.arange(-(N // 2), N // 2 + 1))
    nonres = np.abs(k) != 1
    y_hat = np.zeros_like(c)
    y_hat[nonres] = c[nonres] / (1.0 - k[nonres] ** 2)
    # particular periodic solution, its value and slope at t = 0
    y_p = np.real(np.fft.ifft(y_hat) * N)
    y0 = np.real(np.sum(y_hat))
    dy0 = np.real(np.sum(1j * k * y_hat))
    # add the homogeneous solution that enforces y(0) = y'(0) = 0
    phi = y_p - y0 * np.cos(t) - dy0 * np.sin(t)
    return phi
