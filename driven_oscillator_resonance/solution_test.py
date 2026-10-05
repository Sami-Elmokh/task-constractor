"""Tests for every required function (expected answers verified independently, see verification.md)."""

import numpy as np

from solution import driven_response, periodic_particular_solution, resonant_modes

# The method is exact for the interpolant, so only floating-point error remains (about 1e-13);
# the independent references agree with each other to 6e-12. RTOL/ATOL leave a wide margin
# while still rejecting any modelling error, which is at least of order 1e-3.
RTOL = 1e-8
ATOL = 1e-9


def grid(N):
    return 2 * np.pi * np.arange(N) / N


# Forcings used below:
#   NORMAL   mu = 0.25 + 0.7 cos 2t - 0.4 sin 3t          (no resonant component)
#   HARD     mu = exp(cos(t - 0.4)), N = 49               (smooth, all harmonics, resonant)
def normal_mu(N):
    t = grid(N)
    return 0.25 + 0.7 * np.cos(2 * t) - 0.4 * np.sin(3 * t)


def hard_mu():
    return np.exp(np.cos(grid(49) - 0.4))


# ----------------------------------------------------------------------------- resonant_modes
def test_resonant_modes_normal_case():
    t = grid(21)
    mu = 2 * np.cos(t) - 0.5 * np.sin(t) + np.cos(4 * t)
    assert np.allclose(resonant_modes(mu), [2.0, -0.5], rtol=RTOL, atol=ATOL)


def test_resonant_modes_boundary_case():
    # no component at the natural frequency: both coefficients vanish
    t = grid(15)
    mu = 3.0 + np.cos(2 * t)
    assert np.allclose(resonant_modes(mu), [0.0, 0.0], rtol=RTOL, atol=ATOL)


def test_resonant_modes_difficult_case():
    # a1 = 2 I1(1) cos 0.4, b1 = 2 I1(1) sin 0.4 (modified Bessel function I1)
    expected = [1.0410920121861964, 0.4401666428347915]
    assert np.allclose(resonant_modes(hard_mu()), expected, rtol=RTOL, atol=ATOL)


# -------------------------------------------------------------- periodic_particular_solution
T2 = np.array([0.0, 0.9, 2.5, 4.4, 6.0])


def test_periodic_particular_solution_normal_case():
    expected = [
        [0.01666666666666669, 0.3243828161067452, 0.23071215556398414, 0.4688587123497476, 0.015551380623834703],
        [0.15, 0.3188514064405985, -0.39550269716754405, 0.3938439503288896, -0.1513531888302576],
    ]
    actual = periodic_particular_solution(normal_mu(15), T2)
    assert np.shape(actual) == (2, len(T2))
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_periodic_particular_solution_boundary_case():
    # a purely resonant forcing has no non-resonant part: y_p and y_p' vanish identically
    t = grid(9)
    actual = periodic_particular_solution(np.cos(t) + np.sin(t), T2)
    assert np.allclose(actual, np.zeros((2, len(T2))), rtol=RTOL, atol=ATOL)


def test_periodic_particular_solution_difficult_case():
    expected = [
        [1.2010278617727288, 1.2169486579759008, 1.305091421853389, 1.274896458537104, 1.250604760801054],
        [-0.1469029901028235, 0.17028429046022228, -0.15632428295242937, 0.16982567552723749, -0.19252801857582158],
    ]
    assert np.allclose(periodic_particular_solution(hard_mu(), T2), expected, rtol=RTOL, atol=ATOL)


# ------------------------------------------------------------------------- driven_response
T3 = np.array([0.0, 0.6, 1.9, 3.7, 8.2, 15.0])


def test_driven_response_normal_case():
    # non-resonant forcing: bounded, 2*pi-periodic response starting from rest
    expected = [0.0, 0.11569027424160262, 0.2704680019447615, 0.19155508248722364,
                0.26880269212138347, 0.1716714604436825]
    assert np.allclose(driven_response(normal_mu(15), T3), expected, rtol=RTOL, atol=ATOL)


def test_driven_response_boundary_case():
    # pure resonance mu = 1.5 sin t: y = 0.75 (sin t - t cos t)
    expected = [0.0, 0.052080828336921306, 1.170412698546053, 1.9561004823150123,
                2.7913503125604118, 9.034204899779578]
    assert np.allclose(driven_response(1.5 * np.sin(grid(9)), T3), expected, rtol=RTOL, atol=ATOL)


def test_driven_response_difficult_case():
    # resonant forcing with all harmonics, long times; N = 49 is a grid size for which
    # np.fft.fftfreq(N, 1/N) does not return exactly 1 for the first harmonic
    t = np.array([0.0, 0.5, 2.0, 7.3, 18.0, 33.3, 60.0])
    expected = [0.0, 0.3218720759441649, 3.3189070513421286, 3.304677615392212,
                -9.378180777536853, 20.824607502251148, 5.270604892733836]
    assert np.allclose(driven_response(hard_mu(), t), expected, rtol=RTOL, atol=ATOL)
