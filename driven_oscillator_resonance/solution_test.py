"""Tests for every required function. Repeat this pattern for each function."""

import numpy as np

from solution import driven_response, periodic_particular_solution, resonant_modes

# The method is exact for the trigonometric interpolant, so only rounding error remains
# (about 1e-13); the independent references agree with each other to 6e-12. These tolerances
# leave a wide margin for any correct method while rejecting every modelling error (>= 1e-3).
RTOL = 1e-8
ATOL = 1e-9


def test_resonant_modes_normal_case():
    t = 2 * np.pi * np.arange(21) / 21
    test_input = 2 * np.cos(t) - 0.5 * np.sin(t) + np.cos(4 * t)
    expected = [2.0, -0.5]
    actual = resonant_modes(test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_resonant_modes_boundary_case():
    # no component at the natural frequency: both coefficients vanish
    t = 2 * np.pi * np.arange(15) / 15
    test_input = 3.0 + np.cos(2 * t)
    expected = [0.0, 0.0]
    actual = resonant_modes(test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_resonant_modes_difficult_case():
    # exp(cos(t - 0.4)): a1 = 2 I1(1) cos 0.4, b1 = 2 I1(1) sin 0.4 (modified Bessel function I1)
    test_input = np.exp(np.cos(2 * np.pi * np.arange(49) / 49 - 0.4))
    expected = [1.0410920121861964, 0.4401666428347915]
    actual = resonant_modes(test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_periodic_particular_solution_normal_case():
    t = 2 * np.pi * np.arange(15) / 15
    test_input = (0.25 + 0.7 * np.cos(2 * t) - 0.4 * np.sin(3 * t), [0.0, 0.9, 2.5, 4.4, 6.0])
    expected = [
        [0.01666666666666669, 0.3243828161067452, 0.23071215556398414, 0.4688587123497476, 0.015551380623834703],
        [0.15, 0.3188514064405985, -0.39550269716754405, 0.3938439503288896, -0.1513531888302576],
    ]
    actual = periodic_particular_solution(*test_input)
    assert np.shape(actual) == (2, 5)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_periodic_particular_solution_boundary_case():
    # a purely resonant force has no non-resonant part: y_p and y_p' vanish identically
    t = 2 * np.pi * np.arange(9) / 9
    test_input = (np.cos(t) + np.sin(t), [0.0, 0.9, 2.5, 4.4, 6.0])
    expected = np.zeros((2, 5))
    actual = periodic_particular_solution(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_periodic_particular_solution_difficult_case():
    test_input = (np.exp(np.cos(2 * np.pi * np.arange(49) / 49 - 0.4)), [0.0, 0.9, 2.5, 4.4, 6.0])
    expected = [
        [1.2010278617727288, 1.2169486579759008, 1.305091421853389, 1.274896458537104, 1.250604760801054],
        [-0.1469029901028235, 0.17028429046022228, -0.15632428295242937, 0.16982567552723749, -0.19252801857582158],
    ]
    actual = periodic_particular_solution(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_complete_solution_normal_case():
    # non-resonant force: bounded, 2*pi-periodic response starting from rest
    t = 2 * np.pi * np.arange(15) / 15
    test_input = (0.25 + 0.7 * np.cos(2 * t) - 0.4 * np.sin(3 * t), [0.0, 0.6, 1.9, 3.7, 8.2, 15.0])
    expected = [0.0, 0.11569027424160262, 0.2704680019447615, 0.19155508248722364,
                0.26880269212138347, 0.1716714604436825]
    actual = driven_response(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_complete_solution_boundary_case():
    # pure resonance mu = 1.5 sin t: y = 0.75 (sin t - t cos t)
    t = 2 * np.pi * np.arange(9) / 9
    test_input = (1.5 * np.sin(t), [0.0, 0.6, 1.9, 3.7, 8.2, 15.0])
    expected = [0.0, 0.052080828336921306, 1.170412698546053, 1.9561004823150123,
                2.7913503125604118, 9.034204899779578]
    actual = driven_response(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_complete_solution_difficult_case():
    # resonant force with all harmonics, long times; for N = 49, np.fft.fftfreq(N, 1/N)
    # does not return exactly 1 for the first harmonic
    test_input = (np.exp(np.cos(2 * np.pi * np.arange(49) / 49 - 0.4)), [0.0, 0.5, 2.0, 7.3, 18.0, 33.3, 60.0])
    expected = [0.0, 0.3218720759441649, 3.3189070513421286, 3.304677615392212,
                -9.378180777536853, 20.824607502251148, 5.270604892733836]
    actual = driven_response(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)
