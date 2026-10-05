"""Tests for every required function. Repeat this pattern for each function."""

import numpy as np

from solution import barrier_energy, oscillation_period, period_energy_curve, turning_points

# Turning points and barriers are exact to about 1e-13; reference periods are exact (elliptic
# integrals) or come from ODE integration accurate to about 1e-9 near escape. RTOL = 1e-7 leaves
# a margin of 100 over the least accurate reference and is 10 times looser than the accuracy
# required in problem.md; every modelling error changes the result by at least 1e-3.
RTOL = 1e-7
ATOL = 1e-9


def test_turning_points_normal_case():
    # harmonic well 2 x^2 at energy 3: a, b = -+sqrt(3/2), with an off-center starting point
    test_input = (lambda x: 2.0 * x**2, lambda x: 4.0 * x, 3.0, 0.1)
    expected = [-1.224744871391589, 1.224744871391589]
    actual = turning_points(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_turning_points_boundary_case():
    # pendulum just below the barrier top: the region f >= h around pi is only 2.8e-3 wide
    test_input = (lambda x: -np.cos(x), lambda x: np.sin(x), 1 - 1e-6, 0.0)
    expected = [-3.1401784399095485, 3.1401784399095485]
    actual = turning_points(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_turning_points_difficult_case():
    # trapped lattice 0.01 x^4 - cos 2x at h = 1.5, above its inner barriers (height 1.064 at
    # x = -+1.613): the orbit spans three wells
    test_input = (lambda x: 0.01 * x**4 - np.cos(2 * x), lambda x: 0.04 * x**3 + 2 * np.sin(2 * x), 1.5, 0.0)
    expected = [-3.7167166245096537, 3.7167166245096537]
    actual = turning_points(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_turning_points_asymmetric_well_case():
    # additional regime: tilted potential -cos x + 0.2 x, the well is not symmetric about its bottom
    test_input = (lambda x: -np.cos(x) + 0.2 * x, lambda x: np.sin(x) + 0.2, 0.3, -np.arcsin(0.2))
    expected = [-2.497138225795214, 1.5589967108207259]
    actual = turning_points(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_barrier_energy_normal_case():
    test_input = (lambda x: -np.cos(x), lambda x: np.sin(x), 0.0)
    expected = [-np.pi, np.pi, 1.0]
    actual = barrier_energy(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_barrier_energy_boundary_case():
    # quartic well x^2/2 - x^4/4: maxima at -+1 where the force x - x^3 vanishes exactly
    test_input = (lambda x: x**2 / 2 - x**4 / 4, lambda x: x - x**3, 0.0)
    expected = [-1.0, 1.0, 0.25]
    actual = barrier_energy(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_barrier_energy_difficult_case():
    # tilted potential -cos x + 0.2 x: the left barrier (0.392) is much lower than the right (1.648)
    test_input = (lambda x: -np.cos(x) + 0.2 * x, lambda x: np.sin(x) + 0.2, -np.arcsin(0.2))
    expected = [-2.9402347327994622, 3.342950574380124, 0.39174895055337877]
    actual = barrier_energy(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_oscillation_period_normal_case():
    # harmonic oscillator, omega = 2: T = pi at every energy
    test_input = (lambda x: 2.0 * x**2, lambda x: 4.0 * x, 3.0, 0.1)
    expected = np.pi
    actual = oscillation_period(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_oscillation_period_boundary_case():
    # tiny pendulum oscillation, h = -1 + 1e-8: T = 4 K(5e-9), just above 2 pi; h - f(u) must be
    # evaluated without catastrophic cancellation
    test_input = (lambda x: -np.cos(x), lambda x: np.sin(x), -1 + 1e-8, 0.0)
    expected = 6.283185315033568
    actual = oscillation_period(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_oscillation_period_difficult_case():
    # pendulum 1e-6 below the separatrix: T = 4 K(1 - 5e-7)
    test_input = (lambda x: -np.cos(x), lambda x: np.sin(x), 1 - 1e-6, 0.0)
    expected = 34.56249674156111
    actual = oscillation_period(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_oscillation_period_multiwell_case():
    # additional regime: orbit of the trapped lattice 0.01 x^4 - cos 2x at h = 1.5 that passes over
    # the two inner barriers (height 1.064) and spans three wells
    test_input = (lambda x: 0.01 * x**4 - np.cos(2 * x), lambda x: 0.04 * x**3 + 2 * np.sin(2 * x), 1.5, 0.0)
    expected = 10.953799670554856
    actual = oscillation_period(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_complete_solution_normal_case():
    # quartic well: anharmonic softening, T = 4 K(m) / sqrt(1 - A^2/2)
    test_input = (lambda x: x**2 / 2 - x**4 / 4, lambda x: x - x**3, 0.0, [0.25, 0.5, 0.75])
    expected = [6.626552680946376, 7.124601584453241, 8.00861904364885]
    actual = period_energy_curve(*test_input)
    assert np.shape(actual) == (3,)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_complete_solution_boundary_case():
    # small-oscillation limit of the pendulum: T -> 2 pi / sqrt(f''(0)) = 2 pi
    test_input = (lambda x: -np.cos(x), lambda x: np.sin(x), 0.0, [1e-8])
    expected = [6.28318532288755]
    actual = period_energy_curve(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_complete_solution_difficult_case():
    # tilted potential up to 1e-6 of its (lower, left) escape energy: the period diverges
    test_input = (lambda x: -np.cos(x) + 0.2 * x, lambda x: np.sin(x) + 0.2, -np.arcsin(0.2),
                  [0.5, 0.999, 0.999999])
    expected = [7.175335869853757, 13.61145371126643, 20.589465715929695]
    actual = period_energy_curve(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)
