"""Tests for every required function (expected answers verified independently, see verification.md)."""

import numpy as np

from solution import barrier_energy, oscillation_period, period_energy_curve, turning_points

# Turning points and barriers are exact to about 1e-13. The reference periods are exact
# (elliptic integrals) or come from ODE integration accurate to about 1e-9 near escape.
# RTOL = 1e-7 leaves a margin of 100 over the least accurate reference, while every
# modelling error checked in verification.md changes the result by at least 1e-3.
RTOL = 1e-7
ATOL = 1e-9

# Potentials f and their derivatives fp
HARMONIC = (lambda x: 2.0 * x**2, lambda x: 4.0 * x)
PENDULUM = (lambda x: -np.cos(x), lambda x: np.sin(x))
TILTED = (lambda x: -np.cos(x) + 0.2 * x, lambda x: np.sin(x) + 0.2)
QUARTIC = (lambda x: x**2 / 2 - x**4 / 4, lambda x: x - x**3)
TRAPPED_LATTICE = (lambda x: 0.01 * x**4 - np.cos(2 * x), lambda x: 0.04 * x**3 + 2 * np.sin(2 * x))
X_TILTED = -np.arcsin(0.2)   # bottom of the central well of the tilted potential


# --------------------------------------------------------------------------- turning_points
def test_turning_points_normal_case():
    # harmonic well 2 x^2 at energy 3: a, b = -+sqrt(3/2), independent of the starting point
    expected = [-1.224744871391589, 1.224744871391589]
    assert np.allclose(turning_points(*HARMONIC, 3.0, 0.1), expected, rtol=RTOL, atol=ATOL)


def test_turning_points_boundary_case():
    # just below the barrier top of the pendulum: the region f >= h around pi is only 2.8e-3 wide
    expected = [-3.1401784399095485, 3.1401784399095485]
    assert np.allclose(turning_points(*PENDULUM, 1 - 1e-6, 0.0), expected, rtol=RTOL, atol=ATOL)


def test_turning_points_difficult_case():
    # energy 1.5 is above the inner lattice barriers (height 1.064 at x = -+1.613):
    # the orbit spans three wells of the trapped lattice
    expected = [-3.7167166245096537, 3.7167166245096537]
    assert np.allclose(turning_points(*TRAPPED_LATTICE, 1.5, 0.0), expected, rtol=RTOL, atol=ATOL)


def test_turning_points_asymmetric_case():
    # tilted potential: the well is not symmetric about its bottom
    expected = [-2.497138225795214, 1.5589967108207259]
    assert np.allclose(turning_points(*TILTED, 0.3, X_TILTED), expected, rtol=RTOL, atol=ATOL)


# --------------------------------------------------------------------------- barrier_energy
def test_barrier_energy_normal_case():
    expected = [-np.pi, np.pi, 1.0]
    assert np.allclose(barrier_energy(*PENDULUM, 0.0), expected, rtol=RTOL, atol=ATOL)


def test_barrier_energy_boundary_case():
    # quartic well x^2/2 - x^4/4: maxima at -+1 where the force x - x^3 vanishes exactly
    expected = [-1.0, 1.0, 0.25]
    assert np.allclose(barrier_energy(*QUARTIC, 0.0), expected, rtol=RTOL, atol=ATOL)


def test_barrier_energy_difficult_case():
    # tilted potential: the left barrier (0.392) is much lower than the right one (1.648)
    expected = [-2.9402347327994622, 3.342950574380124, 0.39174895055337877]
    assert np.allclose(barrier_energy(*TILTED, X_TILTED), expected, rtol=RTOL, atol=ATOL)


# ----------------------------------------------------------------------- oscillation_period
def test_oscillation_period_normal_case():
    # harmonic oscillator, omega = 2: T = pi at every energy
    assert np.allclose(oscillation_period(*HARMONIC, 3.0, 0.1), np.pi, rtol=RTOL, atol=ATOL)


def test_oscillation_period_boundary_case():
    # tiny pendulum oscillation, h = -1 + 1e-8: T = 4 K(5e-9), just above 2 pi
    expected = 6.283185315033568
    assert np.allclose(oscillation_period(*PENDULUM, -1 + 1e-8, 0.0), expected, rtol=RTOL, atol=ATOL)


def test_oscillation_period_difficult_case():
    # pendulum 1e-6 below the separatrix: T = 4 K(1 - 5e-7)
    expected = 34.56249674156111
    assert np.allclose(oscillation_period(*PENDULUM, 1 - 1e-6, 0.0), expected, rtol=RTOL, atol=ATOL)


def test_oscillation_period_multiwell_case():
    # orbit of the trapped lattice that crosses two inner barriers
    expected = 10.953799670554856
    assert np.allclose(oscillation_period(*TRAPPED_LATTICE, 1.5, 0.0), expected, rtol=RTOL, atol=ATOL)


# ---------------------------------------------------------------------- period_energy_curve
def test_period_energy_curve_normal_case():
    # quartic well: anharmonic softening, T = 4 K(m) / sqrt(1 - A^2/2)
    expected = [6.626552680946376, 7.124601584453241, 8.00861904364885]
    actual = period_energy_curve(*QUARTIC, 0.0, [0.25, 0.5, 0.75])
    assert np.shape(actual) == (3,)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_period_energy_curve_boundary_case():
    # small-oscillation limit of the pendulum: T -> 2 pi / sqrt(f''(0)) = 2 pi
    expected = [6.28318532288755]
    assert np.allclose(period_energy_curve(*PENDULUM, 0.0, [1e-8]), expected, rtol=RTOL, atol=ATOL)


def test_period_energy_curve_difficult_case():
    # tilted potential up to 1e-6 of its (lower, left) escape energy: the period diverges
    expected = [7.175335869853757, 13.61145371126643, 20.589465715929695]
    actual = period_energy_curve(*TILTED, X_TILTED, [0.5, 0.999, 0.999999])
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)
