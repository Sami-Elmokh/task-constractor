"""Tests for every required function. Repeat this pattern for each function."""

import numpy as np

from solution import driven_response, forcing_coefficients, mode_response

# Expected values come from 50-digit arithmetic, so they are exact to double precision. A correct
# double-precision method reaches about 1e-13 relative error even at t = 200 and 1e-12 from
# resonance; RTOL = 1e-9 leaves a margin of 1e4 and is the accuracy required in problem.md.
# ATOL covers values that are exactly zero (t = 0, sin forcing with k = 0, tiny coefficients).
RTOL = 1e-9
ATOL = 1e-11

T = [0.0, 0.7, 5.3, 31.0, 97.0, 200.0]


def test_forcing_coefficients_normal_case():
    t = 2 * np.pi * np.arange(11) / 11
    test_input = 0.5 + 2 * np.cos(t) - 0.5 * np.sin(t) + 0.3 * np.cos(4 * t)
    expected = [[0.5, 2.0, 0.0, 0.0, 0.3, 0.0],
                [0.0, -0.5, 0.0, 0.0, 0.0, 0.0]]
    actual = forcing_coefficients(test_input)
    assert np.shape(actual) == (2, 6)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_forcing_coefficients_boundary_case():
    # smallest grid, N = 3
    test_input = np.array([1.0, 2.0, 4.0])
    expected = [[2.3333333333333335, -1.3333333333333333],
                [0.0, -1.1547005383792515]]
    actual = forcing_coefficients(test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_forcing_coefficients_difficult_case():
    # exp(cos(t - 0.4)) with N = 49: a_k = 2 I_k(1) cos(0.4 k), b_k = 2 I_k(1) sin(0.4 k), a_0 = I_0(1)
    test_input = np.exp(np.cos(2 * np.pi * np.arange(49) / 49 - 0.4))
    expected = [
        [1.2660658777520084, 1.0410920121861964, 0.18915262460987786, 0.01606580135173126,
         -0.00015984520587153308, -0.00022593706718119485, -3.316599528836625e-05, -3.013638290084386e-06,
         -1.9887277872325843e-07, -9.897317934016424e-09, -3.5988938496167133e-10, -7.677041761795781e-12,
         9.092476382320844e-14, 1.8699734220304465e-14, 1.1042181324520069e-15, 4.552096374884675e-17,
         1.4700925098659199e-18, 3.7818090921282395e-20, 7.345409751149074e-22, 7.97839633896059e-24,
         -1.1543495400997879e-25, -9.803901858737783e-27, -3.478513420892967e-28, -9.084553359900155e-30,
         -1.9109315085693527e-31],
        [0.0, 0.4401666428347915, 0.19475883546560382, 0.04132367700761477, 0.005471906242619154,
         0.0004936814983793065, 3.0380525615384483e-05, 1.0714383139312388e-06, -1.1628857919035804e-08,
         -4.88399711652003e-09, -4.1668758920032325e-10, -2.377060697208949e-11, -1.035166675394809e-12,
         -3.526100194429548e-14, -8.987709326376963e-16, -1.3246882288715108e-17, 1.7251380870318525e-19,
         2.1493532988300278e-20, 9.58297537463796e-22, 3.0734902515175494e-23, 7.849243791144879e-25,
         1.613438644781043e-26, 2.508519085124766e-28, 2.0771078289568653e-30, -3.3830674062473455e-32],
    ]
    actual = forcing_coefficients(test_input)
    assert np.shape(actual) == (2, 25)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_mode_response_normal_case():
    test_input = (1, 2.5, T)
    expected = [[0.0, 0.17963585579694866, -0.0420783112285957, 0.2706897135007949,
                 -0.01872901611505657, 0.2611498949406636],
                [0.0, 0.047737773045355474, -0.20664987377356125, -0.14265945550209458,
                 0.11516081225531676, -0.13070258573047697]]
    actual = mode_response(*test_input)
    assert np.shape(actual) == (2, 6)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_mode_response_boundary_case():
    # exact resonance omega = k = 3
    test_input = (3, 3.0, T)
    expected = [[0.0, 0.10070775944236861, -0.1685917468805332, -4.899457729894728,
                 14.87356331218476, 1.4727482777291065],
                [0.0, 0.10685478812825412, 0.8564922423315181, -1.692730632366789,
                 6.386585483520831, 33.303237208226406]]
    actual = mode_response(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_mode_response_difficult_case():
    # omega = 3 + 1e-12, within 1e-12 of resonance
    test_input = (3, 3.0 + 1e-12, T)
    expected = [[0.0, 0.10070775944233121, -0.1685917468828031, -4.899457729868487,
                 14.873563311874983, 1.4727482743984868],
                [0.0, 0.10685478812823594, 0.856492242330643, -1.692730632441891,
                 6.386585484239069, 33.30323720835704]]
    actual = mode_response(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_mode_response_near_resonance_case():
    # additional regime: omega = 2 + 5e-9, close to resonance but not within 1e-8 of it
    test_input = (2, 2.0 + 5e-9, T)
    expected = [[0.0, 0.17245370258446582, -1.2226774388737154, -5.728649987345405,
                 -17.032864277664277, -42.54598106117764],
                [0.0, 0.0934369661924207, 0.3952261573884994, -5.3120785191390345,
                 -17.348853633385296, 26.15843064108955]]
    actual = mode_response(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_mode_response_zero_frequency_case():
    # additional regime: k = 0, response to a constant force
    test_input = (0, 0.8, T)
    expected = [[0.0, 0.23866388904153726, 2.2734608283489135, 0.08570232164348228,
                 2.484275517828056, 3.086920801242561],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
    actual = mode_response(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_complete_solution_normal_case():
    # non-resonant natural frequency omega = sqrt(2)
    t = 2 * np.pi * np.arange(11) / 11
    test_input = (0.5 + 2 * np.cos(t) - 0.5 * np.sin(t) + 0.3 * np.cos(4 * t), np.sqrt(2.0), T)
    expected = [0.0, 0.5504473854349862, 1.338753313736946, 0.027326172162507333,
                -3.2035439133733723, -0.5118739469951823]
    actual = driven_response(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_complete_solution_boundary_case():
    # exact resonance with the second harmonic of exp(cos(t - 0.4)), omega = 2
    test_input = (np.exp(np.cos(2 * np.pi * np.arange(49) / 49 - 0.4)), 2.0, T)
    expected = [0.0, 0.5502433350539424, 0.5466872525938946, -1.9366124086700387,
                -6.983305572215448, -2.195288233518098]
    actual = driven_response(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_complete_solution_difficult_case():
    # near resonance with the second harmonic, omega = 2 + 5e-9
    test_input = (np.exp(np.cos(2 * np.pi * np.arange(49) / 49 - 0.4)), 2.0 + 5e-9, T)
    expected = [0.0, 0.5502433345955146, 0.5466872284480815, -1.936612496461097,
                -6.983305811413602, -2.1952954045791615]
    actual = driven_response(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_complete_solution_below_resonance_case():
    # additional regime: omega = 3 - 1e-12, just below resonance with the third harmonic
    test_input = (np.exp(np.cos(2 * np.pi * np.arange(49) / 49 - 0.4)), 3.0 - 1e-12, T)
    expected = [0.0, 0.4447508749787031, 0.4598437028911477, 0.029297556376841567,
                0.625575632763717, 1.8086247748613948]
    actual = driven_response(*test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)
