"""Tests for every required function. Repeat this pattern for each function."""

import numpy as np

from solution import doubly_stochastic_scaling, positive_diagonal_support, sinkhorn_limit

# Expected matrices come from 50-digit arithmetic. A correct double-precision method reaches about
# 1e-13 on these cases, including the nearly decomposable ones; ATOL = 1e-10 is the accuracy
# required in problem.md. Masks are compared exactly.
ATOL = 1e-10
RTOL = 0.0

# two dense 3x3 blocks coupled only through two tiny entries (1e-6 and 7e-6): total support
WEAK = np.array([[1.3, 0.7, 1.1, 0, 1e-6, 0],
                 [0.6, 1.9, 0.8, 0, 0, 0],
                 [1.2, 0.5, 1.4, 0, 0, 0],
                 [0, 0, 0, 0.9, 6.0, 1.0],
                 [0, 0, 0, 1.0, 1.1, 8.0],
                 [0, 7e-6, 0, 2.5, 0.8, 1.2]])


# --------------------------------------------------------------------- positive_diagonal_support
def test_positive_diagonal_support_normal_case():
    # every entry of a positive matrix lies on a positive diagonal
    test_input = np.array([[1., 2, 3, 4], [2, 1, 1, 3], [5, 1, 2, 1], [1, 3, 1, 2]])
    expected = np.ones((4, 4), dtype=bool)
    actual = positive_diagonal_support(test_input)
    assert np.array_equal(actual, expected)


def test_positive_diagonal_support_boundary_case():
    # upper triangular with positive diagonal: the identity is the only positive diagonal
    test_input = np.array([[1., 2, 3, 4], [0, 5, 6, 7], [0, 0, 8, 9], [0, 0, 0, 1]])
    expected = np.eye(4, dtype=bool)
    actual = positive_diagonal_support(test_input)
    assert np.array_equal(actual, expected)


def test_positive_diagonal_support_difficult_case():
    # sparse 9x9 pattern: 12 of its 24 positive entries lie on no positive diagonal
    test_input = np.array([[0, 0, 0, 0, 0, 6, 0, 2, 0], [2, 0, 6, 7, 4, 0, 0, 4, 0], [9, 1, 3, 0, 0, 0, 5, 0, 2],
                           [0, 0, 0, 0, 0, 3, 0, 0, 0], [0, 0, 0, 0, 0, 0, 6, 0, 0], [0, 6, 0, 0, 0, 3, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0, 0, 3, 5], [0, 0, 0, 7, 3, 0, 0, 6, 1], [1, 0, 3, 0, 0, 0, 9, 0, 0]],
                          dtype=float)
    expected = np.array([[0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 1, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 1, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0]],
                        dtype=bool)
    actual = positive_diagonal_support(test_input)
    assert np.array_equal(actual, expected)


# --------------------------------------------------------------------- doubly_stochastic_scaling
def test_doubly_stochastic_scaling_normal_case():
    # 2x2 closed form: diagonal entries sqrt(ad) / (sqrt(ad) + sqrt(bc)) = sqrt(8) / (sqrt(8) + sqrt(3))
    test_input = np.array([[2., 1], [3, 4]])
    expected = [[0.6202041028867288, 0.3797958971132712], [0.3797958971132712, 0.6202041028867288]]
    actual = doubly_stochastic_scaling(test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_doubly_stochastic_scaling_boundary_case():
    # already doubly stochastic, with zero entries: the scaling leaves it unchanged
    test_input = np.array([[0.5, 0.5, 0], [0.25, 0.25, 0.5], [0.25, 0.25, 0.5]])
    expected = [[0.5, 0.5, 0.0], [0.25, 0.25, 0.5], [0.25, 0.25, 0.5]]
    actual = doubly_stochastic_scaling(test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_doubly_stochastic_scaling_difficult_case():
    # nearly decomposable matrix with total support
    test_input = WEAK
    expected = [[0.42448151302227605, 0.23689682567400686, 0.33862100530955774, 0.0, 6.559941593492739e-07, 0.0],
                [0.18053478002575526, 0.5925281731641909, 0.2269370468100538, 0.0, 0.0, 0.0],
                [0.3949837069519687, 0.17057434516764286, 0.4344419478803884, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.16378624606618641, 0.7350268013462093, 0.10118695258760428],
                [0.0, 0.0, 0.0, 0.1615867710647205, 0.11965076863382508, 0.7187624603014544],
                [0.0, 6.559941593492739e-07, 0.0, 0.6746269828690931, 0.1453217740258063, 0.1800505871109413]]
    actual = doubly_stochastic_scaling(test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


# ----------------------------------------------------------------------------- sinkhorn_limit
def test_complete_solution_normal_case():
    # positive matrix: the limit is its doubly stochastic scaling
    test_input = np.array([[1., 2, 3], [4, 5, 6], [7, 8, 10]])
    expected = [[0.24602816280696546, 0.3531717489922032, 0.4008000882008313],
                [0.3687690508631263, 0.33085341018893555, 0.3003775389479381],
                [0.38520278632990823, 0.31597484081886124, 0.29882237285123053]]
    actual = sinkhorn_limit(test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_complete_solution_boundary_case():
    # the off-diagonal entry lies on no positive diagonal: the limit is the identity
    test_input = np.array([[1., 1], [0, 1]])
    expected = [[1.0, 0.0], [0.0, 1.0]]
    actual = sinkhorn_limit(test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_complete_solution_difficult_case():
    # 6 of the 17 positive entries lie on no positive diagonal; limit entries sqrt(2) - 1 and 2 - sqrt(2)
    test_input = np.array([[2., 1, 3, 0, 0], [0, 1, 1, 0, 0], [0, 4, 2, 0, 0], [1, 2, 0, 3, 1], [5, 0, 1, 1, 2]])
    expected = [[1.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.41421356237309503, 0.585786437626905, 0.0, 0.0],
                [0.0, 0.585786437626905, 0.41421356237309503, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.7101020514433644, 0.28989794855663564],
                [0.0, 0.0, 0.0, 0.28989794855663564, 0.7101020514433644]]
    actual = sinkhorn_limit(test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)


def test_complete_solution_large_removed_entry_case():
    # additional regime: the entry 1000 couples the two blocks one way only and must vanish in the limit
    test_input = WEAK.copy()
    test_input[5, 1] = 0.0
    test_input[0, 4] = 1000.0
    expected = [[0.42448168869616115, 0.23689714247270457, 0.33862116883113424, 0.0, 0.0, 0.0],
                [0.1805346784157703, 0.5925283868308381, 0.22693693475339158, 0.0, 0.0, 0.0],
                [0.39498363288806854, 0.17057447069645731, 0.43444189641547415, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.16378606211470847, 0.7350270627287487, 0.10118687515654283],
                [0.0, 0.0, 0.0, 0.16158670090854552, 0.1196508936161683, 0.7187624054752861],
                [0.0, 0.0, 0.0, 0.674627236976746, 0.145322043655083, 0.180050719368171]]
    actual = sinkhorn_limit(test_input)
    assert np.allclose(actual, expected, rtol=RTOL, atol=ATOL)
