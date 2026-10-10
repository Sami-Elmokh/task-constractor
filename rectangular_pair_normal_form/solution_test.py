"""Tests for every required function. Repeat this pattern for each function."""

import numpy as np

from solution import alternating_ranks, canonical_form, chain_counts, invertible_blocks

# All outputs are integers (ranks, counts, eigenvalues, sizes, 0/1 and lambda entries), so the
# tests use exact comparisons. Each test pair is built from known blocks (the expected answer)
# and then scrambled by random unimodular integer changes of basis with a fixed seed.


def build_pair(blocks):
    """Block-diagonal pair (A0, B0) from ('chain', L, side) and ('inv', lam, r) blocks, in the given order."""
    m = n = 0
    entries_A, entries_B = [], []
    for b in blocks:
        if b[0] == "inv":
            _, lam, r = b
            for i in range(r):
                entries_A.append((n + i, m + i, 1))
                entries_B.append((m + i, n + i, lam))
                if i + 1 < r:
                    entries_B.append((m + i + 1, n + i, 1))
            m += r
            n += r
        else:
            _, L, side = b
            idx, sides = {}, []
            for j in range(L):
                s = "m" if (j % 2 == 0) == (side == "m") else "n"
                sides.append(s)
                if s == "m":
                    idx[j] = m
                    m += 1
                else:
                    idx[j] = n
                    n += 1
            for j in range(L - 1):
                if sides[j] == "m":
                    entries_A.append((idx[j + 1], idx[j], 1))
                else:
                    entries_B.append((idx[j + 1], idx[j], 1))
    A = [[0] * m for _ in range(n)]
    B = [[0] * n for _ in range(m)]
    for i, j, v in entries_A:
        A[i][j] = v
    for i, j, v in entries_B:
        B[i][j] = v
    return A, B


def _matmul(X, Y):
    return [[sum(X[i][k] * Y[k][j] for k in range(len(Y))) for j in range(len(Y[0]))] for i in range(len(X))]


def _unimodular(dim, rng, ops):
    P = [[int(i == j) for j in range(dim)] for i in range(dim)]
    Pinv = [row[:] for row in P]
    for _ in range(ops):
        i, j = rng.choice(dim, size=2, replace=False)
        c = int(rng.choice([-3, -2, -1, 1, 2, 3]))
        for row in P:
            row[j] += c * row[i]
        Pinv[i] = [a - c * b for a, b in zip(Pinv[i], Pinv[j])]
    return P, Pinv


def scramble(A0, B0, seed, ops_per_dim):
    """(Q A0 P^-1, P B0 Q^-1) for random unimodular integer P, Q, as object arrays of Python ints."""
    n, m = len(A0), len(B0)
    rng = np.random.default_rng(seed)
    P, Pinv = _unimodular(m, rng, ops_per_dim * m)
    Q, Qinv = _unimodular(n, rng, ops_per_dim * n)
    A = _matmul(_matmul(Q, A0), Pinv)
    B = _matmul(_matmul(P, B0), Qinv)
    return np.array(A, dtype=object), np.array(B, dtype=object)


SMALL = [("chain", 3, "m"), ("chain", 2, "n"), ("inv", 2, 1)]
INV_ONLY = [("inv", 3, 2), ("inv", -1, 1), ("inv", 3, 1)]
BIG = [("inv", 2, 2), ("chain", 1, "m"), ("chain", 5, "n"), ("inv", -1, 3), ("chain", 3, "n"),
       ("chain", 4, "m"), ("chain", 2, "m"), ("inv", -1, 1), ("chain", 5, "m"), ("chain", 1, "n")]
SIDES = [("chain", 4, "m"), ("chain", 6, "n"), ("chain", 4, "n"), ("chain", 6, "m"), ("chain", 1, "m"),
         ("inv", 5, 1)]


def small_pair():
    return scramble(*build_pair(SMALL), seed=1, ops_per_dim=2)


def zero_pair():
    return np.zeros((3, 2), dtype=object), np.zeros((2, 3), dtype=object)


def big_pair():
    return scramble(*build_pair(BIG), seed=7, ops_per_dim=8)


# ----------------------------------------------------------------------------- alternating_ranks
def test_alternating_ranks_normal_case():
    test_input = (*small_pair(), 6)
    expected = [[4, 2, 2, 1, 1, 1, 1], [3, 3, 1, 1, 1, 1, 1]]
    actual = alternating_ranks(*test_input)
    assert np.array_equal(actual, expected)


def test_alternating_ranks_boundary_case():
    # zero maps between C^2 and C^3
    test_input = (*zero_pair(), 3)
    expected = [[2, 0, 0, 0], [3, 0, 0, 0]]
    actual = alternating_ranks(*test_input)
    assert np.array_equal(actual, expected)


def test_alternating_ranks_difficult_case():
    # m = 16, n = 17, scrambled integer entries; products of up to 12 factors
    test_input = (*big_pair(), 12)
    expected = [[16, 14, 10, 9, 7, 6, 6, 6, 6, 6, 6, 6, 6],
                [17, 12, 11, 8, 7, 6, 6, 6, 6, 6, 6, 6, 6]]
    actual = alternating_ranks(*test_input)
    assert np.array_equal(actual, expected)


# ---------------------------------------------------------------------------------- chain_counts
def test_chain_counts_normal_case():
    test_input = small_pair()
    expected = {(3, "m"): 1, (2, "n"): 1}
    actual = chain_counts(*test_input)
    assert actual == expected


def test_chain_counts_boundary_case():
    # A and B invertible: no nilpotent part at all
    test_input = scramble(*build_pair(INV_ONLY), seed=5, ops_per_dim=8)
    expected = {}
    actual = chain_counts(*test_input)
    assert actual == expected


def test_chain_counts_difficult_case():
    test_input = big_pair()
    expected = {(1, "m"): 1, (1, "n"): 1, (2, "m"): 1, (3, "n"): 1, (4, "m"): 1, (5, "m"): 1, (5, "n"): 1}
    actual = chain_counts(*test_input)
    assert actual == expected


def test_chain_counts_both_sides_case():
    # additional regime: chains of the same lengths starting on both sides
    test_input = scramble(*build_pair(SIDES), seed=11, ops_per_dim=8)
    expected = {(1, "m"): 1, (4, "m"): 1, (4, "n"): 1, (6, "m"): 1, (6, "n"): 1}
    actual = chain_counts(*test_input)
    assert actual == expected


# ----------------------------------------------------------------------------- invertible_blocks
def test_invertible_blocks_normal_case():
    test_input = small_pair()
    expected = [(2, 1)]
    actual = invertible_blocks(*test_input)
    assert actual == expected


def test_invertible_blocks_boundary_case():
    # zero pair: no invertible part
    test_input = zero_pair()
    expected = []
    actual = invertible_blocks(*test_input)
    assert actual == expected


def test_invertible_blocks_difficult_case():
    # eigenvalue -1 with Jordan blocks of sizes 1 and 3, eigenvalue 2 with a block of size 2,
    # mixed with nilpotent chains
    test_input = big_pair()
    expected = [(-1, 1), (-1, 3), (2, 2)]
    actual = invertible_blocks(*test_input)
    assert actual == expected


def test_invertible_blocks_repeated_eigenvalue_case():
    # additional regime: purely invertible pair, eigenvalue 3 with blocks of sizes 1 and 2
    test_input = scramble(*build_pair(INV_ONLY), seed=5, ops_per_dim=8)
    expected = [(-1, 1), (3, 1), (3, 2)]
    actual = invertible_blocks(*test_input)
    assert actual == expected


# ----------------------------------------------------------------------------- canonical_form
def test_complete_solution_normal_case():
    test_input = small_pair()
    expected = ([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]],
                [[2, 0, 0], [0, 0, 0], [0, 1, 0], [0, 0, 1]])
    actual = canonical_form(*test_input)
    assert np.array_equal(actual[0], expected[0]) and np.array_equal(actual[1], expected[1])


def test_complete_solution_boundary_case():
    # zero pair: five chains of length 1 (two on C^2, three on C^3), canonical matrices are zero
    test_input = zero_pair()
    expected = (np.zeros((3, 2), dtype=int), np.zeros((2, 3), dtype=int))
    actual = canonical_form(*test_input)
    assert np.array_equal(actual[0], expected[0]) and np.array_equal(actual[1], expected[1])


def test_complete_solution_difficult_case():
    test_input = big_pair()
    # expected canonical pair, built from its blocks in canonical order
    expected = build_pair([("inv", -1, 1), ("inv", -1, 3), ("inv", 2, 2),
                           ("chain", 5, "m"), ("chain", 5, "n"), ("chain", 4, "m"), ("chain", 3, "n"),
                           ("chain", 2, "m"), ("chain", 1, "m"), ("chain", 1, "n")])
    actual = canonical_form(*test_input)
    assert np.array_equal(actual[0], expected[0]) and np.array_equal(actual[1], expected[1])
