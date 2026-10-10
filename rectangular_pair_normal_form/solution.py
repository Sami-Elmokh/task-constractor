"""Complete working code for the task."""

from fractions import Fraction

import numpy as np


def _int_rows(M):
    """Matrix as a list of rows of Python integers (exact, no overflow)."""
    return [[int(x) for x in row] for row in np.asarray(M, dtype=object)]


def _image_basis(vectors):
    """Basis (echelon form, exact Fractions) of the span of the given vectors."""
    basis, pivots = [], []
    for v in vectors:
        w = [Fraction(x) for x in v]
        for b, p in zip(basis, pivots):
            if w[p] != 0:
                c = w[p]
                w = [wi - c * bi for wi, bi in zip(w, b)]
        p = next((i for i, x in enumerate(w) if x != 0), None)
        if p is not None:
            c = w[p]
            w = [x / c for x in w]
            # keep the basis fully reduced so later reductions stay exact and short
            for i, b in enumerate(basis):
                if b[p] != 0:
                    d = b[p]
                    basis[i] = [bi - d * wi for bi, wi in zip(b, w)]
            basis.append(w)
            pivots.append(p)
    return basis


def _apply(M, v):
    return [sum(mij * vj for mij, vj in zip(row, v)) for row in M]


def alternating_ranks(A, B, K):
    """Exact ranks of the alternating products, shape (2, K + 1).

    Row 0: rank of I_m, A, BA, ABA, ... (k factors, starting with A on C^m).
    Row 1: rank of I_n, B, AB, BAB, ... (k factors, starting with B on C^n).
    """
    A, B = _int_rows(A), _int_rows(B)
    n, m = len(A), len(B)
    out = np.zeros((2, K + 1), dtype=int)
    for row, (dim, maps) in enumerate([(m, (A, B)), (n, (B, A))]):
        basis = [[int(i == j) for i in range(dim)] for j in range(dim)]
        out[row, 0] = dim
        for k in range(1, K + 1):
            M = maps[(k - 1) % 2]
            basis = _image_basis([_apply(M, v) for v in basis])
            out[row, k] = len(basis)
    return out


def chain_counts(A, B):
    """Number of nilpotent chains of each length L and starting side ('m' or 'n')."""
    A = _int_rows(A)
    n, m = len(A), len(_int_rows(B))
    K = m + n + 2
    r_m, r_n = alternating_ranks(A, B, K).astype(int)
    total = r_m + r_n          # every chain of length L contributes max(L - k, 0)
    diff = r_m - r_n           # invertible blocks cancel here
    s = diff[:-1] + diff[1:]   # (# m-chains longer than k) - (# n-chains longer than k)
    counts = {}
    for L in range(1, m + n + 1):
        n_all = total[L - 1] - 2 * total[L] + total[L + 1]
        n_diff = s[L - 1] - s[L]
        for side, c in (("m", (n_all + n_diff) // 2), ("n", (n_all - n_diff) // 2)):
            if c:
                counts[(L, side)] = int(c)
    return counts


def invertible_blocks(A, B, max_abs_eigenvalue=50):
    """Jordan blocks (lambda, size) of BA for its nonzero eigenvalues, sorted."""
    A, B = _int_rows(A), _int_rows(B)
    m = len(B)
    C = [[sum(B[i][k] * A[k][j] for k in range(len(A))) for j in range(m)] for i in range(m)]
    # dimension of the invertible part = rank of (BA)^m
    d = int(alternating_ranks(A, B, 2 * m)[0, 2 * m])
    blocks, found = [], 0
    for lam in sorted(range(-max_abs_eigenvalue, max_abs_eigenvalue + 1), key=abs):
        if lam == 0 or found == d:
            continue
        shifted = [[C[i][j] - (lam if i == j else 0) for j in range(m)] for i in range(m)]
        ranks = [m]
        basis = [[int(i == j) for i in range(m)] for j in range(m)]
        while True:
            basis = _image_basis([_apply(shifted, v) for v in basis])
            ranks.append(len(basis))
            if ranks[-1] == ranks[-2]:
                break
        if ranks[-1] == m:
            continue          # lam is not an eigenvalue
        ranks.append(ranks[-1])
        for size in range(1, len(ranks) - 1):
            count = ranks[size - 1] - 2 * ranks[size] + ranks[size + 1]
            blocks += [(lam, size)] * count
            found += size * count
    return sorted(blocks)


def canonical_form(A, B):
    """Use the subproblem functions to produce the canonical block-diagonal pair (A0, B0)."""
    A_rows, B_rows = _int_rows(A), _int_rows(B)
    n, m = len(A_rows), len(B_rows)
    inv = invertible_blocks(A, B)
    chains = chain_counts(A, B)
    chain_list = [key for key in sorted(chains, key=lambda c: (-c[0], c[1])) for _ in range(chains[key])]
    A0 = np.zeros((n, m), dtype=int)
    B0 = np.zeros((m, n), dtype=int)
    im = in_ = 0
    for lam, r in inv:
        for i in range(r):
            A0[in_ + i, im + i] = 1
            B0[im + i, in_ + i] = lam
            if i + 1 < r:
                B0[im + i + 1, in_ + i] = 1
        im += r
        in_ += r
    for L, side in chain_list:
        sides = [("m" if (j % 2 == 0) == (side == "m") else "n") for j in range(L)]
        index = {}
        for j, s in enumerate(sides):
            if s == "m":
                index[j] = im
                im += 1
            else:
                index[j] = in_
                in_ += 1
        for j in range(L - 1):
            if sides[j] == "m":
                A0[index[j + 1], index[j]] = 1
            else:
                B0[index[j + 1], index[j]] = 1
    assert (im, in_) == (m, n)
    return A0, B0
