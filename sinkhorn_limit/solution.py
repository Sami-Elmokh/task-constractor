"""Complete working code for the task."""

import numpy as np


def positive_diagonal_support(A):
    """Boolean mask of the positive entries of A that lie on at least one positive diagonal."""
    P = np.asarray(A, dtype=float) > 0
    n = P.shape[0]
    # one perfect matching of the bipartite support graph (augmenting paths)
    row_of_col = [-1] * n

    def augment(i, seen):
        for j in np.flatnonzero(P[i]):
            if not seen[j]:
                seen[j] = True
                if row_of_col[j] == -1 or augment(row_of_col[j], seen):
                    row_of_col[j] = i
                    return True
        return False

    for i in range(n):
        if not augment(i, [False] * n):
            raise ValueError("A has no positive diagonal")
    # an edge (i, j) lies on a positive diagonal iff it lies on an alternating cycle, i.e. iff
    # row i and the row matched to column j are in the same strongly connected component of
    # the digraph i -> row_of_col[j] for every positive entry (i, j)
    adj = [[row_of_col[j] for j in np.flatnonzero(P[i])] for i in range(n)]
    index, low, on_stack, comp = [None] * n, [0] * n, [False] * n, [-1] * n
    stack, counter, n_comp = [], 0, 0
    for s in range(n):
        if index[s] is not None:
            continue
        work = [(s, 0)]
        index[s] = low[s] = counter
        counter += 1
        stack.append(s)
        on_stack[s] = True
        while work:
            v, k = work[-1]
            if k < len(adj[v]):
                work[-1] = (v, k + 1)
                w = adj[v][k]
                if index[w] is None:
                    index[w] = low[w] = counter
                    counter += 1
                    stack.append(w)
                    on_stack[w] = True
                    work.append((w, 0))
                elif on_stack[w]:
                    low[v] = min(low[v], index[w])
            else:
                work.pop()
                if work:
                    low[work[-1][0]] = min(low[work[-1][0]], low[v])
                if low[v] == index[v]:
                    while True:
                        w = stack.pop()
                        on_stack[w] = False
                        comp[w] = n_comp
                        if w == v:
                            break
                    n_comp += 1
    mask = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in np.flatnonzero(P[i]):
            mask[i, j] = comp[i] == comp[row_of_col[j]]
    return mask


def doubly_stochastic_scaling(A):
    """Doubly stochastic matrix diag(r) A diag(c), for A with total support (Newton in log variables)."""
    A = np.asarray(A, dtype=float)
    n = A.shape[0]
    with np.errstate(divide="ignore"):
        logA = np.where(A > 0, np.log(np.where(A > 0, A, 1.0)), -np.inf)

    def scaled(u, v):
        return np.exp(logA + u[:, None] + v[None, :])

    u, v = np.zeros(n), np.zeros(n)
    for _ in range(10):                       # a few Sinkhorn sweeps as a starting point
        u -= np.log(scaled(u, v).sum(1))
        v -= np.log(scaled(u, v).sum(0))
    for _ in range(100):
        S = scaled(u, v)
        F = np.concatenate([S.sum(1) - 1, S.sum(0) - 1])
        err = np.abs(F).max()
        if err < 1e-15:
            break
        J = np.block([[np.diag(S.sum(1)), S], [S.T, np.diag(S.sum(0))]])
        # J is singular along u += t, v -= t on each block: minimum-norm step
        step = np.linalg.lstsq(J, -F, rcond=None)[0]
        t = 1.0
        while t > 1e-8:
            S_new = scaled(u + t * step[:n], v + t * step[n:])
            if np.abs(np.concatenate([S_new.sum(1) - 1, S_new.sum(0) - 1])).max() < err:
                break
            t /= 2
        u, v = u + t * step[:n], v + t * step[n:]
    return scaled(u, v)


def sinkhorn_limit(A):
    """Use the subproblem functions to produce the limit of the Sinkhorn iteration started at A."""
    A = np.asarray(A, dtype=float)
    mask = positive_diagonal_support(A)
    return doubly_stochastic_scaling(np.where(mask, A, 0.0))
