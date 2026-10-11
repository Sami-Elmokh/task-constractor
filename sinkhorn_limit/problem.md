# Limit of the Sinkhorn iteration

## Scientific objective

A nonnegative square matrix $A$ is rescaled into a doubly stochastic matrix, meaning a point of the Birkhoff polytope, by the Sinkhorn iteration: alternately divide every row by its sum and every column by its sum. This operation appears in optimal transport, ranking and contingency-table analysis. Compute the limit of this iteration exactly to the required accuracy, for any nonnegative matrix that admits a positive diagonal.

## Formulas and assumptions

**Sinkhorn iteration.** Starting from $S_0 = A$, for $t \ge 0$:

$$R_{t} = \mathrm{diag}(S_t \mathbf{1})^{-1} S_t, \qquad S_{t+1} = R_t\, \mathrm{diag}(\mathbf{1}^{T} R_t)^{-1},$$

that is, a row normalization followed by a column normalization. For every $A$ with at least one positive diagonal, $S_t$ converges to a doubly stochastic matrix $S_\infty$ (Sinkhorn and Knopp).

**Positive diagonals.** A positive diagonal of $A$ is a permutation $\sigma$ with $A_{i\sigma(i)} > 0$ for all $i$. An entry $(i, j)$ lies on a positive diagonal if some positive diagonal has $\sigma(i) = j$. $A$ has total support if every positive entry lies on a positive diagonal.

**Doubly stochastic scaling.** If $A$ has total support, there are positive diagonal matrices $D_1, D_2$ such that $D_1 A D_2$ is doubly stochastic, and the matrix $D_1 A D_2$ is unique (even though $D_1$ and $D_2$ are not).

**Requirements.**

- $A$ is real, nonnegative, of size $n \times n$ with $1 \le n \le 12$, and has at least one positive diagonal.
- Positive entries may range from $10^{-8}$ to $10^{8}$.
- All returned matrices must be accurate to an absolute error of $10^{-10}$ in every entry.

| Symbol or notation | Meaning | Units, if any |
| --- | --- | --- |
| $A$ | nonnegative input matrix | — |
| $S_t$ | Sinkhorn iterate after $t$ sweeps | — |
| $S_\infty$ | limit of the Sinkhorn iteration | — |
| $\sigma$ | permutation of $\{1, \dots, n\}$ | — |
| $D_1, D_2$ | positive diagonal scaling matrices | — |
| $\mathbf{1}$ | vector of ones | — |

## Required functions

### `positive_diagonal_support(A)`

- Input: `A`, `numpy` array of floats of shape `(n, n)`, nonnegative, with at least one positive diagonal.
- Output: `numpy` boolean array of shape `(n, n)`, `True` exactly at the positive entries of $A$ that lie on a positive diagonal.
- Purpose: Determine which entries of $A$ belong to some positive diagonal.

### `doubly_stochastic_scaling(A)`

- Input: `A`, `numpy` array of floats of shape `(n, n)`, nonnegative, with total support.
- Output: `numpy` array of floats of shape `(n, n)`, the doubly stochastic matrix $D_1 A D_2$.
- Purpose: Scale a matrix with total support to a doubly stochastic matrix.

### `sinkhorn_limit(A)`

- Input: `A`, `numpy` array of floats of shape `(n, n)`, nonnegative, with at least one positive diagonal.
- Output: `numpy` array of floats of shape `(n, n)`, the limit $S_\infty$.
- Purpose: Use `positive_diagonal_support` and `doubly_stochastic_scaling` to produce the limit of the Sinkhorn iteration.

## Allowed libraries

`numpy` and the Python standard library.
