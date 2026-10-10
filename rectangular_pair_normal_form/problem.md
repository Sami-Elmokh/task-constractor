# Normal form of a pair of rectangular matrices

## Scientific objective

A pair of integer matrices $A \in \mathbb{Z}^{n\times m}$ and $B \in \mathbb{Z}^{m\times n}$ defines two linear maps, $A : \mathbb{C}^m \to \mathbb{C}^n$ and $B : \mathbb{C}^n \to \mathbb{C}^m$. Two pairs are **simultaneously equivalent** if

$$A' = Q A P^{-1}, \qquad B' = P B Q^{-1}$$

for some $P \in GL_m(\mathbb{C})$ and $Q \in GL_n(\mathbb{C})$. Every pair is simultaneously equivalent to a direct sum of indecomposable pairs of the two kinds defined below, and this decomposition is unique up to the order of the summands. Compute the decomposition of a given pair and its canonical representative.

## Formulas and assumptions

**Chains.** A chain of length $L \ge 1$ starting on side $s \in \{\texttt{m}, \texttt{n}\}$ has basis vectors $e_0, \dots, e_{L-1}$. The vector $e_j$ lies in $\mathbb{C}^m$ if $j$ is even and $s = \texttt{m}$, or if $j$ is odd and $s = \texttt{n}$; otherwise it lies in $\mathbb{C}^n$. The pair acts by $e_j \mapsto e_{j+1}$ (through $A$ if $e_j \in \mathbb{C}^m$, through $B$ if $e_j \in \mathbb{C}^n$), and $e_{L-1} \mapsto 0$.

**Invertible blocks.** An invertible block $(\lambda, r)$, with $\lambda \neq 0$ and $r \ge 1$, has basis vectors $f_1, \dots, f_r \in \mathbb{C}^m$ and $g_1, \dots, g_r \in \mathbb{C}^n$, with

$$A f_i = g_i, \qquad B g_i = \lambda f_i + f_{i+1} \quad (f_{r+1} = 0).$$

**Alternating products.** $W^{m}_k$ is the product of $k$ factors $A$ and $B$ that alternate and act first on $\mathbb{C}^m$: $W^m_0 = I_m$, $W^m_1 = A$, $W^m_2 = BA$, $W^m_3 = ABA$, and so on. $W^{n}_k$ is defined the same way, starting on $\mathbb{C}^n$: $I_n$, $B$, $AB$, $BAB$, and so on. All ranks are ranks over $\mathbb{Q}$.

**Canonical pair.** The canonical pair $(A_0, B_0)$ is block diagonal. Its blocks are listed in this order:

- first, the invertible blocks, sorted by increasing $(\lambda, r)$;
- then the chains, sorted by decreasing $L$, with chains starting on side $\texttt{m}$ placed before chains starting on side $\texttt{n}$ when $L$ is equal.

Each block occupies the next consecutive coordinates of $\mathbb{C}^m$ and of $\mathbb{C}^n$. Within a block, the basis vectors on each side are taken in this order: $f_1, \dots, f_r$ and $g_1, \dots, g_r$ for an invertible block, and the vectors of that side by increasing $j$ for a chain. $A_0$ and $B_0$ are the matrices of $A$ and $B$ in this basis.

**Assumptions.**

- $1 \le m, n \le 40$.
- The entries are integers of arbitrary size.
- The nonzero eigenvalues of $BA$ are integers with $|\lambda| \le 50$.

| Symbol or notation | Meaning | Units, if any |
| --- | --- | --- |
| $m$, $n$ | dimensions of the two spaces | — |
| $A$, $B$ | the pair, of shapes $n\times m$ and $m \times n$ | — |
| $L$, side | length and starting side of a chain | — |
| $(\lambda, r)$ | eigenvalue and size of an invertible block | — |
| $W^m_k$, $W^n_k$ | alternating products with $k$ factors | — |
| $(A_0, B_0)$ | canonical pair | — |

## Required functions

### `alternating_ranks(A, B, K)`

- Input: `A`, `B`, `numpy` arrays of `dtype=object` holding Python integers, of shapes `(n, m)` and `(m, n)`; `K`, integer $\ge 0$.
- Output: `numpy` integer array of shape `(2, K + 1)`. Entry `[0, k]` is $\operatorname{rank} W^m_k$ and entry `[1, k]` is $\operatorname{rank} W^n_k$.
- Purpose: Rank profile of the alternating products.

### `chain_counts(A, B)`

- Input: `A`, `B` as above.
- Output: `dict` mapping `(L, side)`, with `side` equal to `"m"` or `"n"`, to the number of chains of length `L` starting on that side. Only positive counts are included.
- Purpose: Nilpotent part of the decomposition.

### `invertible_blocks(A, B)`

- Input: `A`, `B` as above.
- Output: `list` of tuples `(lam, r)` of Python integers, one per invertible block, with repetitions, sorted increasingly.
- Purpose: Invertible part of the decomposition.

### `canonical_form(A, B)`

- Input: `A`, `B` as above.
- Output: tuple `(A0, B0)` of `numpy` integer arrays of shapes `(n, m)` and `(m, n)`.
- Purpose: Use `alternating_ranks`, `chain_counts` and `invertible_blocks` to produce the canonical pair.

## Allowed libraries

`numpy` and the Python standard library.
