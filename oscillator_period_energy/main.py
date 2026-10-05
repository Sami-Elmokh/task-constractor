"""Run one example of the completed solution."""

import numpy as np

from solution import period_energy_curve


def main():
    # tilted potential f(x) = -cos x + 0.2 x around its central well
    f = lambda x: -np.cos(x) + 0.2 * x
    fp = lambda x: np.sin(x) + 0.2
    x0 = -np.arcsin(0.2)
    q = [0.1, 0.5, 0.9, 0.999, 0.999999]
    print("q    =", q)
    print("T(q) =", period_energy_curve(f, fp, x0, q))


if __name__ == "__main__":
    main()
