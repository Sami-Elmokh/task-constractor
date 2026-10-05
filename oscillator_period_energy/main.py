"""Run one example of the completed solution."""

import numpy as np

from solution import period_energy_curve


def main():
    # tilted potential f(x) = -cos x + 0.2 x around its central well, from small oscillations to escape
    example_input = (lambda x: -np.cos(x) + 0.2 * x, lambda x: np.sin(x) + 0.2, -np.arcsin(0.2),
                     [0.1, 0.5, 0.9, 0.999, 0.999999])
    result = period_energy_curve(*example_input)
    print(result)


if __name__ == "__main__":
    main()
