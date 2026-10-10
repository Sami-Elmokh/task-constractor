"""Run one example of the completed solution."""

import numpy as np

from solution import driven_response


def main():
    # force exp(cos(t - 0.4)) sampled at N = 49 points, oscillator 5e-9 away from its second harmonic
    example_input = (np.exp(np.cos(2 * np.pi * np.arange(49) / 49 - 0.4)), 2.0 + 5e-9,
                     [0.0, 5.3, 31.0, 97.0, 200.0])
    result = driven_response(*example_input)
    print(result)


if __name__ == "__main__":
    main()
