"""Run one example of the completed solution."""

import numpy as np

from solution import driven_response


def main():
    # force exp(cos(t - 0.4)) sampled at N = 49 points, evaluated at six times
    example_input = (np.exp(np.cos(2 * np.pi * np.arange(49) / 49 - 0.4)), [0.0, 2.0, 7.3, 18.0, 33.3, 60.0])
    result = driven_response(*example_input)
    print(result)


if __name__ == "__main__":
    main()
