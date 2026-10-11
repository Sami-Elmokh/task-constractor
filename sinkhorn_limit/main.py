"""Run one example of the completed solution."""

import numpy as np

from solution import sinkhorn_limit


def main():
    # a 5x5 matrix with a positive diagonal but without total support
    example_input = np.array([[2., 1, 3, 0, 0], [0, 1, 1, 0, 0], [0, 4, 2, 0, 0], [1, 2, 0, 3, 1], [5, 0, 1, 1, 2]])
    result = sinkhorn_limit(example_input)
    print(result)


if __name__ == "__main__":
    main()
