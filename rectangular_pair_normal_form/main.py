"""Run one example of the completed solution."""

import numpy as np

from solution import canonical_form


def main():
    # a pair A : C^4 -> C^3, B : C^3 -> C^4 given in scrambled integer coordinates
    example_input = (np.array([[0, 0, 0, 0], [1, 1, -3, -3], [-1, -1, 3, 4]], dtype=object),
                     np.array([[27, 11, 6], [45, 4, 6], [16, 1, 2], [8, 4, 2]], dtype=object))
    result = canonical_form(*example_input)
    print(result)


if __name__ == "__main__":
    main()
