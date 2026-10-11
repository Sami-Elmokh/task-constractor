"""Run one example of the completed solution."""

import numpy as np

from solution import khintchine_ratio


def main():
    # 1000 coefficients equal to 1, fourth-moment Khintchine ratio
    example_input = (np.ones(1000, dtype=np.int64), 4.0)
    result = khintchine_ratio(*example_input)
    print(result)


if __name__ == "__main__":
    main()
