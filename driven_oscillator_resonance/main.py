"""Run one example of the completed solution."""

import numpy as np

from solution import driven_response


def main():
    # force exp(cos(t - 0.4)) sampled at N = 49 points; it contains a resonant component
    N = 49
    mu = np.exp(np.cos(2 * np.pi * np.arange(N) / N - 0.4))
    t = np.array([0.0, 2.0, 7.3, 18.0, 33.3, 60.0])
    print("t    =", t)
    print("y(t) =", driven_response(mu, t))


if __name__ == "__main__":
    main()
