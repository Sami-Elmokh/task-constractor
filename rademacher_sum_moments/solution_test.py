"""Tests for every required function. Repeat this pattern for each function."""

from fractions import Fraction

import numpy as np

from solution import even_moment, khintchine_ratio, rademacher_law, tail_probability

# Laws, moments and tail probabilities are exact rationals or integers, so they are compared with ==.
# The ratio is irrational in general; the references have 20 correct digits and a correct
# double-precision evaluation reaches about 1e-15, so RTOL = 1e-12 is the accuracy required in
# problem.md.
RTOL = 1e-12

C3 = np.array([1, 2, 3])
C20 = np.array([36, 10, 15, 0, 6, 22, 27, -22, -36, -16, -17, -22, 33, -40, 0, 26, -30, 24, -31, -3])
# 60 coefficients, three magnitudes of order 1e12
C60 = np.array([10 ** 12] * 20 + [-2 * 10 ** 12] * 20 + [3 * 10 ** 12] * 20, dtype=np.int64)
ONES = np.ones(1000, dtype=np.int64)


def big_coefficients():
    # 100000 coefficients, about 78500 distinct magnitudes up to 1e6, about one third negative, 101 zeros
    i = np.arange(100000, dtype=np.int64)
    c = ((i * i * 7919 + i * 104729) % 1000000 + 1) * np.where(i % 3 == 1, -1, 1)
    c[i % 997 == 0] = 0
    return c


def test_rademacher_law_normal_case():
    expected = {-6: Fraction(1, 8), -4: Fraction(1, 8), -2: Fraction(1, 8), 0: Fraction(1, 4),
                2: Fraction(1, 8), 4: Fraction(1, 8), 6: Fraction(1, 8)}
    assert rademacher_law(C3) == expected


def test_rademacher_law_boundary_case():
    # zero coefficients only double the number of sign patterns
    assert rademacher_law(np.array([0, -5, 0, 0])) == {-5: Fraction(1, 2), 5: Fraction(1, 2)}
    assert rademacher_law(np.array([0, 0, 0])) == {0: Fraction(1)}


def test_rademacher_law_cancellation_case():
    # opposite coefficients and different sign patterns giving the same value
    expected = {-12: Fraction(1, 8), -6: Fraction(1, 4), 0: Fraction(1, 4), 6: Fraction(1, 4), 12: Fraction(1, 8)}
    assert rademacher_law(np.array([3, -3, 0, 6])) == expected


def test_rademacher_law_difficult_case():
    law = rademacher_law(C20)
    assert len(law) == 399
    assert sum(law.values()) == 1
    expected = {-416: Fraction(1, 262144), -266: Fraction(21, 65536), -134: Fraction(921, 262144),
                -2: Fraction(479, 65536), 0: Fraction(951, 131072), 130: Fraction(963, 262144),
                262: Fraction(49, 131072), 404: Fraction(1, 262144), 416: Fraction(1, 262144)}
    assert {v: law.get(v) for v in expected} == expected
    assert 1 not in law and 417 not in law


def test_rademacher_law_large_coefficients_case():
    law = rademacher_law(C60)
    assert len(law) == 121
    assert min(law) == -120 * 10 ** 12 and max(law) == 120 * 10 ** 12
    assert law[0] == Fraction(13657533155376919, 288230376151711744)
    assert law[2 * 10 ** 12] == Fraction(3390680379966455, 72057594037927936)
    assert law[120 * 10 ** 12] == Fraction(1, 2 ** 60)
    assert 10 ** 12 not in law


def test_even_moment_normal_case():
    assert even_moment(C3, 2) == 392


def test_even_moment_boundary_case():
    assert even_moment(np.array([7, -3, 0, 2]), 1) == 62
    assert even_moment(np.array([5]), 6) == 5 ** 12


def test_even_moment_equal_coefficients_case():
    # sum of 1000 signs; the Gaussian value 19!! * 1000^10 = 654729075 * 10^30 is wrong
    assert even_moment(ONES, 10) == 635352202681178618159158926784135168000


def test_even_moment_difficult_case():
    c = big_coefficients()
    assert even_moment(c, 2) == 3339175324142214121669644842601845
    expected = int(
        "931456419999085964380909407391336890864387280235128531437925637321954295822232192538658285"
        "643789890695414781704255315787150157435018951056537896867643381769552790780419620095144966"
        "678366438330496046472819605924030349809548737917110964599930992046697634766604974114517987"
        "595214225574816250638736314599774706093899973233259077596133447643554738612647565065")
    assert even_moment(c, 20) == expected


def test_tail_probability_normal_case():
    assert tail_probability(C3, 1) == Fraction(3, 4)


def test_tail_probability_boundary_case():
    # t equal to a value of |S|, t below 0 and t at the maximum
    assert tail_probability(C3, 2) == Fraction(1, 2)
    assert tail_probability(C3, Fraction(7, 2)) == Fraction(1, 2)
    assert tail_probability(C3, -1) == 1
    assert tail_probability(C3, 6) == 0


def test_tail_probability_difficult_case():
    assert tail_probability(C60, 2 * 10 ** 13) == Fraction(60698409131179605, 288230376151711744)
    assert tail_probability(C60, 2 * 10 ** 13 - 1) == Fraction(37126980637086433, 144115188075855872)


def test_khintchine_ratio_normal_case():
    assert np.isclose(khintchine_ratio(C3, 3.0), 1.1118515716615910068, rtol=RTOL, atol=0)


def test_khintchine_ratio_boundary_case():
    assert np.isclose(khintchine_ratio(C20, 2.0), 1.0, rtol=RTOL, atol=0)
    assert np.isclose(khintchine_ratio(np.array([1, 1]), 1.0), 0.7071067811865475244, rtol=RTOL, atol=0)


def test_khintchine_ratio_odd_and_fractional_case():
    assert np.isclose(khintchine_ratio(C20, 7.0), 1.6236854117039805981, rtol=RTOL, atol=0)
    assert np.isclose(khintchine_ratio(C60, 1.5), 0.90524820056160317613, rtol=RTOL, atol=0)


def test_khintchine_ratio_difficult_case():
    assert np.isclose(khintchine_ratio(ONES, 4.0), 1.3158546124259147867, rtol=RTOL, atol=0)
    assert np.isclose(khintchine_ratio(big_coefficients(), 40.0), 3.8689958980990702001, rtol=RTOL, atol=0)
