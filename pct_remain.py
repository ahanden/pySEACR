"""
PCT Remain functions from SEACR.
"""
import warnings

import numpy as np
from scipy.stats import ecdf


def pct_remain_vec(vec_1, vec_2, x):
    """
    Find the ratio of vectors' distributions.

    Parameters:
        vec_1 (collection): Numerator
        vec_2 (collection): Denominator
        x (float): Value to test on

    Returns:
        float
    """
    both = np.concatenate((vec_1, vec_2))
    ecdf_1 = ecdf(vec_1).cdf
    ecdf_both = ecdf(both).cdf
    top = len(vec_1) - ecdf_1.evaluate(x) * len(vec_1)
    bottom = len(both) - ecdf_both.evaluate(x) * len(both)
    with warnings.catch_warnings(action="ignore"):
        return top / bottom


def pct_remain_max(vec_1, vec_2, x):
    """
    Find the ratio of two maxs' distributions.

    Parameters:
        vec_1 (collection): First data set
        vec_2 (collection): Second data set
        x (float): Value to test

    Returns:
        float
    """
    ecdf_1 = ecdf(vec_1).cdf
    ecdf_2 = ecdf(vec_2).cdf
    return 1 - (ecdf_1.evaluate(x) - ecdf_2.evaluate(x))
