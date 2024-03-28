"""PCT Remain functions from SEACR."""
import warnings

import numpy as np
from scipy.stats import ecdf


def pct_remain_vec(exp, ctrl, value):
    """
    Find the ratio of vectors' distributions.

    Parameters:
        exp (collection): Numerator
        ctrl (collection): Denominator
        value (float): Value to test on

    Returns:
        float
    """
    both = np.concatenate((exp, ctrl))
    ecdf_exp = ecdf(exp).cdf
    ecdf_both = ecdf(both).cdf
    top = len(exp) - ecdf_exp.evaluate(value) * len(exp)
    bottom = len(both) - ecdf_both.evaluate(value) * len(both)
    with warnings.catch_warnings(action='ignore'):
        return top / bottom


def pct_remain_max(exp, ctrl, value):
    """
    Find the ratio of two maxs' distributions.

    Parameters:
        exp (collection): First data set
        ctrl (collection): Second data set
        value (float): Value to test

    Returns:
        float
    """
    ecdf_exp = ecdf(exp).cdf
    ecdf_ctrl = ecdf(ctrl).cdf
    return 1 - (ecdf_exp.evaluate(value) - ecdf_ctrl.evaluate(value))
