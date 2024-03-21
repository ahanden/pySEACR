import numpy as np
import warnings
from scipy.stats import gaussian_kde, ecdf


def dist2d(x):
    """Find the two dimensional distance for the first two elements of a vector.

    This is actually just a stupid wrapper for x[0] - x[1]

    Parameters:
        x (collection): Values to calculate by

    Returns:
        float
    """
    a, c = x
    return a - c

def seq(length):
    """
    Create a list with values between 0 and 1 with the given length.

    Parameters:
        length (int): List length

    Returns:
        list
    """
    return [1 - _/(length - 1) for _ in range(length)]

def get_farthest_value(vec):
    """
    Find the value with the largest distance.

    Distance is length / rank - value / max.

    Parameters:
        vec (collection): Numbers to search

    Returns:
        An elment of vec if the distance is > 0.9 * max distance. Else None.
    """
    vec.sort(reverse=True)
    vec_max = max(vec)
    count = seq(len(vec))
    quant = [_ / vec_max for _ in vec]
    distance = [c - q for c, q in zip(count, quant)]

    farthest_index = max(range(len(vec)), key=lambda _: distance[_])
    if distance[farthest_index] <= 0.9 * max(distance):
        return None
    return vec[farthest_index]

def find_farthest(vec):
    """
    Find the value with the largest distance, or make one up.

    Parameters:
        vec (collection): Numbers to search

    Returns:
        The max of get_farthest_value(vec) and vec's 90th quantile.
    """
    value = get_farthest_value(vec)
    ninetieth = np.quantile(vec, 0.9)
    if value is None or value > ninetieth:
        return value
    return ninetieth

def density_subset(vec, max_value):
    """
    Compute kernel density for a subset of a vector.

    Parameters:
        vec (collection): Numbers to estimate by
        max_value (float): Discard members of vec higher than this value

    Returns:
        scipy.stats.guassian_kde
    """
    return gaussian_kde([_ for _ in vec if _ < max_value])

def density_max(vec):
    """
    Estimate input for a kernel density function that has the largest output.

    Parameters:
        vec (collection): Data used to create the kde

    Returns:
        float
    """
    cutoff = find_farthest(vec)
    kde = gaussian_kde([_ for _ in vec if _ < cutoff])
    vec_min = min(vec)
    vec_range = cutoff - vec_min
    values = [vec_min + _ * vec_range for _ in seq(len(vec) // 10000)]
    densities = kde.evaluate(values)
    max_index = max(range(len(values)), key=lambda _: densities[_])
    return values[max_index]

def find_constant(vec_1, vec_2):
    """
    The ratio of kernel density inputs with the highest value.

    Creates a kernel density function for each argument, then finds the
    input for each kde with the highest output (estimated).

    Parameters:
        vec_1 (collection): Numerator values
        vec_2 (collection): Denominator values

    Returns:
        float
    """
    return density_max(vec_1) / density_max(vec_2)

def combine(list_1, list_2):
    """
    Find the sorted union of two lists.

    Parameters:
        list_1 (list): First list
        list_2 (list): Second list

    Returns:
        Sorted list
    """
    return sorted(set(list_1).union(set(list_2)))

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
    return 1 + ecdf_2.evaluate(x) - ecdf_1.evaluate(x)

def find_relaxed_thresh(exp, ctrl):
    """
    Find the relaxed threshold.

    Parameters:
        exp (collection): Experimental data
        ctrl (collection): Control data

    Returns:
        float
    """
    both = combine(exp, ctrl)
    inner = pct_remain_vec(exp, ctrl, both)
    indices = [_ for _ in range(len(both)) if inner[_] < 1]
    outer = pct_remain_vec(exp, ctrl, [both[_] for _ in indices])[:-1]
    return both[outer.argmax()]

def find_stringent_thresh(exp, ctrl, relaxed):
    """
    Find the stringent threshold

    Parameters:
        exp (collection): Experimental data
        ctrl (collection): Control data

    Returns:
        float
    """
    low_values = [_ for _ in combine(exp, ctrl) if _ <= relaxed]
    relaxed_pct = pct_remain_vec(exp, ctrl, relaxed)
    low_pct = pct_remain_vec(exp, ctrl, low_values)
    search_values = abs((relaxed_pct + min(low_pct)) / 2 - low_pct)
    thresh_check = low_values[search_values.argmin()]

    if relaxed == thresh_check:
        return relaxed
    high_values = np.array([_ for _ in low_values if _ > thresh_check])
    high_max = max(high_values)
    high_range = high_max - min(high_values)
    search_values = abs(high_values - high_max + high_range / 2)
    return high_values[search_values.argmin()]

def diff(vec):
    """
    Return suitably lagged and iterated differences.

    Parameters:
        vec (collection): Data to iterate over

    Returns:
        The difference between every neighboring elements as a list
    """
    return [a - b  for a, b in zip(vec[1:], vec[:-1])]

def find_best_quantile(vec):
    """
    Find the lowest quantile above 0.99 that does not have a value of zero.

    Parameters:
        vec (collection): Values to calculate by

    Returns:
        float
    """
    digits = 1
    output = 0
    while output == 0:
        digits += 1
        output = np.nanquantile(vec, 1 - 1 / 10 ** digits)
    return output

def spurious_values(exp, ctrl):
    """
    Find abnormally low values.

    Parameters:
        exp (collection): Experimental data
        ctrl (collection): Control data

    Returns:
        list
    """
    both = combine(exp, ctrl)
    delta = [abs(_) for _ in diff(pct_remain_vec(exp, ctrl, both))]
    best = find_best_quantile(delta)
    indices = [_ for _ in range(len(delta)) if delta[_] != 0 and delta[_] < best]
    return [both[_] for _ in indices]

def thresh_check(exp, ctrl):
    """
    Check for spurious values and conditionally adjust thresholds.

    Parameters:
        exp (collection): Experimental data
        ctrl (collection): Control data
    """
    sp_values = spurious_values(exp, ctrl)
    sp_pct = pct_remain_vec(exp, ctrl, sp_values)
    indices = [_ for _ in range(len(sp_pct)) if sp_pct[_] < 1]
    search_values = pct_remain_vec(
        exp,
        ctrl,
        np.array([sp_values[_] for _ in indices]),
    )
    sp_thresh = sp_values[search_values.argmax()]
    check_thresh = find_stringent_thresh(exp, ctrl, sp_thresh)
    both = combine(exp, ctrl)
    both_pct = pct_remain_vec(exp, ctrl, both)
    indices = [_ for _ in range(len(both_pct)) if both_pct[_] < 1]
    old_values = pct_remain_vec(exp, ctrl, [both[_] for _ in indices])
    if max(search_values) / max(old_values) > 0.95:
        return sp_thresh, check_thresh

def find_global_thresh(exp, ctrl):
    """
    Find the global threshold

    Parameters:
        exp (collection): Experimental data
        ctrl (collection): Control data

    Returns:
        float
    """
    both = combine(exp, ctrl)
    pct = pct_remain_max(exp, ctrl, both)
    if True in (_ > 1 for _ in pct):
        indices = [_ for _ in range(len(pct)) if pct[_] > 1]
        return min(both[_] for _ in indices)
    return 1


def get_static_thresh(vec, ctrl):
    """
    Find the threshold from a vector.

    Parameters:
        vec (collection): Vector to calculate for
        ctrl (float): Control/max values

    Returns:
        float
    """
    model = ecdf(vec).cdf
    quantiles = 1 - model.evaluate(vec)
    indices = [_ for _ in range(len(quantiles)) if quantiles[_] < ctrl]
    return min([vec[_] for _ in indices])

def get_static_thresholds(exp_vec, exp_max, ctrl):
    """
    Find thresholds given an FDR rate.

    Parameters:
        exp_vec (collection): Experimental vector data
        exp_max (collection): Experimental peak data
        ctrl (float): FDR rate

    Returns:
        Relaxed, stringent, and global thresholds as well as FDR.
    """
    relaxed = get_static_thresh(exp_vec, ctrl)
    stringent = get_static_thresh(exp_max, ctrl)
    return (relaxed, stringent, 0, ctrl)
