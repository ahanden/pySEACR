"""
Miscellaneous functions for pySEACR
"""
import warnings

import numpy as np


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

def combine(list_1, list_2):
    """
    Find the sorted union of two lists.

    Parameters:
        list_1 (list): First list
        list_2 (list): Second list

    Returns:
        Sorted list
    """
    return np.unique(np.concatenate((list_1, list_2)))

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
