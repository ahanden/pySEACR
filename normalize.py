import numpy as np
import warnings
from scipy.stats import gaussian_kde, ecdf

class Normalize(object):
    def __init__(self, exp, ctrl):
        self.exp = exp
        self.ctrl = ctrl

    def max_density(self, vec):
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

    def constant(self):
        """
        The ratio of kernel density inputs with the highest value.

        Creates a kernel density function for each argument, then finds the
        input for each kde with the highest output (estimated).

        Returns:
            float
        """
        return max_density(self.exp) / max_density(self.ctrl)
