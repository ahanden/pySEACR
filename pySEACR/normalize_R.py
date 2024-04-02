"""
Find normalization factors.
"""
import numpy as np
from scipy.stats import gaussian_kde
from scipy.signal import decimate
from rpy2.robjects.packages import importr
from rpy2.robjects import vectors

from pySEACR.utils import find_farthest, seq

stats = importr('stats')

class Normalize(object):
    """
    Find normalization factors.
    """

    def __init__(self, exp, ctrl):
        """
        Create a new Normalize object.

        Parameters:
            exp (BDG): Experimental/treatment data
            ctrl (BDG): Control/IgG data
        """
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
        vec_min = min(vec)
        vec_range = cutoff - vec_min
        values = vec[vec <= cutoff]
        values = vectors.FloatVector(values)
        kde = stats.density(values)
        densities = kde[1]
        y_max_index = np.argmax(kde[1])
        x_values = list(kde[0])
        return x_values[y_max_index]

    def constant(self):
        """
        Ratio kernel density inputs with the highest value.

        Creates a kernel density function for each argument, then finds the
        input for each kde with the highest output (estimated).

        Returns:
            float
        """
        return self.max_density(self.exp.vec) / self.max_density(self.ctrl.vec)
