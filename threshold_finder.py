"""
Calculate thresholds for peak identification.
"""
import numpy as np
from scipy.stats import ecdf

from pySEACR.auc_from_bdg import BDG
from pySEACR.pct_remain import pct_remain_max, pct_remain_vec
from pySEACR.utils import combine, diff, find_best_quantile


class ThresholdFinder(object):
    """
    Calculate thresholds for peak identification.
    """
    def __init__(self, exp, ctrl):
        """
        Create a new ThresholdFinder object.

        Parameters:
            exp (BDG): Experimental/treatment BDG data
            ctrl (BDG): Control/IgG BDG data
        """
        self.exp = exp
        self.ctrl = ctrl
        self.both = combine(exp.vec, ctrl.vec)

    def spurious_values(self):
        """
        Find abnormally low values.

        Returns:
            list
        """
        delta = diff(pct_remain_vec(self.exp.vec, self.ctrl.vec, self.both))
        best = find_best_quantile(delta)
        sp = self.both[:-1][delta < best]
        return sp[sp != 0]

    def relaxed(self):
        """
        Find the relaxed threshold.

        Returns:
            float
        """
        if not isinstance(self.ctrl, BDG):
            return self.static(self.exp.vec)
        inner = pct_remain_vec(self.exp.vec, self.ctrl.vec, self.both)
        indices = [_ for _ in range(len(self.both)) if inner[_] < 1]
        outer = pct_remain_vec(self.exp.vec, self.ctrl.vec, [self.both[_] for _ in indices])[:-1]
        return self.both[outer.argmax()]

    def stringent(self, relaxed_thresh):
        """
        Find the stringent threshold

        Parameters:
            relaxed_thresh (float): Relaxed threshold

        Returns:
            float
        """
        if not isinstance(self.ctrl, BDG):
            return self.static(self.exp.max)
        low_values = [_ for _ in combine(self.exp.vec, self.ctrl.vec) if _ <= relaxed_thresh]
        relaxed_thresh_pct = pct_remain_vec(self.exp.vec, self.ctrl.vec, relaxed_thresh)
        low_pct = pct_remain_vec(self.exp.vec, self.ctrl.vec, low_values)
        search_values = abs((relaxed_thresh_pct + min(low_pct)) / 2 - low_pct)
        thresh_check = low_values[search_values.argmin()]

        if relaxed_thresh == thresh_check:
            return relaxed_thresh
        high_values = np.array([_ for _ in low_values if _ > thresh_check])
        high_max = max(high_values)
        high_range = high_max - min(high_values)
        search_values = abs(high_values - high_max + high_range / 2)
        return high_values[search_values.argmin()]

    def genome(self):
        """
        Find the genome threshold

        Returns:
            float
        """
        both = combine(self.exp.max, self.ctrl.max)
        pct = pct_remain_max(self.exp.max, self.ctrl.max, both)
        if any(pct > 1):
            return min(both[pct > 1])
        return 1

    def static(self, vector):
        """
        Find the threshold from a vector.

        Returns:
            float
        """
        model = ecdf(vector).cdf
        quantiles = 1 - model.evaluate(vector)
        indicies = [_ for _ in range(len(quantiles)) if quantiles[_] < self.ctrl]
        return min([vector[_] for _ in indicies])

    def thresh_check(self):
        """
        Check for spurious values and conditionally adjust thresholds.
        """
        sp_values = self.spurious_values()
        sp_pct = pct_remain_vec(self.exp.vec, self.ctrl.vec, sp_values)
        indices = [_ for _ in range(len(sp_pct)) if sp_pct[_] < 1]
        search_values = pct_remain_vec(
            self.exp.vec,
            self.ctrl.vec,
            np.array([sp_values[_] for _ in indices]),
        )
        sp_thresh = sp_values[search_values.argmax()]
        check_thresh = self.stringent(sp_thresh)
        both_pct = pct_remain_vec(self.exp.vec, self.ctrl.vec, self.both)
        old_values = pct_remain_vec(self.exp.vec, self.ctrl.vec, self.both[both_pct < 1])
        if max(search_values) / max(old_values) > 0.95:
            return sp_thresh, check_thresh
        return None
