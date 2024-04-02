"""Calculate thresholds for peak identification."""
from scipy.stats import ecdf

from pySEACR.auc_from_bdg import BDG
from pySEACR.pct_remain import pct_remain_max, pct_remain_vec
from pySEACR.utils import combine, diff, find_best_quantile


class ThresholdFinder():
    """Calculate thresholds for peak identification."""

    def __init__(self, exp, ctrl):
        """
        Create a new ThresholdFinder object.

        Parameters:
            exp (BDG): Experimental/treatment BDG data
            ctrl (BDG): Control/IgG BDG data
        """
        self.exp = exp
        self.ctrl = ctrl
        if isinstance(ctrl, BDG):
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
        outer = pct_remain_vec(
            self.exp.vec,
            self.ctrl.vec,
            self.both[inner < 1],
        )[:-1]
        return self.both[outer.argmax()]

    def stringent(self, relaxed_thresh):
        """
        Find the stringent threshold.

        Parameters:
            relaxed_thresh (float): Relaxed threshold

        Returns:
            float
        """
        if not isinstance(self.ctrl, BDG):
            return self.static(self.exp.max)
        low_values = self.both[self.both <= relaxed_thresh]
        low_pct = pct_remain_vec(self.exp.vec, self.ctrl.vec, low_values)
        search_values = abs(
            (
                pct_remain_vec(
                    self.exp.vec,
                    self.ctrl.vec,
                    relaxed_thresh,
                ) +
                min(low_pct)
            ) /
            2 - low_pct,
        )
        thresh_check = low_values[search_values.argmin()]

        if relaxed_thresh == thresh_check:
            return relaxed_thresh

        high_values = low_values[low_values > thresh_check]
        search_values = abs(
            high_values - max(high_values) +
            (max(high_values) - min(high_values)) / 2,
        )
        return high_values[search_values.argmin()]

    def genome(self):
        """
        Find the genome threshold.

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

        Parameters:
            vector (nparray): Vector for analysis

        Returns:
            float
        """
        model = ecdf(vector).cdf
        quantiles = 1 - model.evaluate(vector)
        return min(vector[quantiles < self.ctrl])

    def filter_by_pct(self, vector):
        """
        Filter pct_remain_vec results by pct_remain_vec.

        Parameters:
            vector (nparray): Vector to calculate by

        Returns:
            nparray
        """
        pct = pct_remain_vec(self.exp.vec, self.ctrl.vec, vector)
        return pct_remain_vec(self.exp.vec, self.ctrl.vec, vector[pct < 1])

    def thresh_check(self):
        """
        Check for spurious values and conditionally adjust thresholds.

        Returns:
            Adjusted thresholds (if they need adjusting)
        """
        sp_values = self.spurious_values()
        search_values = self.filter_by_pct(sp_values)
        old_values = self.filter_by_pct(self.both)
        if max(search_values) / max(old_values) > 0.95:
            sp_thresh = sp_values[search_values.argmax()]
            return sp_values[search_values.argmax()], self.stringent(sp_thresh)
        return None
