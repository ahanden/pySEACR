"""
Command line tool for pySEACR
"""
import argparse
import os

from pySEACR.auc_from_bdg import BDG
from pySEACR.normalize import Normalize
#from pySEACR.pct_remain import pct_remain_vec
from pySEACR.threshold_finder import ThresholdFinder


def parse_args():
    """
    Parse command line arguments.

    Returns:
        namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-e',
        '--exp',
    )
    parser.add_argument(
        '-c',
        '--ctrl',
    )
    parser.add_argument(
        '-n',
        '--norm',
    )
    parser.add_argument(
        '-H',
        '--height',
    )

    return parser.parse_args()

def main(args):
    """
    pySEACR main body.

    Parameters:
        args (namespace): Command line arguments
    """
    exp = BDG(args.exp)

    constant = 1
    if os.path.isfile(args.ctrl):
        ctrl = BDG(args.ctrl)
        if args.norm == "yes":
            norm = Normalize(exp, ctrl)
            constant = norm.constant()
            ctrl.vec = [_ * constant for _ in ctrl.vec]

        thresholds = ThresholdFinder(exp, ctrl)
        genome = thresholds.genome()
        relaxed = thresholds.relaxed()
        check = thresholds.thresh_check()

        if args.height == "relaxed":
            auc_thresh = relaxed if check is None else check[0]
            #fdr = 1 - pct_remain_vec(exp.vec, ctrl.vec, relaxed)
        else:
            stringent = thresholds.stringent(relaxed)
            auc_thresh = stringent if check is None else check[1]
            #fdr = 1 - pct_remain_vec(exp.vec, ctrl.vec, stringent)
    else:
        ctrl = float(args.ctrl)
        thresholds = ThresholdFinder(exp, ctrl)
        if args.height == "relaxed":
            auc_thresh = thresholds.static(exp.vec)
        else:
            auc_thresh = thresholds.static(exp.max)
        genome = 0
        #fdr = (ctrl, ctrl)

    write_output(exp.data, auc_thresh, genome)

def write_output(data, thresh_1, thresh_2):
    """
    Print pySEACR results in BED format.

    Parameters:
        data (list): AUC stretches from a BDG
        thresh_1 (float): Minimum peak height
        thresh_2 (float): Minimum peak width
    """
    for stretch in data:
        print(stretch.peak, stretch.n)
        if stretch.peak > thresh_1 and stretch.n > thresh_2:
            print("\t".join([str(_) for _ in [
                stretch.contig,
                stretch.start,
                stretch.stop,
                stretch.auc,
                stretch.peak,
                stretch.peak_coords(),
                stretch.n
            ]]))

if __name__ == "__main__":
    main(parse_args())
