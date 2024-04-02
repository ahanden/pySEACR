"""Command line tool for pySEACR."""
import argparse
import os

from pySEACR.auc_from_bdg import BDG
from pySEACR.normalize import Normalize
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


def normalize(exp, ctrl):
    """
    Normalize the values in the control data.

    Parameters:
        exp (BDG): BDG file object of experimental/treatment data
        ctrl (BDG): BDG file object of the contorl/IgG data

    Returns:
        Normalized numpy.array
    """
    normer = Normalize(exp, ctrl)
    return ctrl.vec * normer.constant()


def thresh_from_ctrl(exp, ctrl_file, norm, height):
    """
    Compute height thresholds by comparing two files.

    Parameters:
        exp (BDG): Experimental/treatment BDG file object
        ctrl_file (str): Path to the control/IgG bdg file
        norm (str): If yes, will normalize data before calculating thresholds
        height (str): relaxed or stringent

    Returns:
        height threshold and genome threshold
    """
    ctrl = BDG(ctrl_file)
    if norm == 'yes':
        ctrl.vec = normalize(exp, ctrl)

    thresholds = ThresholdFinder(exp, ctrl)
    auc_thresh = thresholds.relaxed()
    check = thresholds.thresh_check()

    if height == 'relaxed':
        if check is not None:
            auc_thresh = check[0]
    else:
        stringent = thresholds.stringent(auc_thresh)
        auc_thresh = stringent if check is None else check[1]

    return auc_thresh, thresholds.genome()


def main(args):
    """
    pySEACR main body.

    Parameters:
        args (namespace): Command line arguments
    """
    exp = BDG(args.exp)

    auc_thresh = 0
    genome = 0
    if os.path.isfile(args.ctrl):
        auc_thresh, genome = thresh_from_ctrl(
            exp,
            args.ctrl,
            args.norm,
            args.height,
        )
    else:
        args.ctrl = float(args.ctrl)
        thresholds = ThresholdFinder(exp, args.ctrl)
        vec_input = exp.vec if args.height == 'relaxed' else exp.max
        auc_thresh = thresholds.static(vec_input)

    print_output(exp.regions, auc_thresh, genome)


def print_output(stretches, height, width):
    """
    Print pySEACR results in BED format.

    Parameters:
        stretches (list): AUC stretches from a BDG
        height (float): Minimum peak height
        width (float): Minimum peak width
    """
    for stretch in stretches:
        if stretch.peak > height and stretch.n > width:
            print('\t'.join([str(_) for _ in (
                stretch.contig,
                stretch.coord[0],
                stretch.coord[1],
                stretch.auc,
                stretch.peak,
                stretch.peak_coords(),
                stretch.n,
            )]))


if __name__ == '__main__':
    main(parse_args())
