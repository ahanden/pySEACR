from pySEACR.auc_bed_file import AUCBEDFile
from pySEACR.threshold_finder import ThresholdFinder
from pySEACR.pct_remain import pct_remain_vec
import argparse
import os

def parse_args():
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
        '-o',
        '--output',
    )

    return parser.parse_args()

def main(args):
    exp = AUCBEDFile(args.exp)

    constant = 1
    if os.path.isfile(args.ctrl):
        ctrl = AUCBEDFile(args.ctrl)
        if args.norm == "yes":
            constant = find_constant(exp.vec, ctrl.vec)
            ctrl.vec = [_ * constant for _ in ctrl.vec]

        thresholds = ThresholdFinder(exp, ctrl)
        relaxed = thresholds.relaxed()
        stringent = thresholds.stringent(relaxed)

        check = thresholds.thresh_check()
        if check is not None:
            relaxed, stringent = check

        genome = thresholds.genome()
        fdr = [
            1 - pct_remain_vec(exp.vec, ctrl.vec, relaxed),
            1 - pct_remain_vec(exp.vec, ctrl.vec, stringent),
        ]
    else:
        ctrl = float(args.ctrl)
        thresholds = ThresholdFinder(exp, ctrl)
        relaxed = thresholds.static(exp.vec)
        stringent = thresholds.static(exp.max)
        genome = 0
        fdr = (ctrl, ctrl)

    with open(f'{args.output}_var.sh', 'w') as stream:
        stream.write(f'''
SEACR_THRESH_R={relaxed}
SEACR_THRESH_S={stringent}
SEACR_THRESH_G={genome}
SEACR_NORM={constant}
SEACR_FDR_S={fdr[0]}
SEACR_FDR_R={fdr[1]}
''')

if __name__ == "__main__":
    args = parse_args()
    main(args)
