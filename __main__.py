from pySEACR.threshold_finder import ThresholdFinder
from pySEACR.pct_remain import pct_remain_vec
from pySEACR.auc_from_bdg import BDG
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
    print(f"Reading {args.exp}")
    exp = BDG(args.exp)
    print(f"Read {len(exp.vec)} auc regions")

    constant = 1
    if os.path.isfile(args.ctrl):
        print(f"Reading {args.ctrl}")
        ctrl = BDG(args.ctrl)
        print(f"Read {len(ctrl.vec)} auc regions")
        if args.norm == "yes":
            print("Normalizing")
            constant = find_constant(exp.vec, ctrl.vec)
            ctrl.vec = [_ * constant for _ in ctrl.vec]

        thresholds = ThresholdFinder(exp, ctrl)
        print("Finding relaxed threshold")
        relaxed = thresholds.relaxed()
        print("Finding stringent threshold")
        stringent = thresholds.stringent(relaxed)

        print("Checking outliers")
        check = thresholds.thresh_check()
        if check is not None:
            relaxed, stringent = check

        print("Finding genome threshold")
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
    print("Writing output")
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
