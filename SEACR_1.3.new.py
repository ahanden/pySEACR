from read_auc_bed import read_auc_bed
from new_seacr import *
import argparse
import os

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

args = parser.parse_args()

exp = read_auc_bed(args.exp)

constant = 1
if os.path.isfile(args.ctrl):
    ctrl = read_auc_bed(args.ctrl)
    if args.norm == "yes":
        constant = find_constant(exp["vec"], ctrl["vec"])
        ctrl["vec"] = [_ * constant for _ in ctrl["vec"]]

    relaxed = find_relaxed_thresh(exp["vec"], ctrl["vec"])
    stringent = find_stringent_thresh(exp["vec"], ctrl["vec"], relaxed)

    check = thresh_check(exp["vec"], ctrl["vec"])
    if check is not None:
        relaxed, stringent = check

    genome = find_global_thresh(exp["max"], ctrl["max"])
    fdr = [
        1 - pct_remain_vec(exp["vec"], ctrl["vec"], relaxed),
        1 - pct_remain_vec(exp["vec"], ctrl["vec"], stringent),
    ]
else:
    ctrl = float(args.ctrl)
    thresholds = get_static_thresholds(exp["vec"], exp["max"], ctrl)
    relaxed, stringent, genome, fdr = thresholds

with open(f'{args.output}_var.sh', 'w') as stream:
    stream.write(f'''
SEACR_THRESH_R={relaxed}
SEACR_THRESH_S={stringent}
SEACR_THRESH_G={genome}
SEACR_NORM={constant}
SEACR_FDR_S={fdr[0]}
SEACR_FDR_R={fdr[1]}
''')
