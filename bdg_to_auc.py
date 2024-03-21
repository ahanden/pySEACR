import sys
import csv 


def calc_auc(region):
    return region["peak"] * (region["stop"] - region["start"])

def frmt_coord(contig, coord):
    return f"{contig}:{coord[0]}-{coord[1]}"

def region_to_str(region):
    return "\t".join([str(_) for _ in [
        region["contig"],
        region["start"],
        region["stop"],
        f"{region['auc']:.2f}",
        region["peak"],
        frmt_coord(region["contig"], region["peak_coord"]),
        region["n"],
    ]])


def parse_bdg_row(row):
    return {
        "contig": row[0],
        "start": int(row[1]),
        "stop": int(row[2]),
        "peak": float(row[3]),
    }

def setup_auc_stretch(region):
    region["peak_coord"] = [region["start"], region["stop"]]
    region["auc"] = calc_auc(region)
    region["n"] = 1
    return region

def extend_stretch(source, to_add):
    source["n"] += 1
    source["stop"] = to_add["stop"]
    source["auc"] += calc_auc(to_add)
    if to_add["peak"] > source["peak"]:
        source["peak"] = to_add["peak"]
        source["peak_coord"] = [to_add["start"], to_add["stop"]]
    elif to_add["peak"] == source["peak"]:
        source["peak_coord"][1] = to_add["stop"]

def bdg_to_auc(bdg_fname):
    with open(bdg_fname, "r") as stream:
        reader = csv.reader(stream, delimiter="\t")
        _ = next(reader)
        row = next(reader)
        while(row[3] == "0"):
            row = next(reader)

        auc_stretch = setup_auc_stretch(parse_bdg_row(row))

        for row in reader:
            region = parse_bdg_row(row)
            if region["peak"] == 0:
                continue
            if auc_stretch["contig"] == region["contig"] and auc_stretch["stop"] == region["start"]:
                extend_stretch(auc_stretch, region)
            else:
                print(region_to_str(auc_stretch))
                auc_stretch = setup_auc_stretch(region)
        print(region_to_str(auc_stretch))

if __name__ == "__main__":
    bdg_to_auc(sys.argv[1])
