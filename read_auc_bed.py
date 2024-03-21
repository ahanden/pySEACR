import csv
def read_auc_bed(*args):
    file_name = args[0]
    with open(file_name, "r") as stream:
        reader = csv.reader(stream, delimiter="\t")
        vec = []
        max = []
        for row in reader:
            vec.append(float(row[3]))
            max.append(float(row[4]))
    return {"vec": vec, "max": max}
