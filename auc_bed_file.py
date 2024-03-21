import csv

class AUCBEDFile(object):
    def __init__(self, file_name):
        with open(file_name, "r") as stream:
            reader = csv.reader(stream, delimiter="\t")
            self.vec = []
            self.max = []
            for row in reader:
                self.vec.append(float(row[3]))
                self.max.append(float(row[4]))
