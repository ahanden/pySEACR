import sys
import csv
import numpy as np

class BDGRow(object):
    def __init__(self, data):
        self.contig = data[0]
        self.start = int(data[1])
        self.stop = int(data[2])
        self.peak = float(data[3])

class Stretch(object):
    def __init__(self, row):
        self.contig, self.start, self.stop = row.contig, row.start, row.stop
        self.peak_start, self.peak_stop = row.start, row.stop
        self.peak = row.peak
        self._auc = self._calc_auc(row.start, row.stop, row.peak)
        self.n = 1

    def _calc_auc(self, start, stop, peak):
        return peak * (stop - start)

    def extend(self, row):
        self.n += 1
        self.stop = row.stop
        self._auc += self._calc_auc(row.start, row.stop, row.peak)
        if row.peak > self.peak:
            self.peak = row.peak
            self.peak_start, self.peak_stop = row.start, row.stop
        elif row.peak == self.peak:
            self.peak_stop = row.stop

    def coords(self):
        return f"{self.contig}:{self.start}-{self.stop}"

    def peak_coords(self):
        return f"{self.contig}:{self.peak_start}-{self.peak_stop}"

    @property
    def auc(self):
        return round(self._auc * 100) / 100

class BDG(object):
    def __init__(self, bdg_fname):
        self.data = []
        regions_gen = self.read_bdg(bdg_fname)
        region = next(regions_gen)
        while region.peak == 0:
            region = next(regions_gen)
        auc_stretch = Stretch(region)
        for region in regions_gen:
            if region.peak == 0:
                continue
            if auc_stretch.contig == region.contig and auc_stretch.stop == region.start:
                auc_stretch.extend(region)
            else:
                self.data.append(auc_stretch)
                auc_stretch = Stretch(region)
        self.data.append(auc_stretch)
        self.vec = np.array([_.auc for _ in self.data])
        self.max = np.array([_.peak for _ in self.data])

    def read_bdg(self, bdg_fname):
        with open(bdg_fname, "r") as stream:
            reader = csv.reader(stream, delimiter="\t")
            _ = next(reader)
            for row in reader:
                yield BDGRow(row)

