"""Calculating area under the curve (AUC) of a BedGraph file."""
import csv

import numpy as np


class BDGRow():
    """Representation of a row from a BedGraph file."""

    def __init__(self, row):
        """
        Create a new BDGRow.

        Parameters:
            row (list): Fields from a BedGraph file.
        """
        self.contig = row[0]
        self.start = int(row[1])
        self.stop = int(row[2])
        self.peak = float(row[3])


class Stretch():
    """SEACR identified window of regions to "clump"."""

    def __init__(self, row):
        """
        Create a new Stretch.

        Parameters:
            row (BDGRow): The first element of the stretch
        """
        self.contig = row.contig
        self.coord = (row.start, row.stop)
        self.peak_coord = (row.start, row.stop)
        self.peak = row.peak
        self._auc = self._calc_auc(row.start, row.stop, row.peak)
        self.n = 1

    def is_contiguous(self, row):
        """
        Determine if a region is contiguous with this stretch.

        Parameters:
            row (BDGRow): Region to check

        Returns:
            True if contigous, else False
        """
        return self.contig == row.contig and self.coord[1] == row.start

    def extend(self, row):
        """
        Extend the stretch by adding another BDG entry.

        Parameters:
            row (BDGRow): The region to add.
        """
        self.n += 1
        self.coord = (self.coord[0], row.stop)
        self._auc += self._calc_auc(row.start, row.stop, row.peak)
        if row.peak > self.peak:
            self.peak = row.peak
            self.peak_coord = (row.start, row.stop)
        elif row.peak == self.peak:
            self.peak_coord = (self.peak_coord[0], row.stop)

    def peak_coords(self):
        """
        String representation of the highest point in the stretch.

        Returns:
            String in the format of contig:peak_start-peak_stop
        """
        return f'{self.contig}:{self.peak_coord[0]}-{self.peak_coord[1]}'

    @property
    def auc(self):
        """
        Area under the curve.

        Returns:
            Float rounded to 2 digits.
        """
        return round(self._auc * 100) / 100

    def _calc_auc(self, start, stop, peak):
        return peak * (stop - start)


class BDG():
    """BedGraph File."""

    def __init__(self, bdg_fname):
        """
        Read a BedGraph file into memory.

        Parameters:
            bdg_fname (str): File name
        """
        self.regions = []
        regions_gen = self._read_bdg(bdg_fname)
        start = next(regions_gen)
        while start.peak == 0:
            start = next(regions_gen)
        auc_stretch = Stretch(start)
        for region in regions_gen:
            if region.peak == 0:
                continue
            if auc_stretch.is_contiguous(region):
                auc_stretch.extend(region)
            else:
                self.regions.append(auc_stretch)
                auc_stretch = Stretch(region)
        self.regions.append(auc_stretch)
        self.vec = np.array([_.auc for _ in self.regions])
        self.max = np.array([_.n for _ in self.regions])

    def _read_bdg(self, bdg_fname):
        with open(bdg_fname, 'r') as stream:
            reader = csv.reader(stream, delimiter='\t')
            next(reader)
            for row in reader:
                yield BDGRow(row)
