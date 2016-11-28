class Hist3D:
    def __init__(self, edges, binsize):
        self.total_frames = 0
        (a, b, c) = edges
        # get bins
        self.Nbins = [math.floor(edge / binsize) for edge in edges]
        self.binSizes = [edge / N for (edge, N) in zip(edges, self.Nbins)]

        self.xedges = [num * self.binSizes[0] for num in range(self.Nbins[0] + 1)]  # coordinates of mesh in x-direction
        self.yedges = [num * self.binSizes[1] for num in range(self.Nbins[1] + 1)]
        self.zedges = [num * self.binSizes[2] for num in range(self.Nbins[2] + 1)]
        self.histogram = np.zeros((self.Nbins[0], self.Nbins[1], self.Nbins[2]), dtype=np.int8)
        if self.xedges[-1] != edges[0]:
            sys.stdout.write('something went wrong with initializing the bins %f != %f' % (self.xedges[-1], edges[0]));
            quit()

    def binCoords(self, Coords, a, b, c):
        '''
        note! this is not used, np.histogramdd does this for you!
        '''
        # test if coords in unit cell
        for (X, sidelength) in zip(Coords, [a, b, c]):
            if X > sidelength:
                sys.stdout.write('Coordinates outside unit cell, no binning will be done\n')
                quit()
        # Bin num starts at 0, so floor
        Bin = [math.floor(X / size) for (X, size) in zip(Coords, self.binSizes)]
        return Bin


class ZeoHist(Hist3D):
    def __init__(self, lattice_vectors, BINSIZE):
        Hist3D.__init__(self, lattice_vectors, BINSIZE)
        self.Vbin = self.binSizes[0] * self.binSizes[1] * self.binSizes[2]

    def Hist(self, xyzLocations):
        '''
        xyzLocations is a nested list ->  [[x1,y1,z1],[x2,y2,z2]...[xN,yN,zN]]
        '''
        Hist, edges = np.histogramdd(np.array(xyzLocations, dtype=np.int8),
                                     bins=(self.xedges, self.yedges, self.zedges), normed=False)
        return Hist, edges, np.sum(Hist)

    def addHist(self, coordinates):
        if coordinates:
            H, edges, total = self.Hist(coordinates)
            self.histogram = self.histogram + H

    def addZeoToHist(self, struc):
        """
        (method for debugging)
        :param struc: from CifParser(cif_file).get_structures(primitive=False)[0]
        :return: self.histogram now has -1000 in each bin occupied by zeolite
        """
        coords = [site['xyz'] for site in struc.as_dict()['sites']]
        Hist, edges = np.histogramdd(np.array(coords, dtype=np.int8),
                                     bins=(self.xedges, self.yedges, self.zedges), normed=False)
        self.histogram = self.histogram + -1000*Hist

    def normalize(self):
        #        self.histogram /= self.totalCount no! dont do this for numpy objects!!
        self.histogram = self.histogram / (np.sum(self.histogram))


class HBondHist(ZeoHist):
    def __init__(self, lattice_vectors, BINSIZE):
        ZeoHist.__init__(self, lattice_vectors, BINSIZE)
        # note, I don't think we should use integer types in this historam, since we divide by the volume
        self.Hsingle_hist = np.array(self.histogram, copy=True, dtype=float)
        self.Hmult_hist = np.array(self.histogram, copy=True, dtype=float)
        self.Hbond_locations = []

    def addHist(self, FRAME_DATA, boxdimensions, mainMlcl, main_bead, box,
                ZEO_DATA=False, Hself=False, Hother=False, Hzeo=False):
        # this is typically called from hbondlocation.py
        # find buoh hbonded with zeo

        frame_locations, *[mult_locations] = Hbond.locationhbond(FRAME_DATA, boxdimensions, mainMlcl, main_bead, box,
                                                               ZEO_DATA=ZEO_DATA, Hself=Hself, Hother=Hother, Hzeo=Hzeo)
        # frame_locations is nested list -> [[x1,y1,z1],[x2,y2,z2]...[xN,yN,zN]]
        if frame_locations:
            H1, H1edges, total = self.Hist(frame_locations)
            self.Hsingle_hist = self.Hsingle_hist + H1
            if set(self.xedges) != set(H1edges[0]) or set(self.yedges) != set(H1edges[1]) or set(self.zedges) != set(
                    H1edges[2]):
                sys.stdout.write('nphistogram is changing the bin edges with out my consent');
                quit()
            if mult_locations:
                HM, HM_edges, total = self.Hist(mult_locations)
                self.Hmult_hist = self.Hmult_hist + HM
        for Coords in frame_locations:
            self.Hbond_locations.append(Coords)

    def normalize(self):
        # to get number density...
        #       self.Hsingle_hist /= (self.total_frames*self.Vbin)
        #       self.Hmult_hist /= (self.total_frames*self.Vbin)
        # else to get probability of finding ...
        self.Hsingle_hist = self.Hsingle_hist / np.sum(self.Hsingle_hist)
        self.Hmult_hist = self.Hmult_hist / np.sum(self.Hmult_hist)


import numpy as np
import sys, math
import StructureAnalysis.Hbond as Hbond
