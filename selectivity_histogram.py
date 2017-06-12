from MCFlow.probhistogram import Hist3D

def hist_norm_height(n,bins,const):
    ''' Function to normalise bin height by a constant.
        Needs n and bins from np.histogram or ax.hist.'''

    n = np.repeat(n,2)
    n = np.float32(n) / const
    new_bins = [bins[0]]
    new_bins.extend(np.repeat(bins[1:],2))
    return n,new_bins[:-1]

class SHist(Hist3D):

    def __init__(self, edges, binsize):
        Hist3D.__init__(self, edges, binsize)
        self.name='S'
        self.valueToIgnore = -1


    def Smap(self, xyzA, xyzB, Kref):
        '''

        :param nFrames: number of frames in data
        :param xyzA: nested list for bead A ->  [[x1,y1,z1],[x2,y2,z2]...[xN,yN,zN]]
        :param xyzB: nested list for bead B ->  [[x1,y1,z1],[x2,y2,z2]...[xN,yN,zN]]
        :param Kref: reference K (xA/xB) in reference phase
        :return:
        '''
        HistA, edgesA = np.histogramdd(np.array(xyzA, dtype=np.int8),
                                       bins=(self.xedges, self.yedges, self.zedges), normed=False)
        HistB, edgesB = np.histogramdd(np.array(xyzB, dtype=np.int8),
                                       bins=(self.xedges, self.yedges, self.zedges), normed=False)
        # might need to always check here that edges are the same
        if Kref == None:
            ref = 1.
        else:
            ref = Kref

        self.edges = edgesA
        self.histogram = np.ones(np.shape(HistA))
        self.histogram = -1*self.histogram
        sum_hist = 0
        for x in range(len(self.histogram[:,0,0])):
            for y in range(len(self.histogram[0,:,0])):
                for z in range(len(self.histogram[0,0,:])):
                    sum_hist = sum_hist + HistB[x,y,z]
                    if HistB[x,y,z] > 0:
                        self.histogram[x,y,z] = (HistA[x,y,z]/HistB[x,y,z]) / ref
                    else:
                        self.histogram[x,y,z] = -1
        assert sum_hist == len(xyzB), 'Histogram sum incorrect; %i %i'%(sum_hist, len(xyzB))



import math
import numpy as np
import matplotlib.pyplot as plt
from MCFlow.parser import Structure
from MCFlow.file_formatting import reader
from MCFlow.file_formatting.writer import vtkRectilinearMesh
from MCFlow.probhistogram import getCoords, getFileName

if __name__ == '__main__':
    my_parser = Structure()
    my_parser.parser.add_argument('-f2','--file2',help='file 2 for selectivity analysis',type=str)
    my_parser.parser.add_argument('-ref','--reference',help='reference density  [for S, Kref]',type=float)
    args = vars(my_parser.parse_args())

    coordsA = reader.xyz(args['file'])
    xyz_dataA = getCoords(coordsA, args['bead'])
    coordsB = reader.xyz(args['file2'])
    xyz_dataB = getCoords(coordsB, args['bead'])

    hist = SHist(args['vectors'],args['bins'])
    hist.Smap(xyz_dataA, xyz_dataB, args['reference'])

    new_file = getFileName(args['file']) + getFileName(args['file2']) + '_bead%s.vtk'%('-'.join(args['bead']))
    vtkRectilinearMesh(new_file, hist.edges, hist.histogram)
