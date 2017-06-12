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
        self.valToIgnore = -1


    def Smap(self, xyzA, xyzB, Kref):
        '''

        :param nFrames: number of frames in data
        :param xyzA: nested list for bead A ->  [[x1,y1,z1],[x2,y2,z2]...[xN,yN,zN]]
        :param xyzB: nested list for bead B ->  [[x1,y1,z1],[x2,y2,z2]...[xN,yN,zN]]
        :param Kref: reference K (xA/xB) in reference phase
        :return:
        '''
        if not isinstance(xyzA,np.ndarray):
            xyzA = np.array(xyzA, dtype=np.int8)
        if not isinstance(xyzB,np.ndarray):
            xyzB = np.array(xyzB, dtype=np.int8)

        HistA, edgesA = np.histogramdd(xyzA,
                                       bins=(self.xedges, self.yedges, self.zedges), normed=False)
        HistB, edgesB = np.histogramdd(xyzB,
                                       bins=(self.xedges, self.yedges, self.zedges), normed=False)
        # might need to always check here that edges are the same
        if Kref == None:
            ref = 1.
        else:
            ref = Kref

        self.edges = edgesA
        self.histogram = np.ones(np.shape(HistA))*self.valToIgnore
        sum_hist = 0
        for x in range(len(self.histogram[:,0,0])):
            for y in range(len(self.histogram[0,:,0])):
                for z in range(len(self.histogram[0,0,:])):
                    sum_hist = sum_hist + HistB[x,y,z]
                    if HistB[x,y,z] > 0:
                        self.histogram[x,y,z] = (HistA[x,y,z]/HistB[x,y,z]) / ref
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
    my_parser.parser.add_argument('-off','--offAxis',help='axis to average along',
                                  default=[],nargs='+')
    args = vars(my_parser.parse_args())
    axis_conf = {'x':0,'y':1,'z':2}

    coordsA = reader.xyz(args['file'])
    xyz_dataA = getCoords(coordsA, args['bead'])
    coordsB = reader.xyz(args['file2'])
    xyz_dataB = getCoords(coordsB, args['bead'])
    for axis in args['offAxis']:
        xyz_dataA = np.array(xyz_dataA, dtype=np.int8)
        xyz_dataB = np.array(xyz_dataB, dtype=np.int8)
        xyz_dataA[:, axis_conf[axis]] = 0.
        xyz_dataB[:, axis_conf[axis]] = 0.

    hist = SHist(args['vectors'],args['bins'])
    hist.Smap(xyz_dataA, xyz_dataB, args['reference'])

    new_file = getFileName(args['file']) + getFileName(args['file2']) + '_bead%s.vtk'%('-'.join(args['bead']))
    vtkRectilinearMesh(new_file, hist.edges, hist.histogram)
