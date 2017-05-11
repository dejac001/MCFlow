class Hist3D:
    def __init__(self, edges, binsize):
        # get bins
        self.Nbins = [math.floor(edge / binsize) for edge in edges]
        self.binSizes = [edge / N for (edge, N) in zip(edges, self.Nbins)]
        self.Vbin = self.binSizes[0]*self.binSizes[1]*self.binSizes[2]
        # coordinates of mesh in each direction
        self.xedges = [num * self.binSizes[0] for num in range(self.Nbins[0] + 1)]
        self.yedges = [num * self.binSizes[1] for num in range(self.Nbins[1] + 1)]
        self.zedges = [num * self.binSizes[2] for num in range(self.Nbins[2] + 1)]
        self.name = 'hist'
        self.valueToIgnore = None
        print('Total number of bins is %i'%(self.Nbins[0]*self.Nbins[1]*self.Nbins[2]))

    def makeHist(self, xyzLocations):
        '''
        xyzLocations is a nested list ->  [[x1,y1,z1],[x2,y2,z2]...[xN,yN,zN]]
        '''
        Hist, edges = np.histogramdd(np.array(xyzLocations, dtype=np.int8),
                                     bins=(self.xedges, self.yedges, self.zedges), normed=False)
        self.histogram = Hist / np.sum(Hist)
        self.edges = edges

    def colorValues(self):
        values = []
        for x in range(len(self.histogram[:,0,0])):
            for y in range(len(self.histogram[0,:,0])):
                for z in range(len(self.histogram[0,0,:])):
                    if self.histogram[x,y,z] != self.valToIgnore:
                        values.append( self.histogram[x,y,z] )
        self.color_values = values

class dG3D(Hist3D):
    def __init__(self, edges, binsize):
        Hist3D.__init__(self,edges,binsize)
        self.name = 'dG'
        self.valueToIgnore = 5000

    def dGmap(self, nFrames, xyzLocations, densFrom, Temp):
        '''
        xyzLocations is a nested list ->  [[x1,y1,z1],[x2,y2,z2]...[xN,yN,zN]]
        densFrom is number density from in units molec/nm**3
        '''
        Hist, edges = np.histogramdd(np.array(xyzLocations, dtype=np.int8),
                                     bins=(self.xedges, self.yedges, self.zedges), normed=False)
        self.edges = edges
        self.histogram = np.ones(np.shape(Hist))
        self.histogram = self.histogram*5000
        for x in range(len(self.histogram[:,0,0])):
            for y in range(len(self.histogram[0,:,0])):
                for z in range(len(self.histogram[0,0,:])):
                    if Hist[x,y,z] > 0:
                        self.histogram[x,y,z] = -8.314/1000*Temp*np.log( ((Hist[x,y,z]/nFrames)/self.Vbin)
                                                                            / densFrom )



def getCoords(xyz_data, bead):
    data = []
    for atom, xyz in zip(xyz_data['atoms'], xyz_data['coords']):
        if atom == bead: data.append(xyz)
    return data

def getFileName(file):
    file = file[(file.find('/')+1):]
    return file.rstrip('.xyz')

import numpy as np
import math
from MCFlow.parser import Structure
from MCFlow.file_formatting import reader
from MCFlow.file_formatting.writer import vtkRectilinearMesh

if __name__ == '__main__':
    my_parser = Structure()

    my_parser.analysis()
    args = vars(my_parser.parse_args())

    coords = reader.xyz(args['file'])
    xyz_data = getCoords(coords, args['bead'])

    if not args['reference']:
        hist = Hist3D(args['vectors'],args['bins'])
        hist.makeHist(xyz_data)
    else:
        hist = dG3D(args['vectors'],args['bins'])
        hist.dGmap(args['numFrames'],xyz_data, args['reference'], args['Temp'])

    new_file = getFileName(args['file']) + '_bead%s_%s.vtk'%(args['bead'],hist.name)
    vtkRectilinearMesh(new_file, hist.edges, hist.histogram)