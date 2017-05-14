from MCFlow.file_formatting.reader import Movie
from MCFlow.structure_analysis import Struc

def writeXY(x,y,name):
    with open(name,'w') as f:
        for i,j in zip(x,y):
            f.write('%e %e\n'%(i,j))

class RDF(Movie):

    def __init__(self, file_name, *args):
        Movie.__init__(self, file_name, *args)
        self.rMax, self.dr, box = args
        if 'box' in box:
            self.box = box
        else:
            self.box = 'box%s'%box
        self.g_bead = []
        self.edges = np.arange(0., self.rMax + 1.1*self.dr, self.dr)
        self.num_bins = len(self.edges) - 1

    def averageG(self):
        self.radii = np.zeros(self.num_bins)
        # convert self.g (nested list) to g (np.matrix)
        g = np.zeros([len(self.g_bead),self.num_bins])
        for i, val in enumerate(self.g_bead):
            g[i,:] = val
        self.g_average = np.zeros(len(self.edges)-1)
        for i in range(len(self.edges)-1):
            rOuter = self.edges[i+1]
            rInner = self.edges[i]
            self.radii[i] = (rInner+rOuter)/2.
            self.g_average[i] = np.mean(g[:,i]) / (4.0/3.0*np.pi*(rOuter**3 - rInner**3))

    def g_frame(self, coords1, coords2, boxLengths ):
        '''
        compute the three-dimensional pair correlsation function for a set of spherical particles contained
        in a rectangular prism
        :param coords1: a matrix of xyz coordinates of bead 1 [[x1,y1,z1]...[xN,yN,zN]]
        :param coords2: a matrix of xyz coordinates of bead 2 [[x1,y1,z1]...[xN,yN,zN]]
        :param S:  a matrix of boxlengths in each direction
        :param dr:  increment for increasing radius of spherical shell
        :return: returns a tuple: (g, radii, interior_indices)
            g(r)        a numpy array containing the correlation function g(r)
            radii       a numpy array containing the radii of the
                        spherical shells used to compute g(r)
            reference_indices   indices of reference particles
        '''
        def minImage(dist, boxl):
            return (
                dist - np.round( np.divide(dist,boxl) )*boxl
            )
        n_coords1, n_coords2 = len(coords1[:,0]) , len(coords2[:,0])
        assert n_coords1 > 3, 'Wrong dimension {}'.format(n_coords1)
        num_particles = n_coords1 + n_coords2
        number_density = num_particles/(boxLengths[0]*boxLengths[1]*boxLengths[2])

        for index in range(n_coords1): #TODO: figure out if this summation is correct
            dx = coords1[index,0] - coords2[:,0]
            dy = coords1[index,1] - coords2[:,1]
            dz = coords1[index,2] - coords2[:,2]
            dx = minImage(dx, boxLengths[0])
            dy = minImage(dy, boxLengths[1])
            dz = minImage(dz, boxLengths[2])
            dist = np.sqrt(dx**2 + dy**2 + dz**2)
            (result, bins) = np.histogram(dist, bins=self.edges, normed=False)
            self.g_bead.append( result / number_density )

    def calculateRDF(self, m1, b1, m2, b2):
        for iframe, FRAME_DATA in enumerate(self.frame_data):
            c1_xyz = []
            c2_xyz = []
            for i, mol1 in enumerate(FRAME_DATA[self.box]['mol%s'%m1]):
                for each_bead in mol1[b1]:
                    c1_xyz.append( list(map(float,each_bead ) ) )
            for i, mol2 in enumerate(FRAME_DATA[self.box]['mol%s'%m2]):
                for each_bead in mol2[b2]:
                    c2_xyz.append( list(map(float,each_bead ) ) )
            c1_xyz = np.array(c1_xyz)
            c2_xyz = np.array(c2_xyz)
            self.g_frame( c1_xyz, c2_xyz, self.boxlengths[iframe][self.box])
        self.averageG()

class G(Struc):
    def __init__(self):
        my_parser = Results()
        my_parser.parser.add_argument('-B1','--bead1',help='bead from for rdf',type = str)
        my_parser.parser.add_argument('-M1','--mol1',help='mol from for rdf',type = str)
        my_parser.parser.add_argument('-B2','--bead2',help='bead to for rdf',type = str)
        my_parser.parser.add_argument('-M2','--mol2',help='mol to for rdf',type = str)
        my_parser.parser.add_argument('-bins','--bins',help='radial bin size for rdf',type = float)
        my_parser.parser.add_argument('-r','--rmax',help='max radius for rdf',type = float)
        args = vars(my_parser.parse_args())
        self.args = args
        self.checks()

    def checks(self):
        assert (self.args['bead1'] and self.args['bead2'] and self.args['mol1']
                 and self.args['mol2'] and self.args['bins'] and self.args['rmax']), 'Necessary input missing'
        assert self.args['box'], 'box needed for rdf'
        self.analysis_class = RDF

    def myCalcs(self, D ):
        D.calculateRDF(self.args['mol1'],self.args['bead1'],
                        self.args['mol2'],self.args['bead2'])
        writeXY(D.radii, D.g_average, '%s/%s/rdf-box%s_mol%s-%s_mol%s-%s.dat'%(self.args['path'],self.feed,self.args['box'],
                                            self.args['mol1'], self.args['bead1'], self.args['mol2'], self.args['bead2']))

    def main(self):
        for feed in self.args['feeds']:
            self.feed = feed
            if self.args['verbosity'] > 0: print('-'*12 + 'Dir is %s'%self.feed + '-'*12)
            analysis = self.read_movies(self.args['rmax'], self.args['bins'], self.args['box'])
            self.myCalcs(analysis)

import numpy as np
from MCFlow.parser import Results

if __name__ == '__main__':
    M = G()
    M.main()
