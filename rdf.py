from MCFlow.file_formatting.reader import Movie
from MCFlow.structure_analysis import Struc

def writeXY(x,y,name):
    with open(name,'w') as f:
        for i,j in zip(x,y):
            f.write('%e %e\n'%(i,j))

class RDF(Movie):

    def __init__(self, file_name, *args):
        Movie.__init__(self, file_name, *args)
        self.rMax, self.dr, dimensions = args
        self.dimension_indices = [['x','y','z'].index(i) for i in dimensions]
        assert len(self.dimension_indices) != 2, '2D RDF not implemented yet'
        self.g_bead = []
        if len(self.dimension_indices) == 3:
            self.edges = np.arange(0., self.rMax + self.dr, self.dr)
        else:
            self.edges = np.arange(-self.rMax, self.rMax + self.dr, self.dr)
        self.boxLengths = []
        self.num_bins = len(self.edges) - 1
        self.numInt = []
        self.unit_cells = [2,2,3] # MFI
        self.sameMol = False

    def getAllBoxlx(self):
        tolerance = 1e-8
        boxlengths = self.boxLengths[0]
        for lengths in self.boxLengths[1:]:
            for i in range(len(lengths)):
                difference = abs(lengths[i]-boxlengths[i])
                assert difference < tolerance, 'Box lengths changed by %e'%difference
        boxlengths_ignore = [boxlengths[i] for i in range(len(boxlengths))
                             if i not in self.dimension_indices]
        uc_ignore = [self.unit_cells[i] for i in range(len(self.unit_cells))
                     if i not in self.dimension_indices]
        self.bin_area = boxlengths_ignore[0]*boxlengths_ignore[1]
        self.n_multiples = uc_ignore[0]*uc_ignore[1]

    def averageG(self):
        if len(self.dimension_indices) == 1:
            self.getAllBoxlx()
        self.radii = np.zeros(self.num_bins)
        # convert self.g (nested list) to g (np.matrix)
        g = np.zeros([len(self.g_bead),self.num_bins])
        n = np.zeros([len(self.numInt), self.num_bins])
        for i in range(len(self.g_bead)):
            g[i,:] = self.g_bead[i]
            n[i,:] = self.numInt[i]
        self.g_average = np.zeros(len(self.edges)-1)
        self.n_average = np.zeros(len(self.edges)-1)
        for i in range(len(self.edges)-1):
            rOuter = self.edges[i+1]
            rInner = self.edges[i]
            self.radii[i] = (rInner+rOuter)/2.
            if len(self.dimension_indices) == 1:
                dV = self.dr
                self.n_average[i] = np.mean(n[:,i])/self.n_multiples
            elif len(self.dimension_indices) == 3:
                self.n_average[i] = np.mean(n[:,i])
                dV = 4./3.*np.pi*(rOuter**3 - rInner**3)
            self.g_average[i] = np.mean(g[:,i]) / dV
        self.g_bead = []
        self.numInt = []

    def g_frame(self, coords1, coords2, boxLengths ):
        '''
        compute the three-dimensional pair correlation function for a set of spherical particles contained
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
        volume = 1.
        for i in self.dimension_indices:
            volume = volume*boxLengths[i]
        if self.sameMol:
            assert n_coords1 == n_coords2, 'Inconsistency in number of atoms'
            n_tot = (n_coords1*(n_coords2-1))/2
        else:
            n_tot = (n_coords1*n_coords2)
        constant = n_tot/volume

        for index in range(n_coords1):
            if len(self.dimension_indices) > 1:
                sum_squared = 0.
                for i in self.dimension_indices:
                    my_dist = coords1[index,i] - coords2[:,i]
                    my_dist = minImage(my_dist, boxLengths[i])
                    sum_squared = sum_squared + np.square(my_dist)
                dist = np.sqrt(sum_squared)
            else:
                i = self.dimension_indices[0]
                dist = coords1[index,i] - coords2[:,i]
                dist = minImage(dist, boxLengths[i])
            if (index < len(dist)) and (abs(dist[index]) < 0.01): dist[index] = 1000.
            (result, bins) = np.histogram(dist, bins=self.edges, normed=False)
            self.g_bead.append( result / constant )
            if abs(self.edges[0]) < 1e-8:
                self.numInt.append( [ sum(result[:i]) for i in range(len(result))] )
            else:
                zero = np.where(abs(self.edges) < 1e-8)[0]
                assert len(zero) == 1, 'More than 1 zero'
                zero = zero[0]
                nint_positive = [ sum(result[zero:i]) for i in range(zero+1,len(result)) ]
                nint_negative = [ sum(result[i:zero]) for i in range(zero) ]
                self.numInt.append( nint_negative + [0.] + nint_positive )

    def calculateRDF(self, m1, b1, m2, b2, sbox):
        self.box = 'box%s'%sbox
        mol1, mol2 = ['mol%s'%m1,'mol%s'%m2]
        if mol1 == mol2:
            self.sameMol = True
        self.ignored_frames = 0
        for iframe, FRAME_DATA in enumerate(self.frame_data):
            try:
                if (iframe+1)%(self.nframes//4) == 0:
                    print('%5.1f %% of RDF frames completed'%(100*(iframe+1)/self.nframes))
            except ZeroDivisionError:
                print('%5.1f %% of RDF frames completed'%(100*(iframe+1)/self.nframes))
            c1_xyz = []
            c2_xyz = []
            for i, each_mol in enumerate(FRAME_DATA[self.box][mol1]):
                for bead, coords in each_mol.items():
                    if bead == b1:
                        for each_bead in coords:
                            c1_xyz.append( list(map(float,each_bead) ) )
            for i, each_mol in enumerate(FRAME_DATA[self.box][mol2]):
                for bead, coords in each_mol.items():
                    if bead == b2:
                        for each_bead in coords:
                            c2_xyz.append( list(map(float,each_bead) ) )
            c1_xyz = np.array(c1_xyz)
            c2_xyz = np.array(c2_xyz)
            if len(c1_xyz) == 0 or len(c2_xyz) == 0:
                self.ignored_frames += 1
                continue
            frame_boxlengths = self.boxlengths[iframe][self.box]
            self.boxLengths.append(frame_boxlengths)
            self.g_frame( c1_xyz, c2_xyz, frame_boxlengths)
        if self.ignored_frames > 0:
            print('Out of %i total frames, %i were ignored b/c of lack of mlcl'%(
                    len(self.frame_data),self.ignored_frames))
        self.averageG()

class G(Struc):
    def __init__(self):
        my_parser = Results()
        my_parser.multi()
        my_parser.parser.add_argument('-B1','--bead1',help='bead from for rdf',
                                      type = str,nargs='+')
        my_parser.parser.add_argument('-M1','--mol1',help='mol from for rdf',
                                      type = str,nargs='+')
        my_parser.parser.add_argument('-B2','--bead2',help='bead to for rdf',
                                      type = str,nargs='+')
        my_parser.parser.add_argument('-M2','--mol2',help='mol to for rdf',
                                      type = str,nargs='+')
        my_parser.parser.add_argument('-bins','--bins',help='radial bin size for rdf',
                                      type = float,default=0.2)
        my_parser.parser.add_argument('-r','--rmax',help='max radius for rdf',
                                      type = float, default=14.0)
        my_parser.parser.add_argument('-d','--dimensions',help='Dimensions for RDF',
                                      type = str,nargs='+',default=['x','y','z'])
        args = vars(my_parser.parse_args())
        self.args = args
        self.checks()

    def checks(self):
        assert (self.args['bead1'] and self.args['bead2'] and self.args['mol1']
                 and self.args['mol2'] and self.args['bins'] and self.args['rmax']), 'Necessary input missing'
        assert self.args['boxes'], 'box needed for rdf'
        self.analysis_class = RDF

    def myCalcs(self, D, m1, b1, m2, b2):
        for box in self.args['boxes']:
            D.calculateRDF(m1, b1, m2, b2, box)
            file_path = '%s/%s/'%(self.args['path'],self.feed)
            file_info = 'box%s_mol%s-%s_mol%s-%s.dat'%(box, m1, b1,m2,b2)
            if len(self.args['dimensions']) == 3:
                writeXY(D.radii, D.g_average, file_path + 'rdf-' + file_info)
                writeXY(D.radii, D.n_average, file_path + 'nint-' + file_info)
            else:
                writeXY(D.radii, D.g_average, file_path +
                        'rdf%s-'%(''.join(self.args['dimensions'])) + file_info)
                writeXY(D.radii, D.n_average, file_path +
                        'nint%s-'%(''.join(self.args['dimensions'])) + file_info)



    def main(self):
        for feed in self.args['feeds']:
            self.feed = feed
            if self.args['verbosity'] > 0: print('-'*12 + 'Dir is %s'%self.feed + '-'*12)
            analysis = self.read_movies(self.args['rmax'], self.args['bins'], self.args['dimensions'])
            for m1, b1, m2, b2 in zip(self.args['mol1'], self.args['bead1'],
                                        self.args['mol2'],self.args['bead2']):
                self.myCalcs(analysis, m1, b1, m2, b2)

import numpy as np
from MCFlow.parser import Results

if __name__ == '__main__':
    M = G()
    M.main()
