from MCFlow.file_formatting.reader import Movie
from MCFlow.structure_analysis import Struc

def writeXY(x,y,name):
    with open(name,'w') as f:
        for i,j in zip(x,y):
            f.write('%e %e\n'%(i,j))

class S_k(Movie):

    def __init__(self, file_name, *args):
        Movie.__init__(self, file_name, *args)
        self.tol = 1e-10 # tolerance for |q|
        self.maxN = 10 # max for nx, ny, nz
        self.maxQ = 1.
        self.s_by_frame = {'S(q)':[],'|q|':[]}

    def averageS(self):
        self.s_average = {}
        num_q_mag = set(len(i) for i in self.s_by_frame['|q|'])
        assert len(num_q_mag) == 1, 'All frames did not have same |q|; set = {}'.format(num_q_mag)
        self.s_average['|q|'] = self.s_by_frame['|q|'][0]
        self.s_average['S(q)'] = np.mean(self.s_by_frame['S(q)'],axis=0)
        self.s_by_frame = {'S(q)':[],'|q|':[]} # initialize for next box

    def get_hkl(self, boxLengths, q):
        my_hkl = []
        Lx, Ly, Lz = boxLengths
        for nx in range(-self.maxN, self.maxN + 1):
            for ny in range(-self.maxN, self.maxN + 1):
                for nz in range(-self.maxN, self.maxN + 1):
                    q_vect = np.array([2*np.pi*nx/Lx, 2*np.pi*ny/Ly, 2*np.pi*nz/Lz])
                    q_mag = np.linalg.norm(q_vect, 2, 0)
                    if (((nx == 0) and (ny == 0) and (nz == 0))
                        or (q_mag > self.maxQ)):
                        continue
                    if abs(q_mag - q) < self.tol:
                        #magnitude is the same
                        my_hkl.append([nx, ny, nz])
        assert len(my_hkl) > 0, 'No hkl found'
        return my_hkl


    def s_frame(self, coords, boxLengths ):
        '''
        coords: [[x1,y1,z1],,,[xN,yN,zN]]
        '''
        Lx, Ly, Lz = boxLengths
        self.boxLengths = Lx, Ly, Lz
        n_site = len(coords)
        S_func = {'|q|':[],'S(q)':[]}
        for nx in range(-self.maxN, self.maxN + 1):
            for ny in range(-self.maxN, self.maxN + 1):
                for nz in range(-self.maxN, self.maxN + 1):
                    q_vect = np.array([2*np.pi*nx/Lx, 2*np.pi*ny/Ly, 2*np.pi*nz/Lz])
                    q_mag = np.linalg.norm(q_vect, 2, 0)
                    if (((nx == 0) and (ny == 0) and (nz == 0))
                        or (q_mag > self.maxQ)):
                        continue
                    q_dot_r = np.dot(q_vect, coords.T)
                    my_cos = np.sum(np.cos(q_dot_r))
                    my_sin = np.sum(np.sin(q_dot_r))
                    S_q = (my_cos*my_cos + my_sin*my_sin) / n_site
                    mag_found = False
                    for i, q in enumerate(S_func['|q|']):
                        if abs(q_mag - q) < self.tol:
                            #magnitude is the same
                            mag_found = True
                            S_func['S(q)'][i].append( S_q )
                    if mag_found == False:
                        S_func['|q|'].append( q_mag )
                        S_func['S(q)'].append( [S_q] )
        for i in range(len(S_func['|q|'])): #TODO: speedup here with map
            S_func['S(q)'][i] = np.mean(S_func['S(q)'][i])
        for key, value in S_func.items():
            self.s_by_frame[key].append( value )

    def calculate(self, sbox, beads):
        self.box = 'box%s'%sbox
        for iframe, FRAME_DATA in enumerate(self.frame_data):
            try:
                if (iframe+1)%(self.nframes//4) == 0:
                    print('%5.1f %% of frames analyzed'%(100*(iframe+1)/self.nframes))
            except ZeroDivisionError:
                print('%5.1f %% of RDF frames completed'%(100*(iframe+1)/self.nframes))
            xyz = []
            for mol, molInfo in FRAME_DATA[self.box].items():
                for chain in molInfo:
                    for bead, coords in chain.items():
                        if bead not in beads: continue
                        for each_bead in coords:
                            xyz.append( list(map(float,each_bead)) )
            self.s_frame( np.matrix(xyz), self.boxlengths[iframe][self.box])
        self.averageS()

class S(Struc):
    def __init__(self):
        my_parser = Results()
        my_parser.multi()
        my_parser.parser.add_argument('-E','--beads',help='beads for structure analysis',default=['62'],
                                     type = str, nargs = '+')
        my_parser.parser.add_argument('-abc','--vectors', help='unit cell vectors. For folding coordinates into a unit cell',
                                      type = float, nargs = '+', default = [20.022,19.899,13.383])
        args = vars(my_parser.parse_args())
        self.args = args
        self.checks()

    def checks(self):
        assert self.args['boxes'], 'box needed for structure factor'
#       assert self.args['box'] == None, 'Use multiple boxes argument'
        self.analysis_class = S_k

    def myCalcs(self, D):
        for box in self.args['boxes']:
            if box == '1':
                print('folding...')
                D.foldMovieToUC(self.args['vectors'])
            D.calculate(box, self.args['beads'])
            x, y = D.s_average['|q|'], D.s_average['S(q)']
            writeXY(x, y, '%s/%s/Sfact-box%s-%s.dat'%(self.args['path'],self.feed,box,''.join(self.args['beads'])))


    def main(self):
        for feed in self.args['feeds']:
            self.feed = feed
            if self.args['verbosity'] > 0: print('-'*12 + 'Dir is %s'%self.feed + '-'*12)
            analysis = self.read_movies()
            print('Done reading movie files')
            self.myCalcs(analysis)

import numpy as np
from MCFlow.parser import Results

if __name__ == '__main__':
    M = S()
    M.main()
