from MCFlow.structure_analysis import Struc
from MCFlow.file_formatting.reader import Movie
from MCFlow.calc_tools import fold
from MCFlow.file_formatting.reader import convertMovieCoordsXYZ
from MCFlow.file_formatting import writer

class FindRegion(Struc):
    def __init__(self,vectors,lmn):
        self.abc = vectors
        self.lmn = lmn
        assert len(vectors) == len(lmn), '%i != %i'%(len(vectors), len(lmn))
        self.ABC = [i*j for i,j in zip(self.abc, self.lmn)]

    def init_grids(self):
        unit_cells = []
        nL, nM, nN = self.lmn
        for l in range(nL):
            for m in range(nM):
                for n in range(nN):
                    unit_cells.append('%i%i%i'%(l,m,n))
        self.regionNames = unit_cells

    def getRegion(self, xyz):
        (x,y,z) = fold(xyz, self.ABC)
        (a,b,c) = self.abc
        x_unit_cell = math.floor(x/a)
        y_unit_cell = math.floor(y/b)
        z_unit_cell = math.floor(z/c)
        region = '%i%i%i'%(x_unit_cell,y_unit_cell,z_unit_cell)
        assert region in self.regionNames, 'region not found {}'.format(region)
        return region

    def move_to_one_uc(self, xyz, region):
        (x,y,z) = fold(xyz, self.ABC)
        a,b,c = self.abc
        xuc, yuc, zuc = map(int,region)
        return [x-xuc*a, y-yuc*b, z-zuc*c]

    def getChannel(self, xyz):
        return self.getRegion(xyz)

class Channel(FindRegion):

    def __init__(self):
        my_parser = Results()
        my_parser.parser.add_argument('-ID','--name',help='Name of db for molecule number counting',
                               type=str)
        my_parser.parser.add_argument('-B','--bead',help='bead for structure analysis',default=['COM'],
                                     type = str, nargs = '+')
        my_parser.parser.add_argument('-abc','--vectors', help='unit cell vectors. For folding coordinates into a unit cell',
                                      type = float, nargs = '+', default = [20.022,19.899,13.383])
        my_parser.parser.add_argument('-lmn','--replications', help='unit cell vectors. For folding coordinates into a unit cell',
                                      type = int, nargs = '+')
        my_args = vars(my_parser.parse_args())
        if my_args['verbosity'] > 0:
            print(my_args)
        self.args = my_args
        assert len(self.args['replications']) == 3, 'No replications provided: {}'.format(self.args['replications'])
        FindRegion.__init__(self, self.args['vectors'],self.args['replications'])
        self.numIndep = self.args['indep']
        self.checks()
        self.init_grids()
        self.count = {name:[0 for i in self.numIndep]
                            for name in self.regionNames}
        self.data = {region:[] for region in self.regionNames}
        self.averages = {region:{} for region in self.regionNames}

    def checks(self):
        assert self.args['name'], 'Output ID name needed to output local structure info'
        assert 'box' in self.args['box'], 'Need to specify box for structure analysis'
        assert self.args['replications'] != None, 'Supercell construction needs to be specifice'
        self.analysis_class = Movie

    def addToGrid(self, xyz, iseed):
        region = self.getRegion(xyz)
        self.count[region][iseed] += 1
    #
    def __str__(self):
        rep = 'REGION      LOADING (mlcl)\n'
        for region, values in self.averages.items():
            rep += '%12s%12.4f +/-%12.4f\n'%(region, values['mean'],values['stderr'])
        return rep

    def __repr__(self):
        Rep = 'REGION      LOADING (mlcl)\n'
        for region, values in self.averages.items():
            Rep += '   %-9s & %s\n'%(region, tabularNumbers(values['mean'],
                                                      values['stderr']))
        return Rep

    def avgIndep(self):
        self.averages[self.feed] = {}
        for mol in self.N.keys():
            self.averages[self.feed][mol] = {}
            N_total = 0.
            for channel, values in self.N[mol].items():
                mean = np.mean(values)
                self.averages[self.feed][mol][channel] = {}
                self.averages[self.feed][mol][channel]['mean'] = mean
                self.averages[self.feed][mol][channel]['stdev'] = np.std(values)
                self.averages[self.feed][mol][channel]['raw'] = values
                N_total = N_total + mean
            if N_total < 1e-8:
                self.averages[self.feed].pop(mol)
            else:
                self.averages[self.feed][mol]['sum'] = N_total

    def normalize(self, frame_seeds):
        '''
        normalize by total amount of frames
        '''
        indepRange = self.args['indep']
        frame_by_seed = [frame_seeds.count(i) for i in indepRange]
        for mol in self.N.keys():
            for region in self.N[mol].keys():
                for i, val in enumerate(frame_by_seed):
                    self.N[mol][region][i] = self.N[mol][region][i]/val

    def countBins(self, frame_data, frame_seeds):
        box, indepRange = self.args['box'], self.args['indep']
        total_frames = -1
        for iframe, FRAME_DATA in enumerate(frame_data):
            if total_frames == -1:
                N = {mol:{region:[0 for i in indepRange]
                                    for region in self.regionNames}
                                        for mol in FRAME_DATA[box].keys()}
                N_frame =  {mol:{region:[   [0 for j in range(frame_seeds.count(i))] for i in indepRange]
                                    for region in self.regionNames}
                                        for mol in FRAME_DATA[box].keys()}

            total_frames += 1
            try:
                if (iframe+1)%(len(frame_data)//4) == 0:
                    print('%5.1f %% of frames analyzed'%(100*(iframe+1)/len(frame_data)))
            except ZeroDivisionError:
                print('%5.1f %% of frames analyzed'%(100*(iframe+1)/len(frame_data)))
            seed = frame_seeds[total_frames]
            seed_index = self.numIndep.index(seed)
            for mol in FRAME_DATA[box].keys():
                molecule_frame = {region:[] for region in self.regionNames}
                for imol, each_molecule in enumerate(FRAME_DATA[box][mol]):
                    beads = list(set(each_molecule.keys()) &
                                set(self.args['bead']))
                    print('beads are', beads)
                    assert len(beads) == 1, 'Ambiguous beads'
                    for bead in beads:
                        for each_coord in each_molecule[bead]:
                            region = self.getRegion(each_coord)
                            try:
                                N[mol][region][seed_index] += 1
                                N_frame[mol][region][seed_index][iframe] += 1
                            except  KeyError:
                                print(region)
                                print(each_coord)
                                raise KeyError
        self.N = N
        self.normalize(frame_seeds)
        self.avgIndep()
        for mol in self.averages[self.feed].keys():
            assert mol in N_frame.keys(), 'Inconsistent dict: {} not in {}'.format(mol, self.averages[self.feed].keys())
            for region, vals in N_frame[mol].items():
                assert region in self.averages[self.feed][mol].keys(), 'Inconsistent dict'
                self.averages[self.feed][mol][region]['raw data'] = vals
        self.averages[self.feed]['total frames'] = total_frames

    def myCalcs(self, D ):
        self.countBins(D.frame_data, D.frame_seed)
        outputDB(self.args['path'],[self.feed],self.args['type'],{self.args['name']: self } )

    def main(self):
        for feed in self.args['feeds']:
            self.feed = feed
            if self.args['verbosity'] > 0: print('-'*12 + 'Dir is %s'%self.feed + '-'*12)
            analysis = self.read_movies()
            self.myCalcs(analysis)

from PythonCode.PaperWriting.makeHTable import tabularNumbers
from MCFlow.parser import Results
from MCFlow.getData import outputDB
import numpy as np
from MCFlow.UCDens import rho_v_r
import math

if __name__ == '__main__':
    M = Channel()
    M.main()
