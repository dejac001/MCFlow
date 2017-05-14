def findHydroxylHydrogen(Oxyz, Hcoords, abc):
    for Hxyz in Hcoords:
        if calculate_distance2(Oxyz, Hxyz, abc) < 1.0:
            return Hxyz

def findHB(beadsFrom, molsTo, abc, criteria):
    aOHO_min = 150./180*3.1415926535897931 # degrees
    rOH_max = 6.25 # 2.5*2.5
    rOO_max = 10.89 # 3.3*3.3
    nHB = 0
    # # iterate through all Oxygens
    # for O1 in beadsFrom['O']:
    #     # iterate through all H on molecule that are bonded to O1
    #     for H1 in [i for i in beadsFrom['H'] if calculate_distance2(O1,i,abc) < 1.0]:
    #         # iterate through all molecules of those specified in box
    #         for mol in molsTo:
    #             # iterate through Oxygens of selected molecule if selected molecule is not the same as O1
    #             for O2 in [j for j in mol['O'] if calculate_distance2(O1,j,abc) > 0.1]:
    #                 # do not consider same molecule
    #                 for H2 in [i for i in mol['H'] if calculate_distance2(O2,i,abc) < 1.0]:
    # iterate through all Oxygens
    for O1 in beadsFrom['O']:
        # iterate through all molecules of those specified in box
        for mol in molsTo:
            # iterate through Oxygens of selected molecule if selected molecule is not the same as O1
            for O2 in [j for j in mol['O'] if calculate_distance2(O1,j,abc) > 0.1]:
                # iterate through all H on molecule 1 that are bonded to O1
                our_HB = 0
                for H1 in [i for i in beadsFrom['H'] if calculate_distance2(O1,i,abc) < 1.0]:
                    # iterate through all H on molecule 2 that are bonded to O2
                    for H2 in [i for i in mol['H'] if calculate_distance2(O2,i,abc) < 1.0]:
                        if criteria == 'loose':
                            if ((calculate_distance2(O1, O2, abc) < rOO_max) and
                                ((calculate_distance2(O1, H2, abc) < rOH_max) or
                                (calculate_distance2(O2, H1, abc) < rOH_max))):
                                nHB += 1
                                our_HB += 1
                                break # can only have 1 hydrogen bond between two molecules
                        elif criteria == 'strict':
                            if calculate_distance2(H1,O2,abc) < rOH_max:
                                # 1 can be charge acceptor
                                if calculate_angle(O1,H1,O2,abc) > aOHO_min:
                                    our_HB += 1
                                    nHB +=1
                                    break # can only have 1 hydrogen bond between two molecules
                            elif calculate_distance2(H2,O1,abc) < rOH_max:
                                # 2 can be charge acceptor
                                if calculate_angle(O2,H2,O1,abc) > aOHO_min:
                                    nHB += 1
                                    our_HB += 1
                                    break # can only have 1 hydrogen bond between two molecules
                    if our_HB == 1:
                        break # can only have 1 hydrogen bond between two molecules
    return nHB

from MCFlow.file_formatting.reader import Movie

class HydrogenBond(Movie):
    def __init__(self, file_name):
        Movie.__init__(self, file_name)


    def getBeads(self, box):
        '''
        '''
        my_box = 'box%s'%box
        self.HB_info = []
        OTypes = {'62':'alkanol','114':'water','178':'zeo','181':'silanol'}
        HTypes ={'61':'alkanol','115':'water','182':'silanol'}
        for iframe, FRAME_DATA in enumerate(self.frame_data):
            # store H and O for all other mols
            HB_mols = {}
            for molType in FRAME_DATA[my_box].keys():
                if molType not in HB_mols.keys():
                    HB_mols[molType] = []
                for imol, each_molecule in enumerate(FRAME_DATA[my_box][molType]):
                    beads = {'H':[],'O':[]}
                    for bead in each_molecule.keys():
                        if bead in OTypes.keys():
                            for coord in each_molecule[bead]:
                                beads['O'].append( list(map(float,coord)) )
                        elif bead in HTypes.keys():
                            for coord in each_molecule[bead]:
                                beads['H'].append( list(map(float,coord)) )
                    HB_mols[molType].append( beads )
            self.HB_info.append( HB_mols )

    def calcHB(self, mols, box):
        self.getBeads(self, box)
        my_box = 'box%s'%box
        self.HB = []
        for iframe, HB_data in enumerate(self.HB_info):
            try:
                if (iframe+1)%(self.nframes//4) == 0:
                    print('%5.1f %% of frames analyzed for HB'%(100*(iframe+1)/self.nframes))
            except ZeroDivisionError:
                print('%5.1f %% of frames analyzed for HB'%(100*(iframe+1)/self.nframes))
            self.HB.append( {my_box:{}  } )
            HB_to_mols = []
            for mol1 in mols:
                HB_to_mols.append( HB_data[mol1] )
            for mol2 in HB_data.keys():
                pair = mol1 + '-' + mol2
                self.HB[iframe][my_box][pair] = []
                for HB_from_beads in HB_data[mol2]:
                    # determine hydrogen bonds with molecules of interest
                    my_nHB = findHB(HB_from_beads, HB_to_mols, self.boxlengths[iframe][my_box], self.criteria)
                    self.HB[iframe][my_box][pair].append( my_nHB )
                    
    def filterHB(self):
        for iframe, FRAME_DATA in enumerate(self.frame_data):
            for box, values in self.HB[iframe].items():
                for mol in values.keys():
                    indices_to_keep = []
                    for imol, nHB in enumerate(values[mol]):
                        if nHB > 0:
                            indices_to_keep.append( imol )
                    FRAME_DATA[box][mol] = [value for i, value in enumerate(FRAME_DATA[box][mol])
                                                   if i in indices_to_keep]
    def countHB(self, nIndep, feed, data):
        new = self.__class__('HB-data')
        for attr, value in self.__dict__.items():
            new.__dict__[attr] = value
        new.countMols(nIndep, feed, data)
        return new

class HB_format:
    def __init__(self):
        my_parser = MultMols()
        #TODO: make able to do multiple boxes at same time
        my_parser.parser.add_argument('-ID','--name',help='Name of db for molecule number counting',
                               type=str,default = '')
        my_parser.parser.add_argument('-H','--htype',help='hydrogen bonding criteria type',
                                      type=str, choices = ['loose','strict'],default='strict')
        my_args = vars(my_parser.parse_args())
        self.args = my_args
        self.checks()

    def checks(self):
        assert self.args['name'], 'Output ID name needed to output local structure info'
        self.HB_class = HydrogenBond

    def read_movies(self):
        for seed in self.args['indep']:
            my_dir = '%s/%s/%i'%(self.args['path'],self.feed,seed)
            (old_begin, nfiles) = what2Analyze(my_dir, self.args['type'],
                                                       self.args['guessStart'],
                                               self.args['interval'])
            if self.args['verbosity'] > 0:
                print('old_begin = {}, nfiles = {}'.format(old_begin, nfiles))
            for fileNum in range(old_begin, old_begin+nfiles):
                movie_file = '%s/%s%i/movie.%s%i'%(my_dir, self.args['type'],fileNum,
                                                   self.args['type'],fileNum)
                if (fileNum == old_begin) and (seed == self.args['indep'][0]):
                    # keep track of info only for each feed
                    I = self.HB_class(movie_file)
                    I.read_header()
                    I.read_movie_frames(seed)
                else:
                    F = self.HB_class(movie_file)
                    F.read_header()
                    F.read_movie_frames(seed)
                    I = I + F
        return I

    def myCalcs(self, D ):
        D.criteria = self.args['htype']
        D.calcHB(self.args['mol'], self.args['box'])
        HB = D.countHB(self.args['indep'],feed, D.HB)
        D.filterHB()
        D.countMols(self.args['indep'],feed, D.frame_data)
        for mol_num in D.averages[feed].keys():
            xyz_data = D.getCoords(mol_num, self.args['box'], ['COM'])
            xyz('%s/%s/movie_HB_coords_mol%s_box%s.xyz'%(self.args['path'], feed,
                                                            mol_num, self.args['box']), xyz_data)
            #TODO: make so updates recent file--does not overwrite it
        outputDB(self.args['path'],[feed],self.args['type'],{self.args['name']+'-' + self.args['htype'] + '-box'+self.args['box']:HB})
    #                self.args['name']: D })
                    #{self.args['name']: D
    #                       'HB':HB } )

    def main(self):
        for feed in self.args['feeds']:
            self.feed = feed
            if self.args['verbosity'] > 0: print('-'*12 + 'Dir is %s'%self.feed + '-'*12)
            HB_class = self.read_movies()
            self.myCalcs(HB_class)
    

from MCFlow.runAnalyzer import what2Analyze
from MCFlow.file_formatting.writer import xyz
from MCFlow.calc_tools import calculate_distance2, calculate_angle
from MCFlow.parser import MultMols
from MCFlow.getData import outputDB

if __name__ == '__main__':
    M = HB_format()
    M.main()
