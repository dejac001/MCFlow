def findHydroxylHydrogen(Oxyz, Hcoords, abc):
    for Hxyz in Hcoords:
        if calculate_distance2(Oxyz, Hxyz, abc) < 1.0:
            return Hxyz

def findHB(beadsFrom, molsTo, abc, criteria):
    aOHO_min = 150./180*3.1415926535897931 # degrees
    rOH_max = 6.25 # 2.5*2.5
    rOO_max = 10.89 # 3.3*3.3
    nHB = 0
    for O1 in beadsFrom['O']:
        for H1 in [i for i in beadsFrom['H'] if calculate_distance2(O1,i,abc) < 1.0]:
            for mol in molsTo:
                for O2 in [j for j in mol['O'] if calculate_distance2(O1,j,abc) > 0.1]:
                    # do not consider same molecule
                    for H2 in [i for i in mol['H'] if calculate_distance2(O2,i,abc) < 1.0]:
                        if criteria == 'loose':
                            if ((calculate_distance2(O1, O2, abc) < rOO_max) and
                                ((calculate_distance2(O1, H2, abc) < rOH_max) or
                                (calculate_distance2(O2, H1, abc) < rOH_max))):
                                nHB += 1
                        elif criteria == 'strict':
                            if calculate_distance2(H1,O2,abc) < rOH_max:
                                # 1 can be charge acceptor
                                if calculate_angle(O1,H1,O2,abc) > aOHO_min:
                                    nHB +=1
                            elif calculate_distance2(H2,O1,abc) < rOH_max:
                                # 2 can be charge acceptor
                                if calculate_angle(O2,H2,O1,abc) > aOHO_min:
                                    nHB += 1
    return nHB

from MCFlow.file_formatting.reader import Movie

class HydrogenBond(Movie):
    def __init__(self, file_name):
        Movie.__init__(self, file_name)

    def calcHB(self, mols, box):
        '''
        In this case, hydrogen bonding will be found for anything hydrogen bonding with mol
        Only keep molecules in movie files that are hydrogen bonding to mol
        '''
        my_box = 'box%s'%box
        self.HB = []
        OTypes = {'62':'alkanol','114':'water','178':'zeo','181':'silanol'}
        HTypes ={'61':'alkanol','115':'water','182':'silanol'}
        for iframe, FRAME_DATA in enumerate(self.frame_data):
            try:
                if (iframe+1)%(self.nframes//4) == 0:
                    print('%5.1f %% of frames analyzed for HB'%(100*(iframe+1)/self.nframes))
            except ZeroDivisionError:
                print('%5.1f %% of frames analyzed for HB'%(100*(iframe+1)/self.nframes))
            self.HB.append( {my_box:{}  } )
            # store H and O for mols that you want to calculate hydrogen bonds with
            # i.e., to_mols = mols you specified
            HB_to_mols = []
            for mlcl in mols:
                HB_to_beads = {'H':[],'O':[]}
                for each_molecule in FRAME_DATA[my_box]['mol%s'%mlcl]:
                    for bead in each_molecule.keys():
                        if bead in OTypes.keys():
                            for coord in each_molecule[bead]:
                                HB_to_beads['O'].append( list(map(float,coord)) )
                        elif bead in HTypes.keys():
                            for coord in each_molecule[bead]:
                                HB_to_beads['H'].append( list(map(float,coord)) )
                HB_to_mols.append( HB_to_beads )
            # store H and O for all other mols
            for molType in FRAME_DATA[my_box].keys():
                self.HB[iframe][my_box][molType] = []
                for imol, each_molecule in enumerate(FRAME_DATA[my_box][molType]):
                    HB_from_beads = {'H':[],'O':[]}
                    for bead in each_molecule.keys():
                        if bead in OTypes.keys():
                            for coord in each_molecule[bead]:
                                HB_from_beads['O'].append( list(map(float,coord)) )
                        elif bead in HTypes.keys():
                            for coord in each_molecule[bead]:
                                HB_from_beads['H'].append( list(map(float,coord)) )
                    # determine hydrogen bonds with molecules of interest
                    my_nHB = findHB(HB_from_beads, HB_to_mols, self.boxlengths[iframe][my_box], self.criteria)
                    self.HB[iframe][my_box][molType].append( my_nHB )
                    
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

def main():
    my_parser = MultMols()
    #TODO: make able to do multiple boxes at same time
    my_parser.parser.add_argument('-ID','--name',help='Name of db for molecule number counting',
                           type=str,default = '')
    my_parser.parser.add_argument('-H','--htype',help='hydrogen bonding criteria type',
                                  type=str, choices = ['loose','strict'],default='strict')
    args = vars(my_parser.parse_args())

    assert args['name'], 'Output ID name needed to output local structure info'

    for feed in args['feeds']:
        if args['verbosity'] > 0: print('-'*12 + 'Dir is %s'%feed + '-'*12)
        for seed in args['indep']:
            my_dir = '%s/%s/%i'%(args['path'],feed,seed)
            (old_begin, nfiles) = what2Analyze(my_dir, args['type'],
                                                       args['guessStart'],
                                               args['interval'])
            if args['verbosity'] > 0:
                print('old_begin = {}, nfiles = {}'.format(old_begin, nfiles))
            for fileNum in range(old_begin, old_begin+nfiles):
                movie_file = '%s/%s%i/movie.%s%i'%(my_dir, args['type'],fileNum,
                                                   args['type'],fileNum)
                if (fileNum == old_begin) and (seed == args['indep'][0]):
                    # keep track of info only for each feed
                    D = HydrogenBond(movie_file)
                    D.read_header()
                    D.read_movie_frames(seed)
                else:
                    F = HydrogenBond(movie_file)
                    F.read_header()
                    F.read_movie_frames(seed)
                    D = D + F
        D.criteria = args['htype']
        D.calcHB(args['mol'], args['box'])
        HB = D.countHB(args['indep'],feed, D.HB)
        D.filterHB()
        D.countMols(args['indep'],feed, D.frame_data)
        for mol_num in D.averages[feed].keys():
            xyz_data = D.getCoords(mol_num, args['box'], ['COM'])
            xyz('%s/%s/movie_HB_coords_mol%s_box%s.xyz'%(args['path'], feed,
                                                        mol_num, args['box']), xyz_data)
        #TODO: make so updates recent file--does not overwrite it
        outputDB(args['path'],[feed],args['type'],{args['name']+'-' + args['htype'] + '-box'+args['box']:HB})
#                args['name']: D })
                #{args['name']: D
#                       'HB':HB } )
    

from MCFlow.runAnalyzer import what2Analyze
from MCFlow.file_formatting.writer import xyz
from MCFlow.calc_tools import calculate_distance2, calculate_angle
from MCFlow.parser import MultMols
from MCFlow.getData import outputDB

if __name__ == '__main__':
    main()
