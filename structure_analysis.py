class Region:
    def __init__(self, func):
        '''
        func: function that returns for a given xyz coordinate whether or not
                a bead is in a defined region
        '''
        self.my_region = func

    def COM(self, molecule_data):
        '''
        :param molecule_data: {'atoms': [list,of,bead,nums],'coords':[xyz1,xyz2...]}
        determine if COM within region
        '''
        xyz = molecule_data['coords'][molecule_data['atoms'].index('COM')]
        return self.my_region(xyz)

    def allBeads(self, molecule_data):
        '''
        determine if all beads are in region
        '''
        number_in_region = 0
        for xyz in molecule_data['coords']:
            if self.my_region(xyz):
                number_in_region += 1
        if number_in_region == len(molecule_data['coords']):
            return True
        else:
            return False

def intersection(coords):
    rcut = 5.0/2.
    rcut_2 = rcut*rcut
    a, b, c = [20.022, 19.899, 13.383] # MFI
    center = [10.011,4.9748,0.] # (a/2, b/4, 0)
    x, y, z = center
    sphere = [
            [x,y,z],
            [a/2-x,b/2+y,c/2+z],
            [x,b/2-y,z],
            [a/2+x,y,c/2-z],
            [a-x,b-y,c-z],
            [a/2+x,b/2-y,c/2-z],
            [a-x,b/2+y,c-z],
            [a/2-x,b-y,c/2+z]
        ]
    for my_center in sphere:
        if calculate_distance2(my_center, coords, [a,b,c]) < rcut_2:
            return True
    return False

def SnSPP(coords):
    x,y,z = coords
    if ((z > 2.73) and (z < 32.12)):
        return True
    else:
        return False

class Struc:
    def __init__(self, filterFunc=None):
        my_parser = Results()
        my_parser.parser.add_argument('-ID','--name',help='Name of db for molecule number counting',
                               type=str)
        my_parser.parser.add_argument('-B','--bead',help='bead for structure analysis',default=['COM'],
                                     type = str, nargs = '+')
        my_parser.parser.add_argument('-abc','--vectors', help='unit cell vectors. For folding coordinates into a unit cell',
                                      type = float, nargs = '+', default = [20.022,19.899,13.383])
        my_args = vars(my_parser.parse_args())
        self.args = my_args
        self.filter_function = filterFunc
        self.checks()

    def checks(self):
        if self.args['vectors']:
            assert self.args['box'] == '1', 'It doesnt make sense to fold into a non replicated cell'
        if self.filter_function:
            assert self.args['name'], 'Output ID name needed to output local structure info'
        assert self.args['box'], 'Need to specify box for structure analysis'
        self.analysis_class = Movie

    def read_movies(self, *args):
        for seed in self.args['indep']:
            my_dir = '%s/%s/%i'%(self.args['path'],self.feed,seed)
            (old_begin, nfiles) = what2Analyze(my_dir, self.args['type'],
                                                       self.args['guessStart'],
                                               self.args['interval'])
            if self.args['verbosity'] > 0:
                print('old_begin = {}, nfiles = {}'.format(old_begin, nfiles))
            for fileNum in range(old_begin, old_begin+nfiles):
                movie_file = fo.read(my_dir,'movie.',self.args['type'], fileNum)
                if (fileNum == old_begin) and (seed == self.args['indep'][0]):
                    # keep track of info only for each feed
                    I = self.analysis_class(movie_file, *args)
                    I.read_header(self.args['ignore'])
                    I.read_movie_frames(seed)
                else:
                    F = self.analysis_class(movie_file, *args)
                    F.read_header(self.args['ignore'])
                    F.read_movie_frames(seed)
                    I = I + F
        return I

    def myCalcs(self, D ):
        if self.filter_function:
            D.filterCoords(self.filter_function, self.args['box'])
        D.countMols(self.args['indep'],self.feed, D.frame_data)
        if self.args['vectors']: D.foldMovieToUC(self.args['vectors'])
        for mol_num in D.averages[self.feed].keys():
            if mol_num == 'number of frames': continue
            xyz_data = D.getCoords(mol_num, self.args['box'], self.args['bead'])
            flag=''
            if self.filter_function:
                flag = '_filtered'
            xyz('%s/%s/movie_coords_mol%s_box%s%s_sim%s_%s.xyz'%(self.args['path'], self.feed,
                                                        mol_num, self.args['box'],flag, ''.join(map(str,self.args['indep'])),
                                                        '-'.join(self.args['bead'])), xyz_data)
        if self.args['name']:
            outputDB(self.args['path'],[self.feed],self.args['type'],{self.args['name']: D } )

    def main(self):
        for feed in self.args['feeds']:
            self.feed = feed
            if self.args['verbosity'] > 0: print('-'*12 + 'Dir is %s'%self.feed + '-'*12)
            analysis = self.read_movies()
            self.myCalcs(analysis)


from MCFlow.runAnalyzer import what2Analyze
from MCFlow.file_formatting.reader import Movie
from MCFlow.file_formatting.writer import xyz
from MCFlow.parser import Results
from MCFlow.getData import outputDB
from MCFlow.calc_tools import calculate_distance2
import MCFlow.file_organization as fo


if __name__ == '__main__':
    M = Struc()
    M.main()
