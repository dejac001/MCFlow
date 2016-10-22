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


def main(filter_function=None):
    my_parser = Results()
    my_parser.parser.add_argument('-ID','--name',help='Name of db for molecule number counting',
                           type=str,default = '')
    args = vars(my_parser.parse_args())

    N_database = {}
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
                    D = Movie(movie_file)
                    D.read_header()
                    D.read_movie_frames()
                else:
                    F = Movie(movie_file)
                    F.read_header()
                    F.read_movie_frames()
                    D = D + F
        if filter_function:
            D.filterCoords(filter_function, args['box'])
        D.countMols(len(args['indep']),feed)
        xyz_data = D.getCoords(args['mol'], args['box'], ['COM'])
        xyz('%s/%s/movie_coords_mol%s_box%s.xyz'%(args['path'], feed, args['mol'], args['box']), xyz_data)
    if args['name']:
        outputDB(args['path'],args['feeds'],args['type'],{args['name']: D } )

from MCFlow.runAnalyzer import what2Analyze
from MCFlow.file_formatting.reader import Movie
from MCFlow.file_formatting.writer import xyz
from MCFlow.parser import Results
from MCFlow.getData import outputDB

if __name__ == '__main__':
    main()
