def filmUC(molecule_data):
    '''
    :param molecule_data: {'atoms': [list,of,bead,nums],'coords':[xyz1,xyz2...]}
    :return: whether or not molecule_data is in region
    '''
    number_in_region = 0
    '''
    all atoms in position
    for xyz in molecule_data['coords']:
        if (xyz[0] > 29.953) and (xyz[0] < 49.975):
            number_in_region += 1
    if number_in_region == len(molecule_data['coords']):
        return True
    else:
        return False
    '''
    xyz = molecule_data['coords'][molecule_data['atoms'].index('COM')]
    if (xyz[0] > 29.953) and (xyz[0] < 49.975):
        return True
    else:
        return False

from runAnalyzer import what2Analyze
from file_formatting.reader import Movie
from file_formatting.writer import xyz

if __name__ == '__main__':
    from parser import Results
    from getData import outputDB
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
        D.filterCoords(filmUC, args['box'])
        D.countMols(len(args['indep']),feed)
        xyz_data = D.getCoords(args['mol'], args['box'], ['COM'])
        xyz('%s/%s/movie_coords_mol%s_box%s.xyz'%(args['path'], feed, args['mol'], args['box']), xyz_data)
    if args['name']:
        outputDB(args['path'],args['feeds'],args['type'],{args['name']: D } )

