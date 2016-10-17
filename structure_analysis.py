

from runAnalyzer import what2Analyze
from file_formatting.reader import Movie
from file_formatting.writer import xyz

if __name__ == '__main__':
    from parser import Results
    my_parser = Results()
    args = vars(my_parser.parse_args())

    for feed in args['feeds']:
        for seed in args['indep']:
            my_dir = '%s/%s/%i'%(args['path'],feed,seed)
            (old_begin, nfiles) = what2Analyze(my_dir, args['type'],
                                                       args['guessStart'],
                                               args['interval'])
            for fileNum in range(old_begin, old_begin+nfiles):
                movie_file = '%s/%s%i/movie.%s%i'%(my_dir, args['type'],fileNum,
                                                   args['type'],fileNum)
                if (fileNum == old_begin) and (seed == args['indep'][0]):
                    M = Movie(movie_file)
                    M.read_header()
                    M.read_movie_frames()
                else:
                    N = Movie(movie_file)
                    N.read_header()
                    N.read_movie_frames()
                    M = M + N
        xyz_data = M.getCoords(args['mol'], args['box'], ['COM'])
        xyz('%s/%s/movie_coords_mol%s.xyz'%(args['path'], feed, mol), xyz_data)