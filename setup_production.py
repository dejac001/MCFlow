'''
For production, to satisfy microscopic reversibility, we need to
      - set iratio >= nstep
      - set iratv >= nstep
in addition, it is good practice to set the maximum x,y,z displacements in
file_restart to the average of these 3 values
      - for both translations and rotations
we can calulate the pressure more frequently, so set iratp = 5
'''
def avg_displ(iline, *maxvalue):
    x, y, z = [float(x) for x in iline.split()]
    avgxyz = sum([x, y, z])/3
    if maxvalue and maxvalue[0] < avgxyz:
        print('average value too large: {}; changing to {}'.format(avgxyz, maxvalue[0]))
        avgxyz = maxvalue[0]
    myline = '  %f       %f       %f'%((avgxyz,)*3)
    return myline

def makeProdFiles(path, lastEquilNum, nstep, imv, total_time):
    'Note: chdir to correct directory first'
    input_data = read_fort4(path + 'fort.4')
    nbox = int(input_data['&mc_shared']['nbox'])
    nmolty = int(input_data['&mc_shared']['nmolty'])
    new_input_data = copy.deepcopy(input_data)
    new_input_data['&analysis']['iratp'] = ' 5'
    new_input_data['&mc_shared']['iratio'] = ' %i'%(nstep+100)
    new_input_data['&mc_volume']['iratv'] = ' %i'%(nstep+100)
    new_input_data['&mc_shared']['nstep'] = ' %i'%nstep
    new_input_data['&analysis']['imv'] = ' %i'%imv
    new_input_data['&mc_shared']['time_limit'] = ' %i'%total_time
    for var in ['iprint','iblock']:
        new_input_data['&analysis'][var] = ' %i'%(math.ceil(nstep/10))
    new_input_data['&mc_volume'].pop( 'allow_cutoff_failure' )
    rcut = {}
    for box in new_input_data['SIMULATION_BOX'].keys():
        rcut[box] = float(new_input_data['SIMULATION_BOX'][box]['rcut'].rstrip('d0'))

    restart_data = read_restart('%s/equil-%i/config.equil-%i'%(path,
                                                               lastEquilNum,lastEquilNum),
                                nmolty, nbox)
    new_restart_data = copy.deepcopy(restart_data)
    for box, value in restart_data['max displacement']['translation'].items():
        for mol, line in value.items():
            new_restart_data['max displacement']['translation'][box][mol] = avg_displ(line, 2*rcut[box]) + '\n'
    write_fort4(new_input_data, path + 'fort.4.prod-1')
    write_restart(new_restart_data, path +'config.prod-1')
    shutil.copy(path +'config.prod-1', path +'fort.77')

import math, copy, shutil
from file_formatting.reader import read_restart, read_fort4
from file_formatting.writer import write_restart, write_fort4
from runUtil import what2Analyze

if __name__ == '__main__':
    from parser import Change

    my_parser = Change()
    my_parser.parser.add_argument('-imv','--imovie', help='Movie output frequency',
                                 type=int, default=1000)

    args = vars(my_parser.parse_args())

    nstep = args['nstep']
    imv = args['imovie']
    time = args['time']

    for feed in args['feeds']:
        if args['verbosity'] > 0: print('-'*12 + 'Dir is %s'%feed + '-'*12)
        for seed in args['indep']:
            my_dir = '%s/%s/%i/'%(args['path'],feed,seed)
            (old_begin, nfiles) = what2Analyze(my_dir, args['type'],
                                                   args['guessStart'],args['interval'])
            if args['verbosity'] > 0:
                print('old_begin = {}, nfiles = {}'.format(old_begin, nfiles))
            makeProdFiles(my_dir, old_begin+nfiles-1, nstep, imv, time)
