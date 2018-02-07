def addAtoms(input_dat, restart_dat, coords, mol_num, charge):
    for xyz in coords:
        restart_dat['mol types'].append( mol_num )
        restart_dat['box types'].append( '1' ) # always box 1
        restart_dat['coords'].append( [
                {'xyz': '%f %f %f\n'%(xyz[0],xyz[1],xyz[2]),
                 'q': '%f\n'%charge}])
        restart_dat['nchain'] = restart_dat['nchain'].replace(
            restart_dat['nchain'].split()[0], '%i'%(int(restart_dat['nchain'].split()[0])+1)
        )
        input_dat['SIMULATION_BOX']['box1']['mol%s'%mol_num] = '%i'%(
            int(input_dat['SIMULATION_BOX']['box1']['mol%s'%mol_num]) + 1
        )
    input_dat['&mc_shared']['nchain'] = restart_dat['nchain'].split()[0]
    return copy.deepcopy(input_dat), copy.deepcopy(restart_dat)

from MCFlow.file_formatting import reader, writer
import copy
import numpy as np

if __name__ == '__main__':
    import argparse, os
    parser = argparse.ArgumentParser(description='add explicit to unit cell')
    parser.add_argument('-f','--file',type=str)
    parser.add_argument('-lmn','--vectors',type=int,nargs='+',default=[2,2,3])
    parser.add_argument('-i','--input',default='fort.4')
    parser.add_argument('-r','--restart',default='fort.77')

    args = vars(parser.parse_args())

    base_dir = os.getcwd() + '/'
    data = reader.PDB(args['file'])
    input_data = reader.read_fort4('%s%s'%(base_dir,args['input']))
    nmolty, nbox = (int(input_data['&mc_shared']['nmolty']),
                int(input_data['&mc_shared']['nbox']))
    restart_data = reader.read_restart('%s%s'%(base_dir,args['restart']),nmolty, nbox)

    si_atoms = [data['coords'][i] for i in range(len(data['coords'])) if data['atoms'][i] == 'Si']
    o_atoms = [data['coords'][i] for i in range(len(data['coords'])) if data['atoms'][i] == 'O']
    assert len(si_atoms) + len(o_atoms) == len(data['coords']), 'Miss counted'
    input_data, restart_data = addAtoms(input_data, restart_data, si_atoms, '2', 1.5)
    input_data, restart_data = addAtoms(input_data, restart_data, o_atoms, '3', -0.75)
    a,b,c = data['box info']['a'], data['box info']['b'], data['box info']['c']
    alpha, gamma, beta = data['box info']['alpha'], data['box info']['gamma'], data['box info']['beta']
    if alpha == 90.0 and gamma == 90.0 and beta == 90.0:
        restart_data['box dimensions']['box1'] = '%f %f %f\n'%(a,b,c)
    else:
        alp = alpha/180*np.pi
        bet = beta/180*np.pi
        gam = gamma/180*np.pi
        restart_data['box dimensions']['box1'] = (
        '  %f %f %f  \n'%(a*np.sin(bet), b*np.sin(alp)*np.cos(gam), 0.) +
        '  %f %f %f  \n'%(0., b*np.sin(alp)*np.sin(gam), 0.) +
        '  %f %f %f  \n'%(a*np.cos(bet), b*np.cos(alp), c)
            )
        old_def = input_data['SIMULATION_BOX']['box1']['defaults'].split() # '3.5 0.000 0 F F F F'
        old_def[-4] = 'T'
        input_data['SIMULATION_BOX']['box1']['defaults'] = ' '.join(k for k in old_def)

    writer.write_fort4(input_data,'fort.4.new')
    writer.write_restart(restart_data,'fort.77.new')
