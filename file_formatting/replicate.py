def get_H_matrix(a,b,c,alpha,beta,gamma):
    return np.matrix([[a*np.sin(beta), b*np.sin(alpha)*np.cos(gamma),0.],
                    [0., b*np.sin(alpha)*np.sin(gamma), 0.],
                    [a*np.cos(beta), b*np.cos(alpha), c]])

def replicate(coordinates, abc, lmn, alpha=90.,gamma=90.,beta=90.):
    '''

    :param coordinates: data type read in by xyz or pdb
    :param abc: box dimensions
    :param lmn: integer number of cells in each dimension
    :return:
    '''
    a, b, c = abc
    replicated = {'atoms':[], 'coords':[]}
    replicated['box info'] = {'a':a*lmn[0],'b':b*lmn[1],'c':c*lmn[2],
                                'alpha':alpha,'gamma':gamma,'beta':beta}
    alpha = alpha/180.*np.pi
    beta = beta/180.*np.pi
    gamma = gamma/180.*np.pi
    H = get_H_matrix(a,b,c,alpha,beta,gamma)
    for xcell in range(lmn[0]):
        for ycell in range(lmn[1]):
            for zcell in range(lmn[2]):
                for atom, xyz in zip(coordinates['atoms'], coordinates['coords']):
                    xyz_cryst = np.linalg.inv(H)*np.vstack(xyz)
                    xyz_cryst_new = xyz_cryst + np.matrix([[xcell],[ycell],[zcell]])
                    new_xyz = H*xyz_cryst_new
                    new_xyz = new_xyz.tolist()
                    new_xyz = [new_xyz[0][0], new_xyz[1][0], new_xyz[2][0]]
                    replicated['atoms'].append( atom )
                    replicated['coords'].append( new_xyz )
    return replicated

def replicate_mol(molecule_coords, abc, lmn):
    a,b,c = abc
    replicated = []
    for xcell in range(lmn[0]): # start from 1 so don't duplicate
        for ycell in range(lmn[1]):
            for zcell in range(lmn[2]):
                # make ONE new molecule each time here
                new_mol_coords = []
                for bead in molecule_coords:
                    x, y, z = map(float,bead['xyz'].split())
                    new_xyz = [x + xcell*a, y + ycell*b, z + zcell*c]
                    new_mol_coords.append({'xyz':' '.join(['%e'%k for k in new_xyz]),'q':bead['q']})
                replicated.append(new_mol_coords)
    return replicated

def replicate_file(parent_parser):
    parser = argparse.ArgumentParser(description='replicate struxture file',parents=[parent_parser])
    parser.add_argument('-abc','--boxlengths',help='dimensions of orthorhombic box',
                        type=float,nargs='+')
    parser.add_argument('-f','--file',help='file to replicate',
                        type=str)
    parser.add_argument('-w','--writeFile',help='new file name', type=str,
                        default='replicated.xyz')
    args = vars(parser.parse_args())

    assert len(args['boxlengths']) == 3, '3 dimensions needed for box lengths'

    if '.xyz' in args['file']:
        coordinates = reader.xyz(args['file'])
    elif '.pdb' in args['file']:
        coordinates = reader.PDB(args['file'])
    else:
        raise TypeError('File type not known')

    my_coords = replicate(coordinates, args['boxlengths'] ,args['replicate'])
    writer.PDB(args['writeFile'], my_coords)

def replicate_box(parent_parser):
    import copy
    parser = argparse.ArgumentParser(description='replicate box',parents=[parent_parser])
    parser.add_argument('-b','--box',help='box to replicate',type=str)
    parser.add_argument('-r','--restart',help='restart file',type=str)
    parser.add_argument('-p','--path',help='path to main directories',
                        type=str,default=os.getcwd())
    parser.add_argument('-f','--input',help='input file (fort.4)',type=str)
    parser.add_argument('-w','--extension',help='new file extension', type=str,
                        default='replicated')
    args =  vars(parser.parse_args())
    assert args['input'], parser.print_help()
    input_data = reader.read_fort4(args['input'])
    nmolty = int(input_data['&mc_shared']['nmolty'])
    restart_data = reader.read_restart(args['restart'],
                                       nmolty,
                                       int(input_data['&mc_shared']['nbox']))
    mol_added = {'%i'%i:0 for i in range(1,nmolty+1)}
    box_dimensions = list(map(float,restart_data['box dimensions']['box%s'%args['box']].split()))
    new_restart_data = copy.deepcopy(restart_data)
    new_input_data = copy.deepcopy(input_data)
    for key in 'mol types', 'box types', 'coords':
        new_restart_data[key] = []
    for c, ibox in enumerate(restart_data['box types']):
        imolty = restart_data['mol types'][c]
        if ibox == args['box']:
            mol_added[imolty] -= 1 # dont double count
            # we need to replicate these coordinates
            new_mol = replicate_mol(restart_data['coords'][c], box_dimensions, args['replicate'])
            # keep track of added mols
            for m in range(len(new_mol)):
                new_restart_data['mol types'].append(imolty)
                new_restart_data['box types'].append(ibox)
                new_restart_data['coords'].append(new_mol[m])
                mol_added[imolty] += 1
        else:
            # keep
            for key in 'mol types', 'box types', 'coords':
                new_restart_data[key].append( restart_data[key][c] )
    # tell how many mols added
    total_added = sum(mol_added.values())
    old_nchain = int(input_data['&mc_shared']['nchain'])
    new_nchain = int(new_input_data['&mc_shared']['nchain']) + total_added
    new_input_data['&mc_shared']['nchain'] = input_data['&mc_shared']['nchain'].replace(
        '%i'%old_nchain,'%i'%new_nchain
    )
    new_restart_data['nchain'] = restart_data['nchain'].replace(
        '%i'%old_nchain,'%i'%new_nchain
    )
    for mol, nAdded in mol_added.items():
        nOld = int(new_input_data['SIMULATION_BOX']['box%s'%args['box']]['mol%s'%mol])
        new_input_data['SIMULATION_BOX']['box%s'%args['box']]['mol%s'%mol] = '%i'%(nOld+nAdded)
        molID = input_data['MOLECULE_TYPE']['mol%s'%mol].split('\n')[1].split()[0]
        print('For input dir %s , %i molecules of type %s added '%(
            args['input'][:args['input'].rfind('/')], nAdded,molID)
              )
    for func, inputName, newData, fortNum in zip( [writer.write_restart,writer.write_fort4],
                                            [args['restart'],args['input']],
                                            [new_restart_data, new_input_data],
                                            ['fort.77.','fort.4.']):
        new_file_name = '%s'%args['path']
        if '/' in inputName:
            new_file_name += '/%s'%inputName[:inputName.rfind('/')]
        new_file_name += '/%s%s'%(fortNum, args['extension'])
        func(newData, new_file_name)

from MCFlow.file_formatting import reader
from MCFlow.file_formatting import writer
import os, argparse
import numpy as np

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-t','--type',choices=['replicate_box','replicate_file'])

    parent_parser.add_argument('-lmn','--replicate',help='integer number of cells in each dimension '
                                                  '(1 does no replication)',
                        type=int, nargs='+')
    replicate_box(parent_parser)
