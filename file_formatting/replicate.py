def replicate(coordinates, abc, lmn):
    '''

    :param coordinates: data type read in by xyz or pdb
    :param abc: box dimensions
    :param lmn: integer number of cells in each dimension
    :return:
    '''
    a, b, c = abc
    replicated = {'atoms':[], 'coords':[]}
    for xcell in range(lmn[0]):
        for ycell in range(lmn[1]):
            for zcell in range(lmn[2]):
                for atom, xyz in zip(coordinates['atoms'], coordinates['coords']):
                    x, y, z = xyz
                    new_xyz = [x + xcell*a, y + ycell*b, z + zcell*c]
                    replicated['atoms'].append( atom )
                    replicated['coords'].append( new_xyz )
    return replicated

def replicate_mol(molecule_coords, abc, lmn):
    a,b,c = abc
    replicated = []
    for xcell in range(lmn[0]):
        for ycell in range(lmn[1]):
            for zcell in range(lmn[2]):
                # make ONE new molecule each time here
                new_mol_coords = []
                for bead in molecule_coords:
                    x, y, z = bead['xyz']
                    new_xyz = [x + xcell*a, y + ycell*b, z + zcell*c]
                    new_mol_coords.append({'xyz':new_xyz,'q':bead['q']})
                replicated.append(new_mol_coords)
    return replicated

def replicate_file(parser):
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
    writer.xyz(args['writeFile'], my_coords)

def replicate_box(parser):
    import copy
    parser.add_argument('-b','--box',help='box to replicate',type=str)
    parser.add_argument('-r','--restart',help='restart file',type=str)
    parser.add_argument('-p','--path',help='path to main directories',
                        type=str,default=os.getcwd())
    parser.add_argument('-f','--input',help='input file (fort.4)',type=str)
    parser.add_argument('-w','--extension',help='new file extension', type=str,
                        default='replicated')
    args =  vars(parser.parse_args())
    input_data = reader.read_fort4(args['input'])
    nmolty = int(input_data['&mc_shared']['nmolty'])
    restart_data = reader.read_restart(args['restart'],
                                       nmolty,
                                       int(input_data['&mc_shared']['nbox']))
    mol_added = {'%i'%i:0 for i in range(1,nmolty+1)}
    box_dimensions = list(map(float,restart_data['box dimensions']['box%s'%args['box']]))
    new_restart_data = copy.deepcopy(restart_data)
    new_input_data = copy.deepcopy(input_data)
    for c, ibox in enumerate(restart_data['box types']):
        if ibox == args['box']:
            imolty = restart_data['mol types'][c]
            # we need to replicate these coordinates
            new_mol = replicate_mol(restart_data['coords'][c], box_dimensions, args['replicate'])
            # keep track of added mols
            for m in range(len(new_mol)):
                new_restart_data['mol types'].append(imolty)
                new_restart_data['box types'].append(ibox)
                new_restart_data['coords'].append(new_mol[m])
                mol_added[imolty] += 1
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
        print('For input dir %s , %i molecules of type %s added '%(
            args['input'][:args['input']].rfind('/'), nAdded,
                        input_data['MOLECULE_TYPE']['mol%s'%mol].split()[0])
              )
    writer.write_restart(new_restart_data, '%s/%s/fort.77.%s'%(
        args['path'],args['restart'][:args['restart'].rfind('/')],args['extension']))
    writer.write_fort4(new_input_data, '%s/%s/fort.4.%s'%(
        args['path'],args['input'][:args['input'].rfind('/')],args['extension']))

from MCFlow.file_formatting import reader
from MCFlow.file_formatting import writer
import os

if __name__ == '__main__':
    import argparse
    parent_parser = argparse.ArgumentParser(description='replicate xyz or pdb file or box')
    parent_parser.add_argument('-t','--type',choices=['replicate_box','replicate_file'])

    parent_parser.add_argument('-lmn','--replicate',help='integer number of cells in each dimension '
                                                  '(1 does no replication)',
                        type=int, nargs='+')
    parent_args = vars(parent_parser.parse_args())

    assert args['file'], parent_parser.print_help()
    if parent_args['type'] == 'replicate_box':
        replicate_box(parent_parser)
    elif parent_args['type'] == 'replicate_file':
        replicate_file(parent_parser)