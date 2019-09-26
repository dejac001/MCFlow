'''
Add molecules to a fort.4 and restart file.
You will need to use this, for example, if
too many of a given sorbate adsorb and not
enough molecules are left in vapor phase.
'''
def findMolNum(mol_type_data, molName):
    for mol in sorted(mol_type_data.keys()):
        my_name = mol_type_data[mol].split('\n')[1].split()[0]
        if molName == my_name:
            return mol.strip('mol')
    assert molName == my_name, 'Molecule type %s not found'%molName


def makeRandomStruc(boxlengths, coordinates, previous_coords):
    boxlx, boxly, boxlz = boxlengths
    rand_x = random.random()*boxlx*0.80 + 0.10*boxlx
    rand_y = random.random()*boxly*0.80 + 0.10*boxly
    rand_z = random.random()*boxlz*0.80 + 0.10*boxlz
    # place molecule from first bead and fold back into box
    first_bead = [rand_x, rand_y, rand_z]
    coords = [first_bead]
    for ii in range(1,len(coordinates)):
        xyz = coordinates[ii]
        vector = [coordinates[0][i] - xyz[i] for i in range(len(xyz))]
        coords.append(
            [first_bead[i] - vector[i] for i in range(3)]
                       )
    # test that vector is the same
    for bead in range(1,len(coordinates)):
        for dir in range(3):
            dist_old = coordinates[bead][dir] - coordinates[0][dir]
            dist_new = coords[bead][dir] - coords[0][dir]
            assert abs(dist_old-dist_new) < 1e-6, 'incorrect struc added'

    # test if there is an overlap
    for xyz in coords:
        for xyz_old in previous_coords:
            dist = calculate_distance(xyz, xyz_old, [boxlx, boxly, boxlz])
            if dist < 1.0:
                print('distance less than 1.0, using recursion')
                return makeRandomStruc(boxlengths, coordinates, previous_coords)
    return coords

def initialize(input, restart,  molID, box):

    try:
        boxlx, boxly, boxlz = map(float,restart['box dimensions']['box%s'%box].split())
    except ValueError:
        print(restart['box dimensions']['box%s'%box])
        print('non-orthorhombic box not work')
        quit()

    # find mol number
    mol_number = findMolNum(input['MOLECULE_TYPE'], molID)

    # get list of coordinates previously in box
    old_coords = []
    for i, ibox in enumerate(restart['box types']):
        if ibox == box:
            for bead in restart['coords'][i]:
                old_coords.append( list(map(float,bead['xyz'].split())) )

    # get charges from previous molecule
    q = []
    for bead in restart['coords'][restart['mol types'].index(mol_number)]:
        q.append( float(bead['q'].rstrip('\n')) )
    return (boxlx, boxly, boxlz), mol_number, old_coords, q

def addMolecules(input_dat, restart_dat, nAddorRemove, box, molID, changeVol):
    def getXYZCoords(mol_coords):
        coords = []
        for bead in mol_coords:
            coords.append(list(map(float,bead['xyz'].split())))
        return coords
    # initialize volume to ideal gas volume in new box
    p = float(input_dat['SIMULATION_BOX']['box%s'%box]['pressure'])
    N = restart_dat['box types'].count(box) + nAddorRemove
    temperature = input_dat['SIMULATION_BOX']['box%s'%box]['temperature']
    if 'd0' in temperature:
        temperature = temperature.rstrip('d0')
    T = float(temperature)
    V = N/N_av*R['\AA**3*MPa/(mol*K)']*T/p
    boxlx = pow(V, 1/3)
    boxlx_old = next(map(float,restart_dat['box dimensions']['box%s'%box].split()))
    if (boxlx > boxlx_old) and changeVol == 'True':
        restart_dat['box dimensions']['box%s'%box] = '{} {} {}\n'.format(boxlx, boxlx, boxlx)
        input_dat['SIMULATION_BOX']['box%s'%box]['rcut'] = '%e'%(boxlx/2*0.99)
    input_dat['&mc_shared']['iratio'] = '500'
    input_dat['&mc_shared']['time_limit'] = '1'
    if '&analysis' in input_dat.keys():
        input_dat['&analysis']['imv'] = '%i'%(int(input_dat['&mc_shared']['nstep']) + 10)
    input_dat['&mc_volume']['iratv'] = '500'
    if float(restart_dat['max displacement']['volume']['box%s'%box]) < 100000.0:
        restart_dat['max displacement']['volume']['box%s'%box] = '100000.0'
    (boxlengths, mol_num, old_coordinates, charges) = initialize(input_dat, restart_dat, molID, box)
    # get info of old structures of molecules to choose from in making new structures
    mol_coord_data = {}
    for key in ['box types','mol types', 'coords']:
        mol_coord_data[key] = []
    for my_mol, my_box, my_coords in zip(restart_dat['mol types'], restart_dat['box types'],
                                         restart_dat['coords']):
        if (my_mol == mol_num):
            for key, value in zip(['mol types','box types','coords'], [my_mol, my_box, my_coords]):
                mol_coord_data[key].append(value)
    assert len(mol_coord_data['mol types']) > 0, 'No mols found for moltype in boxtype provided'
    for newMol in range(nAddorRemove):
        try:
            if (newMol+1)%(nAddorRemove//4) == 0:
                print('%5.1f %% of molecules added'%(100*(newMol+1)/nAddorRemove))
        except ZeroDivisionError:
            print('%5.1f %% of molecules added'%(100*(newMol+1)/nAddorRemove))
        # find random molecule to copy configuration of
        random_mol = random.randint(0, len(mol_coord_data['mol types'])-1)
        new_config = getXYZCoords(mol_coord_data['coords'][random_mol])
        #  get new coordinates
        my_coordinates = makeRandomStruc(boxlengths, new_config, old_coordinates)
        old_coordinates += my_coordinates
        # add to fort.77 file
        restart_dat['mol types'].append( mol_num )
        restart_dat['box types'].append( box )
        new_coordinates = []
        for i in range(len(my_coordinates)):
            new_coordinates.append(
                {'xyz': ' '.join(['%16.8f'%j for j in my_coordinates[i]]) + '\n',
                 'q': '%f\n'%charges[i]}
            )
        restart_dat['coords'].append( new_coordinates )
        restart_dat['nchain'] = restart_dat['nchain'].replace(
            restart_dat['nchain'].split()[0], '%i'%(int(restart_dat['nchain'].split()[0])+1)
        )
        input_dat['SIMULATION_BOX']['box%s'%box]['mol%s'%mol_num] = '%i'%(
            int(input_dat['SIMULATION_BOX']['box%s'%box]['mol%s'%mol_num]) + 1
        )
    input_dat['&mc_shared']['nchain'] = restart_dat['nchain'].split()[0]
    return copy.deepcopy(input_dat), copy.deepcopy(restart_dat)

def removeMolecules(input_dat, restart_dat, nAddorRemove, box, molID):
    def keepMol():
        for key in ['box types', 'mol types', 'coords']:
            new_restart_data[key].append(restart_dat[key][i])
    # find mol number
    mol_num = findMolNum(input_dat['MOLECULE_TYPE'], molID)
    taken_out = 0
    new_restart_data = copy.deepcopy(restart_dat)
    # initialize new restart data
    for key in ['box types', 'mol types', 'coords']:
        new_restart_data[key] = []
    for i in range(len(restart_dat['box types'])):
        if (restart_dat['mol types'][i] == mol_num) and (restart_dat['box types'][i] == box):
            if taken_out > nAddorRemove:
                # don't keep
                taken_out -= 1
            elif taken_out == nAddorRemove:
                # correct molec but have already taken enough out
                keepMol()
        else:
            keepMol()
    assert taken_out == nAddorRemove, 'Not enough mols of type {} taken out; only {}'.format(mol_num, taken_out)

    new_restart_data['nchain']  = restart_dat['nchain'].replace(
                restart_dat['nchain'].split()[0], '%i'%(int(restart_dat['nchain'].split()[0])+taken_out)
                )
    input_dat['&mc_shared']['nchain'] = new_restart_data['nchain'].split()[0]
    assert int(new_restart_data['nchain']) ==  int(input_dat['&mc_shared']['nchain']), 'Conflicting info in rest & control files'
    input_dat['SIMULATION_BOX']['box%s'%box]['mol%s'%mol_num] = '%i'%(
        int(input_dat['SIMULATION_BOX']['box%s'%box]['mol%s'%mol_num]) + taken_out
    )
    return copy.deepcopy(input_dat), copy.deepcopy(new_restart_data)

def removeChains(input_dat, restart_dat, chains_to_remove):
    taken_out = 0
    new_restart_data = copy.deepcopy(restart_dat)
    # initialize new restart data
    for key in ['box types', 'mol types', 'coords']:
        new_restart_data[key] = []
    for ichain in range(1,len(restart_dat['box types'])+1):
        if ichain in chains_to_remove:
            # don't keep
            molty = 'mol%s'%restart_dat['mol types'][ichain-1]
            boxty = 'box%s'%restart_dat['box types'][ichain-1]
            input_dat['SIMULATION_BOX'][boxty][molty] = '%i'%(
                int(input_dat['SIMULATION_BOX'][boxty][molty]) -1
            )
            taken_out -= 1
        else:
            for key in ['box types', 'mol types', 'coords']:
                new_restart_data[key].append(restart_dat[key][ichain-1])
    new_restart_data['nchain']  = restart_dat['nchain'].replace(
                restart_dat['nchain'].split()[0], '%i'%(int(restart_dat['nchain'].split()[0])+taken_out)
                )
    input_dat['&mc_shared']['nchain'] = new_restart_data['nchain'].split()[0]
    return copy.deepcopy(input_dat), copy.deepcopy(new_restart_data)

from MCFlow.file_formatting import reader, writer
from MCFlow.calc_tools import fold, calculate_distance
from MCFlow.chem_constants import N_av, R
import random, copy
import numpy as np

if __name__ == '__main__':
    from analysis_parsers import ChangeInput
    my_parser = ChangeInput()
    my_parser.molecules()
    my_parser.parser.add_argument('-cv','--changeVol',help='whether to change volume',type=str,default='False',choices=['True','False'])
    args = vars(my_parser.parse_args())
    assert args['nAddorRemove'] != 0, 'Cannot add or remove 0 molecules'

    for feed in args['feeds']:
        for seed in args['indep']:
            if (feed == '.'):
                base_dir = ''
            else:
                base_dir = '%s/%s/%i/'%(args['path'],feed,seed)
            input_data = reader.read_fort4('%s%s'%(base_dir,args['input']))
            nmolty, nbox = (int(input_data['&mc_shared']['nmolty']),
                    int(input_data['&mc_shared']['nbox']))
            restart_data = reader.read_restart('%s%s'%(base_dir,args['restart']),nmolty, nbox)
            if args['boxAdd']:
                new_input_data, new_restart_data = addMolecules(input_data, restart_data,
                                                args['nAddorRemove'], args['boxAdd'],args['molID'], args['changeVol'])
                if args['boxRemove']:
                    new_input_data, new_restart_data = removeMolecules(new_input_data, new_restart_data,
                        -1*args['nAddorRemove'], args['boxRemove'],args['molID'])
            elif args['boxRemove']:
                new_input_data, new_restart_data = removeMolecules(input_data, restart_data,
                        -1*args['nAddorRemove'], args['boxRemove'],args['molID'])
            elif len(args['chainsToRemove']) > 0:
                new_input_data, new_restart_data = removeChains(input_data, restart_data, args['chainsToRemove'])
            writer.write_fort4(new_input_data, '%sfort.4.newMols'%base_dir)
            writer.write_restart(new_restart_data, '%sfort.77.newMols'%base_dir)
