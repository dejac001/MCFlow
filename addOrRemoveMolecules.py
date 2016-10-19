'''
Add molecules to a fort.4 and restart file.
You will need to use this, for example, if
too many of a given sorbate adsorb and not
enough molecules are left in vapor phase.
'''
def fold(xyz, boxlengths):
    def pbc(coord, boxlength):
        if coord > boxlength:
            folded = coord - boxlength
        elif coord < 0:
            folded = coord + boxlength
        else:
            folded = coord
        return folded
    new_coords = []
    for i, boxl in zip(xyz, boxlengths):
        new_coords.append( pbc(i,boxl) )
    return new_coords

def calculate_distance(xyz1, xyz2, abc):
    vector = [xyz1[i] - xyz2[i] for i in range(len(xyz1))]
    for i in range(len(abc)):
        # implement minimum image
        if vector[i] > abc[i]/2.:
            # positive
            vector[i] = abc[i] - vector[i]
        elif vector[i] < -1.*abc[i]/2.:
            vector[i] = abc[i] + vector[i]
    return np.linalg.norm(vector, 2, 0)

def makeRandomStruc(boxlengths, coordinates, previous_coords):
    boxlx, boxly, boxlz = boxlengths
    rand_x = random.random()*boxlx
    rand_y = random.random()*boxly
    rand_z = random.random()*boxlz
    # place molecule from first bead and fold back into box
    coords = []
    for xyz in coordinates:
        coords.append(
            fold([xyz[0] + rand_x, xyz[1]+rand_y, xyz[2]+rand_z],
                 [boxlx, boxly, boxlz])
                       )
    # test if there is an overlap
    for xyz in coords:
        for xyz_old in previous_coords:
            dist = calculate_distance(xyz, xyz_old, [boxlx, boxly, boxlz])
            if dist < 1.0:
                print('distance less than 1.0, using recursion')
                return makeRandomStruc(boxlengths, coordinates, previous_coords)
    return coords

def initialize(path, feed, seed, restart_name, input_name,  molID, box):
    input_data = reader.read_fort4('%s/%s/%i/%s'%(path,feed,seed,input_name))
    nmolty, nbox = (int(input_data['&mc_shared']['nmolty']),
                    int(input_data['&mc_shared']['nbox']))
    restart_data = reader.read_restart('%s/%s/%i/%s'%(path,feed,seed,restart_name),nmolty, nbox)

    boxlx, boxly, boxlz = map(float,restart_data['box dimensions']['box%s'%box].split())

    # find mol number
    for mol in input_data['MOLECULE_TYPE'].keys():
        if molID in input_data['MOLECULE_TYPE'][mol]:
            mol_number = mol.strip('mol')

    # get list of coordinates previously in box
    old_coords = []
    for i, ibox in enumerate(restart_data['box types']):
        if ibox == box:
            for bead in restart_data['coords'][i]:
                old_coords.append( list(map(float,bead['xyz'].split())) )

    # get charges from previous molecule
    q = []
    for bead in restart_data['coords'][restart_data['mol types'].index(mol_number)]:
        q.append( float(bead['q'].rstrip('\n')) )
    return restart_data, input_data, (boxlx, boxly, boxlz), mol_number, old_coords, q

def addMolecules(feed, seed, path, nAdd, box, restart, input, molID, **kwargs):
    def getXYZCoords(mol_coords):
        coords = []
        for bead in mol_coords:
            coords.append(list(map(float,bead['xyz'].split())))
        return coords
    (restart_data, input_data, boxlengths,
     mol_num, old_coordinates, charges) = initialize(path, feed, seed, restart,input, molID, box)
    # get info of old structures of molecules to choose from in making new structures
    mol_coord_data = {}
    for key in ['box types','mol types', 'coords']:
        mol_coord_data[key] = []
    for my_mol, my_box, my_coords in zip(restart_data['mol types'], restart_data['box types'],
                                         restart_data['coords']):
        if (my_mol == mol_num) and (my_box == box):
            for key, value in zip(['mol types','box types','coords'], [my_mol, my_box, my_coords]):
                mol_coord_data[key].append(value)
    for newMol in range(nAdd):
        # find random molecule to copy configuration of
        random_mol = random.randint(0, len(mol_coord_data['mol types'])-1)
        new_config = getXYZCoords(mol_coord_data['coords'][random_mol])
        #  get new coordinates
        my_coordinates = makeRandomStruc(boxlengths, new_config, old_coordinates)
        # add to fort.77 file
        restart_data['mol types'].append( mol_num )
        restart_data['box types'].append( box )
        new_coordinates = []
        for i in range(len(my_coordinates)):
            new_coordinates.append(
                {'xyz': ' '.join(['%f'%j for j in my_coordinates[i]]) + '\n',
                 'q': '%f\n'%charges[i]}
            )
        restart_data['coords'].append( new_coordinates )
        restart_data['nchain'] = restart_data['nchain'].replace(
            restart_data['nchain'].split()[0], '%i'%(int(restart_data['nchain'].split()[0])+1)
        )
        input_data['&mc_shared']['nchain'] = restart_data['nchain'].split()[0]
        input_data['SIMULATION_BOX']['box%s'%box]['mol%s'%mol_num] = '%i'%(
            int(input_data['SIMULATION_BOX']['box%s'%box]['mol%s'%mol_num]) + 1
        )
    return copy.deepcopy(restart_data), copy.deepcopy(input_data)

def removeMolecules(feed, seed, path, box, nAdd,restart, input, molID, **kwargs):
    def keepMol():
        for key in ['box types', 'mol types', 'coords']:
            new_restart_data[key].append(restart_data[key][i])
    (restart_data, input_data, boxlengths,
     mol_num, old_coordinates, charges) = initialize(path, feed, seed, restart,input, molID, box)
    taken_out = 0
    new_restart_data = copy.deepcopy(restart_data)
    # initialize new restart data
    for key in ['box types', 'mol types', 'coords']:
        new_restart_data[key] = []
    for i in range(len(restart_data['box types'])):
        if (restart_data['mol types'][i] == mol_num) and (restart_data['box types'][i] == box):
            # don't keep
            if taken_out > nAdd:
                taken_out -= 1
            elif taken_out == nAdd:
                keepMol()
        else:
            keepMol()

    new_restart_data['nchain']  = restart_data['nchain'].replace(
                restart_data['nchain'].split()[0], '%i'%(int(restart_data['nchain'].split()[0])+taken_out)
                )
    input_data['&mc_shared']['nchain'] = restart_data['nchain'].split()[0]
    input_data['SIMULATION_BOX']['box%s'%box]['mol%s'%mol_num] = '%i'%(
        int(input_data['SIMULATION_BOX']['box%s'%box]['mol%s'%mol_num]) + taken_out
    )
    return new_restart_data, copy.deepcopy(input_data)


from file_formatting import reader, writer
import random, copy
import numpy as np

if __name__ == '__main__':
    from parser import ChangeInput
    my_parser = ChangeInput()
    my_parser.molecules()
    args = vars(my_parser.parse_args())
    assert args['nAdd'] != 0, 'Cannot add or remove 0 molecules'

    for feed in args['feeds']:
        for seed in args['indep']:
            if args['nAdd'] > 0:
                new_restart_data, new_input_data = addMolecules(feed, seed, **args)
            elif args['nAdd'] < 0:
                new_restart_data, new_input_data = removeMolecules(feed, seed, **args)
            writer.write_fort4(new_input_data, '%s/%s/%i/fort.4.newMols'%(args['path'],feed,seed))
            writer.write_restart(new_restart_data, '%s/%s/%i/fort.77.newMols'%(args['path'],feed,seed))
