import file_organization as fo


def convertMovieCoordsXYZ(each_molecule):
    conv = {'61': 'H', '62': 'O', '64': 'C', '5': 'C', 'COM': 'F'}
    data = {'atoms': [], 'coords': []}
    for beadType, beads in each_molecule.items():
        for bead_xyz in beads:
            data['atoms'].append(conv[beadType])
            data['coords'].append([float(i) for i in bead_xyz])
    return data


class BeadNotFound(BaseException):
    def __init__(self):
        BaseException.__init__(self)


class Movie:
    def __init__(self, file_name, *args):
        self.file_name = file_name
        self.anal_args = args

    def read_header(self, frac_of_frames=1):
        self.file = open(self.file_name)
        (self.nframes, self.nchain, self.nmolty,
         self.nbox, nBeadType) = [int(i) for i in self.file.readline().split()]
        rcuts = [float(i) for i in self.file.readline().split()]
        typeIDs = [float(i) for i in self.file.readline().split()]
        for mol in range(self.nmolty):  # this nested loops ignores all parameters
            nbeads = int(self.file.readline())
            for paramtype in 'bonds', 'torsions':
                for bead in range(nbeads):
                    self.file.readline()
        self.frame_data = []
        self.boxlengths = []
        self.frame_seed = []
        self.fraction_read = frac_of_frames

    def read_movie_frames(self, seed):
        for frame_count in range(1, self.nframes + 1):
            if (frame_count - 1) % self.fraction_read != 0:
                continue
            try:
                NumPrevCycles = int(self.file.readline())  # num cycles since last frame
            except ValueError:
                self.nframes = frame_count - 1
                break
            self.frame_seed.append(seed)
            Nframe = {};
            boxlxFrame = {};
            FRAME_DATA = {'box%s' % j: {'mol%i' % k: []
                                        for k in range(1, self.nmolty + 1)} for j in range(1, self.nbox + 1)}
            for box in range(1, self.nbox + 1):
                line1 = self.file.readline()
                line2 = self.file.readline()
                # total number of molecules and volume in frame
                Nframe['box%s' % box] = {'mol%i' % j: int(line1.split()[j - 1])
                                         for j in range(1, self.nmolty + 1)}
                boxlxFrame['box%s' % box] = [float(j) for j in line2.split()]
            self.boxlengths.append(copy.deepcopy(boxlxFrame))
            for molnumber in range(1, self.nchain + 1):
                line = self.file.readline()
                (molID, moltype, nunit,
                 cbox, *COM) = [int(j) for j in line.split()[:4]] + [float(j) for j in line.split()[4:]]
                mlcl_info = {'COM': [COM]}
                for bead in range(nunit):
                    line = self.file.readline()
                    xyz = list(map(float, line.split()[:3]))
                    BeadType = line.split()[-1]  # as defined in topmon.inp
                    if BeadType in list(mlcl_info.keys()):
                        mlcl_info[BeadType].append(xyz)
                    else:
                        mlcl_info[BeadType] = [xyz]
                FRAME_DATA['box%s' % cbox]['mol%i' % moltype].append(mlcl_info)
            self.frame_data.append(copy.deepcopy(FRAME_DATA))
        # TODO: assert that all seeds at least have 1 frame
        self.file.close()

    def __add__(self, other):
        '''
        add frame datafrom two movie files
        :param other: other movie file
        '''
        new = self.__class__('masterData', *self.anal_args)
        new.frame_data = self.frame_data + other.frame_data
        new.frame_seed = self.frame_seed + other.frame_seed
        new.nframes = self.nframes + other.nframes
        assert self.nchain == other.nchain, 'nchain not equal for adding movie files'
        new.nchain = self.nchain
        assert self.nmolty == other.nmolty, 'nmolty not equal for adding movie files'
        new.nmolty = self.nmolty
        assert self.nbox == other.nbox, 'nbox not equal for adding movie files'
        new.nbox = self.nbox
        assert self.fraction_read == other.fraction_read, 'fraction read not equal'
        new.fraction_read = self.fraction_read
        new.boxlengths = self.boxlengths + other.boxlengths
        return new

    def filterCoords(self, region, box):
        '''
        :param region: function that determines whether or not in a specific region (bool)
        :param box: box number as string
        :param beadNames: list of bead numbers as strings as in topmon.inp
        '''
        if 'box' in box:
            my_box = box
        else:
            my_box = 'box%s' % box
        for FRAME_DATA in self.frame_data:
            for molType in FRAME_DATA[my_box].keys():
                indices_to_keep = []
                for imol, each_molecule in enumerate(FRAME_DATA[my_box][molType]):
                    if region(convertMovieCoordsXYZ(each_molecule)):
                        indices_to_keep.append(imol)
                FRAME_DATA[my_box][molType] = [value for i, value in enumerate(FRAME_DATA[my_box][molType])
                                               if i in indices_to_keep]

    def countMols(self, indepRange, feed, frame_data):
        '''
        :param box: String box number.
        :param indepRange: range of indep simulations
        '''

        self.averages = {}
        N = {}
        raw_data = {}
        total_frames = -1
        frame_by_seed = [self.frame_seed.count(i) for i in indepRange]
        for FRAME_DATA in frame_data:
            if total_frames == -1:
                N = {box: {mol: {'raw data': [[] for i in indepRange]}
                           for mol in FRAME_DATA[box].keys()}
                     for box in FRAME_DATA.keys()}
                raw_data = {box: {mol: {'raw data': [[] for i in indepRange]}
                                  for mol in FRAME_DATA[box].keys()}
                            for box in FRAME_DATA.keys()}
            # determine which independent simulation this corresponds to
            total_frames += 1
            seed_index = self.frame_seed[total_frames] - 1
            if len(frame_by_seed) == 1: seed_index = 0
            for box in FRAME_DATA.keys():
                for mol in FRAME_DATA[box].keys():
                    if FRAME_DATA[box][mol]:
                        if type(FRAME_DATA[box][mol][0]) == type(1):
                            N[box][mol]['raw data'][seed_index].append(np.mean(FRAME_DATA[box][mol]))
                            raw_data[box][mol]['raw data'][seed_index] += FRAME_DATA[box][mol]
                        elif type(FRAME_DATA[box][mol]) == type([]):
                            #                           print(N.keys(), FRAME_DATA.keys())
                            #                           print(N[box].keys(), FRAME_DATA[box].keys())
                            #                           print(N[box][mol].keys(), FRAME_DATA[box][mol].keys())
                            print(seed_index)
                            N[box][mol]['raw data'][seed_index].append(len(FRAME_DATA[box][mol]))
                            if not raw_data[box][mol]['raw data'][seed_index]:
                                raw_data[box][mol]['raw data'][seed_index] = len(FRAME_DATA[box][mol])
                            else:
                                raw_data[box][mol]['raw data'][seed_index] += len(FRAME_DATA[box][mol])
                        else:
                            print(FRAME_DATA[box][mol])
                            raise NotImplementedError
        #                       if len(N[box][mol]['raw data'][seed_index]) == self.frame_seed.count(seed_index + 1):
        #                           # change list of data to floating point average
        #                           N[box][mol]['raw data'][seed_index] = np.mean(N[box][mol]['raw data'][seed_index])
        #                           if (seed_index > 0) and (type(N[box][mol]['raw data'][seed_index-1]) == type([])):
        #                               print(len(N[box][mol]['raw data'][seed_index-1]))
        #                               print('somethign went wronng')
        for box in N.keys():
            for mol in N[box].keys():
                for i in range(len(N[box][mol]['raw data'])):
                    N[box][mol]['raw data'][i] = np.mean(N[box][mol]['raw data'][i])
        num_molec_data = {}
        for box in N.keys():
            for mol in N[box].keys():
                mol_num = mol.strip('mol')
                if not N[box][mol]:
                    continue
                if mol_num not in num_molec_data.keys():
                    num_molec_data[mol_num] = {}
                if box not in num_molec_data[mol_num].keys():
                    num_molec_data[mol_num][box] = {}
                if len(indepRange) > 1:
                    try:
                        mean, stdev = weighted_avg_and_std(N[box][mol]['raw data'], frame_by_seed)
                    except TypeError:
                        print('no molecules of type %s in %s! -->' % (mol, box), N[box][mol]['raw data'])
                        num_molec_data.pop(mol_num)
                        continue
                else:
                    mean, stdev = np.mean(N[box][mol]['raw data'][0]), np.std(N[box][mol]['raw data'][0])
                num_molec_data[mol_num][box]['mean'] = mean.tolist()
                num_molec_data[mol_num][box]['stdev'] = stdev
                # make histogram
                try:
                    all_data = []
                    for seed_data in raw_data[box][mol]['raw data']:
                        if type(seed_data) == type(1):
                            all_data.append(seed_data)
                        else:
                            all_data += seed_data
                    if len(all_data) > 0:
                        histogram, edges = np.histogram(all_data, bins=list(range(max(all_data) + 2)))
                        num_molec_data[mol_num][box]['histogram'] = histogram.tolist()
                        num_molec_data[mol_num][box]['edges'] = edges.tolist()
                    else:
                        raise ValueError
                except ValueError:
                    if len(num_molec_data[mol_num].keys()) > 1:
                        num_molec_data[mol_num].pop(box)
                    else:
                        num_molec_data.pop(mol_num)
        self.averages[feed] = num_molec_data
        self.averages[feed]['number of frames'] = self.nframes

    def foldMovieToUC(self, uc_vectors):
        '''
        :param box: String box number. Folds all coordinates of all beads
        of all molecules in this box
        '''
        my_box = 'box1'  # no volume moves
        a, b, c = uc_vectors
        for iframe, FRAME_DATA in enumerate(self.frame_data):
            # changing FRAME_DATA will change self.frame_data indexes
            self.boxlengths[iframe][my_box] = uc_vectors
            for molType in FRAME_DATA[my_box].keys():
                for each_molecule in FRAME_DATA[my_box][molType]:
                    for beadType, beadCoords in each_molecule.items():
                        foldedBeads = []
                        for bead_xyz in beadCoords:
                            beadFolded = [float(superCoord) % lattice for (superCoord, lattice) in
                                          zip(bead_xyz, [a, b, c])]
                            foldedBeads.append(beadFolded)
                        each_molecule[beadType] = foldedBeads

    def getCoords(self, mlcl, box, beadNames=['COM']):
        '''
        :param mlcl: String molecule number.
        :param box: String box number.
        of all molecules in this box
        '''
        i, j = math.floor(self.nframes / self.fraction_read), len(self.frame_data)
        assert i == j, 'Error in adding frames %i != %i' % (i, j)
        print('Total amount of frames analyzed was %i' % len(self.frame_data))
        xyz_data = {'atoms': [], 'coords': []}
        for FRAME_DATA in self.frame_data:
            for each_molecule in FRAME_DATA['box%s' % box]['mol%s' % mlcl]:
                for beadType in beadNames:
                    try:
                        for each_coord in each_molecule[beadType]:
                            # for all beads of given bead type in molecule
                            xyz_data['atoms'].append(beadType)
                            xyz_data['coords'].append(each_coord)
                    except KeyError:
                        pass
        #                       print('No beads of type %s found for mol%s'%(beadType, mlcl))
        return xyz_data

    def getAngles(self, uc_vectors, beadOrder=['62', 'COM', '62']):
        '''
        :param mlcl: String molecule number.
        :param box: String box number.
        of all molecules in this box
        '''
        from mcflow.calc_tools import calculate_angle
        angle_histogram = {}
        for iframe, FRAME_DATA in enumerate(self.frame_data):
            for box in FRAME_DATA.keys():
                if box not in angle_histogram.keys(): angle_histogram[box] = {}
                for mlcl, data in FRAME_DATA[box].items():
                    try:
                        for each_molecule in data:
                            angle_coords = {}
                            for bead in set(beadOrder):
                                if (bead not in each_molecule.keys() or
                                        len(each_molecule[bead]) != beadOrder.count(bead)):
                                    raise BeadNotFound
                                angle_coords[bead] = (i for i in each_molecule[bead])
                            coords = []
                            for bead in beadOrder:
                                coords.append([float(i) for i in next(angle_coords[bead])])
                            if box == 'box1':
                                abc = uc_vectors
                            else:
                                abc = self.boxlengths[iframe][box]
                            if mlcl not in angle_histogram[box].keys(): angle_histogram[box][mlcl] = []
                            angle_histogram[box][mlcl].append(calculate_angle(coords[0], coords[1], coords[2], abc))
                    except BeadNotFound:
                        continue
        return angle_histogram

    def getTorsions(self, uc_vectors):
        '''
        :param mlcl: String molecule number.
        :param box: String box number.
        of all molecules in this box
        '''
        from mcflow.calc_tools import get_tors
        torsion_histogram = {}
        for iframe, FRAME_DATA in enumerate(self.frame_data):
            for box in FRAME_DATA.keys():
                if box not in torsion_histogram.keys(): torsion_histogram[box] = {}
                for mlcl, data in FRAME_DATA[box].items():
                    for each_molecule in data:
                        if box == 'box1':
                            abc = uc_vectors
                        else:
                            abc = self.boxlengths[iframe][box]
                        torsion_data = get_tors(each_molecule, abc)
                        if mlcl not in torsion_histogram[box].keys(): torsion_histogram[box][mlcl] = []
                        # combine torsion types
                        all_tors_data = []
                        for torsion_type, tdata in torsion_data.items():
                            all_tors_data += tdata
                        torsion_histogram[box][mlcl].append(
                            all_tors_data)  # do we want to append here, or would it be better to +=?
        return torsion_histogram


def go_through_runs(path, ncycle_total, start_of_runs, num_files, tag='equil-'):
    def initVars(nbox, nmolty):
        chemical_potential = properties.MolProperty(nbox, nmolty)
        number_density = properties.MolProperty(nbox, nmolty)
        volume = properties.Property(nbox)
        return chemical_potential, number_density, volume

    swap_info = {}  # note, swap info will include swatches
    cbmc_info = {}  # includes SAFE-CBMC
    avg_weights = [i / sum(ncycle_total) for i in ncycle_total]
    totalComposition = {}
    zeolite = {}
    runBegin = False
    for j in range(start_of_runs, start_of_runs + num_files):
        f = open(fo.read(path, 'run.', tag, j))
        swap_section = False  # need to find a better way to read the files than this
        swatch_section = False
        cbmc_section = False
        safe_cbmc = False
        for line in f:
            if line.startswith('number of boxes in the system'):
                nBox = int(line.split()[-1])
            elif line.startswith('number of molecule types'):
                nMolTy = int(line.split()[-1])
                if j == start_of_runs:
                    chemical_potential, number_density, volume = initVars(nBox, nMolTy)
            elif line.startswith('number of chains:'):
                nchains = int(line.split()[-1])
            elif ('number of chains of type' in line) and runBegin and j == start_of_runs:
                molName = ' '.join(line.split()[5:-1])
                if molName not in totalComposition.keys():
                    totalComposition[molName] = 0
                if (sum(i for i in totalComposition.values()) < nchains):
                    totalComposition[molName] += int(line.split()[-1])
            elif line.startswith('   temperature:'):
                T = float(line.split()[-2])
            elif 'start of markov chain' in line:
                runBegin = True
            elif line.startswith(' number density      itype'):
                mol = line.split()[3]
                box = int(line.split()[5])
                number_density.data[mol]['box%i' % box].append(float(line.split()[7]))
            elif line.startswith(' chemical potential  itype'):
                mol = line.split()[3]
                box = int(line.split()[5])
                chemical_potential.data[mol]['box%i' % box].append(float(line.split()[7]))
            elif line.startswith('Molecule type, biasing potential 1 through nbox [K]:'):
                biasPot = {}
                for mol in range(1, nMolTy + 1):
                    biasPot[str(mol)] = {}
                    line = f.readline().split()
                    while nBox > len(line):
                        line += f.readline().split()
                    print(len(line), nBox)
                    for box, val in enumerate(line):
                        if '*****' not in val:
                            bp = float(val)
                        else:
                            # set to dummy val
                            bp = 10 ** 5
                        biasPot[str(mol)]['box%i' % (box + 1)] = bp
            elif 'number of unit cells' in line:
                zeolite['unit cells'] = (int(line.split()[-3]) *
                                         int(line.split()[-2]) * int(line.split()[-1]))
            elif 'framework mass' in line:
                zeolite['mass (g)'] = float(line.split()[3])
            elif 'framework volume' in line:
                zeolite['volume [Angst.^3]'] = float(line.split()[3])
            elif 'one adsorbed molecule in sim box' in line:
                zeolite[' mol/kg / 1 mlcl adsorbed'] = float(line.split()[-2])
            elif '###' in line:  # find what section we are in
                if 'Configurational-bias' in line:
                    cbmc_section = True
                elif 'Volume change' in line:
                    cbmc_section = False
                elif 'Molecule swap' in line:
                    swap_section = True
                elif 'Molecule swatch' in line:
                    swatch_section = True
                    swap_section = False
                elif 'Charge Fluctuation' in line:
                    swatch_section = False
                elif 'problem' in line:
                    print('PROBLEM PROBLEM PROBLEM PROBLEM PROBLEM !!!!')
                    print(line)
                    print('file is', f)
                    print('path is', os.getcwd())
                    print(next(f))
            # find what molecule we are looking at
            elif line.startswith('molecule typ') or line.startswith('moltyps'):
                if cbmc_section:
                    # note: not getting box specific information here
                    mol_type = ' '.join(line.split()[3:5])
                    box_type = line.split()[-1]
                    safe_cbmc = False
                elif swap_section:
                    mol_type = ' '.join(line.split()[3:])
                elif swatch_section:
                    mol_type = ' '.join(line.split()[2:])
            # cbmc section analysis
            elif cbmc_section and (len(line) > 2):
                if line.split()[0].isdecimal():
                    if mol_type not in cbmc_info.keys(): cbmc_info[mol_type] = {}  # only want mols that do CBMC
                    if box_type not in cbmc_info[mol_type].keys(): cbmc_info[mol_type][box_type] = {}
                    length = line.split()[0]
                    if safe_cbmc:
                        if 'SAFE' not in cbmc_info[mol_type].keys():
                            cbmc_info[mol_type]['SAFE'] = {length: {'accepted': 0, 'attempted': 0}}
                        elif length not in cbmc_info[mol_type]['SAFE'].keys():
                            cbmc_info[mol_type]['SAFE'][length] = {'accepted': 0, 'attempted': 0}
                        cbmc_info[mol_type]['SAFE'][length]['accepted'] += int(float(line.split()[3]))
                        cbmc_info[mol_type]['SAFE'][length]['attempted'] += int(float(line.split()[1]))
                    else:
                        if length not in cbmc_info[mol_type][box_type].keys():
                            cbmc_info[mol_type][box_type][length] = {'accepted': 0, 'attempted': 0}
                        # the float call here is needed because our code currently
                        # outputs att and acct for CBMC as floats (i.e. 64546525.0)
                        cbmc_info[mol_type][box_type][length]['accepted'] += int(float(line.split()[3]))
                        cbmc_info[mol_type][box_type][length]['attempted'] += int(float(line.split()[1]))
                elif 'SAFE-CBMC' in line:
                    safe_cbmc = True
            # swap section analysis
            elif line.startswith('between box') and swap_section:
                my_attempt, my_accept = int(line.split()[10]), int(line.split()[-1])
                if my_attempt == 0: continue
                if mol_type not in swap_info.keys(): swap_info[mol_type] = {}
                boxpair = line.split()[2] + line.split()[4]  # does not account for box from and box to
                if boxpair not in swap_info[mol_type].keys():
                    swap_info[mol_type][boxpair] = {'accepted': 0, 'attempted': 0}
                swap_info[mol_type][boxpair]['accepted'] += int(line.split()[-1])
                swap_info[mol_type][boxpair]['attempted'] += int(line.split()[10])
            # swatch section analysis
            elif line.startswith('pair typ') and swatch_section:
                swatchPairNum = int(line.split()[-1])
            elif line.startswith('between box') and swatch_section:
                if mol_type not in swap_info.keys(): swap_info[mol_type] = {}
                boxpair = line.split()[2] + line.split()[4]
                if boxpair not in swap_info[mol_type].keys():
                    swap_info[mol_type][boxpair] = {'accepted': 0, 'attempted': 0}
                    swap_info[mol_type][boxpair]['swatchNum'] = swatchPairNum
                swap_info[mol_type][boxpair]['accepted'] += int(line.split()[-1])
                swap_info[mol_type][boxpair]['attempted'] += int(line.split()[10])
            # volume of box section
            elif line.startswith(' boxlength'):
                offset = 2
                boxlx = line
                boxly = next(f)
                boxlz = next(f)
                for boxNum in range(1, nBox + 1):
                    x, y, z = [float(direction.split()[offset + boxNum])
                               for direction in (boxlx, boxly, boxlz)]
                    volume.data['box%i' % boxNum].append([x * y * z])
    if (nchains != sum(i for i in totalComposition.values())):
        print('{} {}'.format(nchains, totalComposition))
    return (number_density.avgOverRuns(avg_weights), chemical_potential.avgOverRuns(avg_weights),
            swap_info, biasPot, volume.avgOverRuns(avg_weights), totalComposition, cbmc_info, T,
            zeolite)


def read_fort12(path, start_of_runs, num_files, tag='equil-'):
    """Read for12 file and output data

    .. warning:
        does not support non-orthorhombic boxes

    """
    os.chdir(path)
    ncycles = [0 for i in range(start_of_runs, start_of_runs + num_files)]
    for j in range(start_of_runs, start_of_runs + num_files):
        run_index = j - start_of_runs
        try:
            f = open(fo.read(path, 'fort12.', tag, j))
        except:
            print(fo.read(path, 'fort12.', tag, j), 'not found')
            quit()
        nline = 0
        line1 = next(f)
        if j == start_of_runs:
            nbox = int(line1.split()[2])
            nmolty = int(line1.split()[3])
            N_mlcls = properties.MolProperty(nbox, nmolty)
            Pressure = properties.Property(nbox)
            InternalEnergy = properties.Property(nbox)
            # Volume = properties.Property(nbox)
            box_length_x = properties.Property(nbox)
            box_length_y = properties.Property(nbox)
            box_length_z = properties.Property(nbox)
            # note that boxlength should only be used for cubic boxes--or should I store as volume?
            molWeights = {}
            mol_num = 0
            for MW in list(map(float, line1.split()[4:])):
                mol_num += 1
                molWeights[str(mol_num)] = MW

        for line in f:  # starting on 2nd line
            nline += 1
            box = nline % nbox
            if box == 0:
                box = nbox
                ncycles[run_index] += 1
            my_split = line.split()
            offset = len(my_split) - nmolty - 1
            for mol in range(1, nmolty + 1):
                N_mlcls.data[str(mol)]['box%i' % box].append(int(my_split[offset + mol]))
            if 'E' in my_split[4]:
                Pressure.data['box%i' % box].append(float(my_split[4]))
            else:
                # pressure has not been calculated yet
                Pressure.data['box%i' % box].append(np.nan)
            InternalEnergy.data['box%i' % box].append(float(my_split[3]))
            box_length_x.data['box%i' % box].append(float(my_split[0]))
            box_length_y.data['box%i' % box].append(float(my_split[1]))
            box_length_z.data['box%i' % box].append(float(my_split[2]))
    return (N_mlcls.data, Pressure.data, box_length_x.data, box_length_y.data,
            box_length_z.data, ncycles, molWeights, InternalEnergy.data)


def read_hbond(path, start_of_runs, num_files):
    n_hBond_oneSim = {j: {'Soln': [], 'Zeo': []} for j in ('Water', 'BuOH')}
    n_hType_oneSim = {mol: {typ: [] for typ in ('Self', 'Other', 'Zeo')} for mol in ('1', '2')}
    for run in range(start_of_runs, num_files + 1):
        movie_nHBond = open(path + '/prod-%i/movie_Hbonds.txt' % run)
        movie_nHBondByType = open(path + '/prod-%i/movie_Hbond_moltype.txt' % run)
        next(movie_nHBond)  # skip first line
        next(movie_nHBondByType)  # skip first line
        for line in movie_nHBond:
            n_hBond_oneSim['Water']['Soln'].append(float(line.split()[1]))
            if line.split()[2] != 'NaN': n_hBond_oneSim['BuOH']['Soln'].append(float(line.split()[2]))
            if line.split()[3] != 'NaN': n_hBond_oneSim['Water']['Zeo'].append(float(line.split()[3]))
            if line.split()[4] != 'NaN': n_hBond_oneSim['BuOH']['Zeo'].append(float(line.split()[4]))
        for line in movie_nHBondByType:
            if 'NaN' not in line:
                mlcl = int(line.split()[1])
                (Hself, Hother, Hzeo) = [float(i) for i in line.split()[2:]]
                n_hType_oneSim[str(mlcl)]['Self'].append(Hself)
                n_hType_oneSim[str(mlcl)]['Other'].append(Hother)
                n_hType_oneSim[str(mlcl)]['Zeo'].append(Hzeo)
    return n_hBond_oneSim, n_hType_oneSim


def getNumMolList(composition, nmolty):
    molTypeCount = []
    for mol in range(1, nmolty + 1):
        for molName in composition.keys():
            if str(mol) in molName:
                molTypeCount.append(composition['molName'])
    return molTypeCount


def read_restart(file, nmolty, nbox):
    config_data = {'max displacement': {}}
    moltyp = 0
    ibox = 0
    charge_moves = []
    nunit = []
    mol_types = []
    box_types = []
    f = open(file)
    nline = 0
    for line in f:
        nline += 1
        if 'number of cycles' not in config_data.keys():
            config_data['number of cycles'] = line
        elif 'atom translation' not in config_data['max displacement'].keys():
            config_data['max displacement']['atom translation'] = line
        elif (moltyp != nmolty) and (ibox != nbox):
            for ibox in range(1, nbox + 1):
                for moltyp in range(1, nmolty + 1):
                    if 'translation' not in config_data['max displacement'].keys():
                        # first time
                        config_data['max displacement']['translation'] = {'box%i' % i: {'mol%i' % j: ''
                                                                                        for j in range(1, nmolty + 1)}
                                                                          for i in range(1, nbox + 1)}
                        config_data['max displacement']['rotation'] = {'box%i' % i: {'mol%i' % j: ''
                                                                                     for j in range(1, nmolty + 1)}
                                                                       for i in range(1, nbox + 1)}
                        config_data['max displacement']['translation']['box%i' % ibox]['mol%i' % moltyp] = line
                        config_data['max displacement']['rotation']['box%i' % ibox]['mol%i' % moltyp] = next(f)
                    elif not config_data['max displacement']['translation']['box%i' % ibox]['mol%i' % moltyp]:
                        config_data['max displacement']['translation']['box%i' % ibox]['mol%i' % moltyp] = next(f)
                        config_data['max displacement']['rotation']['box%i' % ibox]['mol%i' % moltyp] = next(f)
        elif len(charge_moves) < nbox * nmolty:
            if 'fluctuating charge' not in config_data['max displacement'].keys():
                config_data['max displacement']['fluctuating charge'] = {'box%i' % i: {'mol%i' % j: ''
                                                                                       for j in range(1, nmolty + 1)}
                                                                         for i in range(1, nbox + 1)}
            charge_moves = charge_moves + line.split()
            if len(charge_moves) == nbox * nmolty:
                for i in range(len(charge_moves)):
                    c_box = i // nmolty + 1
                    c_mol = i % nmolty + 1
                    config_data['max displacement']['fluctuating charge']['box%i' % c_box]['mol%i' % c_mol] = \
                        charge_moves[i]
            elif len(charge_moves) > nbox * nmolty:
                print('error reading max displ for translation and rotation')
                print(config_data['max displacement'])
                print(nline)
                quit()
        elif 'volume' not in config_data['max displacement'].keys():
            config_data['max displacement']['volume'] = {}
            for index, value in enumerate(line.split()):
                config_data['max displacement']['volume']['box%i' % (index + 1)] = value
            if len(config_data['max displacement']['volume'].keys()) != nbox:
                line = next(f)
                cbox = len(config_data['max displacement']['volume'].keys())
                for index, value in enumerate(line.split()):
                    config_data['max displacement']['volume']['box%i' % (cbox + index + 1)] = value
        elif 'box dimensions' not in config_data.keys():
            if (len(line.split()) == 9):
                config_data['box dimensions'] = {'box1': next(f)}
            else:
                config_data['box dimensions'] = {'box1': line}
            for i in range(2, nbox + 1):
                config_data['box dimensions']['box%i' % i] = next(f)
        elif 'nchain' not in config_data.keys():
            if len(line.split()) == 3:
                box_info = ''
                missing_lines = 5
                for key, value in config_data['box dimensions'].items():
                    box_info += value
                for i in range(missing_lines):
                    box_info += line
                    line = next(f)
                config_data['box dimensions'] = {'box%i' % i: '' for i in range(1, nbox + 1)}
                for c, box_d in enumerate(box_info.split('\n')):
                    if box_d:
                        ibox = c + 1 - missing_lines
                        if ibox < 1: ibox = 1
                        config_data['box dimensions']['box%i' % ibox] += box_d + '\n'
            config_data['nchain'] = line
            nchain = int(line.split()[0])
        elif 'nmolty' not in config_data.keys():
            config_data['nmolty'] = line
            nmolty = int(line.split()[0])
        elif len(nunit) < nmolty:
            nunit += line.split()
            if 'nunit' not in config_data.keys(): config_data['nunit'] = {}
            if len(nunit) == nmolty:
                for i in range(1, len(nunit) + 1):
                    config_data['nunit']['mol%i' % i] = nunit[i - 1]
        elif len(mol_types) < nchain:
            mol_types += line.split()
            if len(mol_types) == nchain:
                config_data['mol types'] = mol_types
        elif len(box_types) < nchain:
            box_types += line.split()
            if len(box_types) == nchain:
                config_data['box types'] = box_types
        else:
            config_data['coords'] = []
            icoords = ''
            for mol_number in config_data['mol types']:
                nbead_mol = int(config_data['nunit']['mol%s' % mol_number])
                molecule_coords = []
                for bead in range(1, nbead_mol + 1):
                    if icoords == '':
                        icoords = line
                    else:
                        icoords = next(f)
                    if len(icoords.split()) == 4:
                        q = icoords.split()[3] + '\n'
                        xyz = ' '.join(icoords.split()[:3]) + ' '
                    else:
                        xyz = icoords
                        q = next(f)
                    molecule_coords.append({'xyz': xyz, 'q': q})
                config_data['coords'].append(molecule_coords)
    return config_data


def read_fort4(file):
    input_data = {}
    f = open(file)
    nmolty = -1
    nbox = -1
    for line in f:
        if len(line.split()) == 0:
            continue
        elif line.split()[0].startswith('&'):
            # in namelist section
            namelist = line.split()[0]
            input_data[namelist] = {}
        elif line.split()[0].startswith('/'):
            namelist = ''
        elif namelist:
            variable = line.split()[0]
            if '=' == variable[-1]: variable = variable.rstrip('=')
            values = line[(line.find('=') + 1):]
            my_val = 'None'
            if (variable == 'pmsatc'):
                my_val = {'swatch%i' % i: values.split()[i - 1]
                          for i in range(1, len(values.split()) + 1)}
            elif variable == 'pmswmt':
                my_val = {'mol%i' % i: values.split()[i - 1] for i in range(1, nmolty + 1)}
            elif (len(values.split()) == 1) and (variable != 'pmvlmt'):
                my_val = values.rstrip('\n')
                if variable == 'nmolty':
                    nmolty = int(values)
                elif variable == 'nbox':
                    nbox = int(values)
            elif (variable == 'pmvlmt'):
                my_val = {'box%i' % i: values.split()[i - 1] for i in range(1, nbox + 1)}
            elif len(values.split()) == nmolty:
                my_val = {'mol%i' % i: values.split()[i - 1] for i in range(1, nmolty + 1)}
            elif len(values.split()) == nbox:
                my_val = {'box%i' % i: values.split()[i - 1] for i in range(1, nbox + 1)}
            elif len(values.split()) > nbox and len(values.split()) < nmolty:
                my_val = {'box%i' % i: values.split()[i - 1] for i in range(1, nbox + 1)}
            elif len(values.split()) > nmolty:
                my_val = {'box%i' % i: values.split()[i - 1] for i in range(1, nbox + 1)}
            else:
                print('error in input file for variable', variable)
                # there was an error in previous input file
                if variable == 'pmswmt':
                    print('there were too many swap types in pmswmt')
                    my_val = {'mol%i' % i: values.split()[i - 1] for i in range(1, nmolty + 1)}
            input_data[namelist][variable] = my_val
        elif (len(line.split()) == 1) and line.split()[0].isupper():
            section = line.split()[0]
            itype = 0
            input_data[section] = {}
        elif line.split()[0] == 'END':
            section = ''
        elif section == 'SIMULATION_BOX':
            while ((itype != nbox) or (not line.split())):
                if 'box%i' % (itype + 1) not in input_data[section].keys():
                    input_data[section]['box%i' % (itype + 1)] = {}
                line = next(f)
                if 'END' in line:
                    print(itype, input_data[section])
                if not line.startswith('!'):
                    if (len(line.split()) == 13) and ('F' in line):
                        (boxlx, boxly, boxlz, rcut, kalp, rcutnn,
                         NDII, lsolid, lrect, lideal, ltwice, T, P) = line.split()
                        input_data[section]['box%i' % (itype + 1)]['dimensions'] = '%s %s %s' % (boxlx, boxly, boxlz)
                        input_data[section]['box%i' % (itype + 1)]['rcut'] = rcut
                        input_data[section]['box%i' % (itype + 1)]['defaults'] = '{} {} {} {} {} {} {}'.format(kalp,
                                                                                                               rcutnn,
                                                                                                               NDII,
                                                                                                               lsolid,
                                                                                                               lrect,
                                                                                                               lideal,
                                                                                                               ltwice)
                        input_data[section]['box%i' % (itype + 1)]['temperature'] = T
                        input_data[section]['box%i' % (itype + 1)]['pressure'] = P
                    elif (len(line.split()) == 9) and ('F' in line):
                        input_data[section]['box%i' % (itype + 1)]['initialization data'] = line
                        itype += 1
                    elif len(line.split()) == nmolty + 1:
                        nmols = line.split()[:-1]
                        nghost = line.split()[-1]
                        for i in range(1, len(nmols) + 1):
                            input_data[section]['box%i' % (itype + 1)]['mol%i' % i] = nmols[i - 1]
                        input_data[section]['box%i' % (itype + 1)]['nghost'] = nghost
                    elif len(line.split()) == nmolty + 2:
                        nmols = line.split()[:-2]
                        nghost = line.split()[-2]
                        for i in range(1, len(nmols) + 1):
                            input_data[section]['box%i' % (itype + 1)]['mol%i' % i] = nmols[i - 1]
                        input_data[section]['box%i' % (itype + 1)]['nghost'] = nghost
        elif section == 'MOLECULE_TYPE':
            if 'nunit' in line:
                itype += 1
                input_data[section]['mol%i' % itype] = ''
            if line.split() and itype > 0:
                input_data[section]['mol%i' % itype] += line
        elif section == 'MC_SWAP':
            if line.split() and line[0].isdigit():
                if len([i for i in line.split() if i.isdigit()]) == 2:
                    # line is box1, box2
                    swap_pair += 1
                    input_data[section]['mol%i' % itype]['box1 box2'][swap_pair] = list(map(int, line.split()))
                else:
                    itype += 1
                    input_data[section]['mol%i' % itype] = {'nswapb': 0, 'pmswapb': []}
                    nswapb = int(line.split()[0])
                    input_data[section]['mol%i' % itype]['nswapb'] = nswapb
                    input_data[section]['mol%i' % itype]['pmswapb'] = list(
                        map(float, [i.rstrip('d0') for i in line.split()[1:]]))
                    input_data[section]['mol%i' % itype]['box1 box2'] = [0 for i in range(nswapb)]
                    swap_pair = -1
        elif section == 'MC_SWATCH':
            if 'moltyp1<->moltyp2' in line:
                itype += 1
                input_data[section]['swatch%i' % itype] = ''
            if (itype > 0) and line.split():
                input_data[section]['swatch%i' % itype] += line
        elif section == 'INTERMOLECULAR_EXCLUSION':
            if (line.rstrip('\n').replace(' ', '').isdigit()):
                mol1, unit1, mol2, unit2 = map(int, line.split())
                if 'mol%i' % mol1 not in input_data[section].keys():
                    input_data[section]['mol%i' % mol1] = {}
                if 'unit%i' % unit1 not in input_data[section]['mol%i' % mol1].keys():
                    input_data[section]['mol%i' % mol1]['unit%i' % unit1] = {}
                if 'mol%i' % mol2 not in input_data[section]['mol%i' % mol1]['unit%i' % unit1].keys():
                    input_data[section]['mol%i' % mol1]['unit%i' % unit1]['mol%i' % mol2] = []
                input_data[section]['mol%i' % mol1]['unit%i' % unit1]['mol%i' % mol2].append('unit%i' % unit2)
        elif section == 'INTRAMOLECULAR_OH15':
            if (line.rstrip('\n').replace(' ', '').isdigit()):
                mol = int(line.split()[0])
                if 'mol%i' % mol not in input_data[section].keys():
                    input_data[section]['mol%i' % mol] = []
                input_data[section]['mol%i' % mol] += [line[line.find(' '):]]
        elif section == 'INTRAMOLECULAR_SPECIAL':
            params = line.split()
            if params[0].isdigit():
                mol, i, j, logic, sLJ, sQ = list(map(int, params[:4])) + list(map(float, params[-2:]))
                if 'mol%i' % mol not in input_data[section].keys():
                    input_data[section]['mol%i' % mol] = []
                input_data[section]['mol%i' % mol].append('%i %i %i %2.1f %2.1f' % (i, j, logic, sLJ, sQ))
        elif section == 'UNIFORM_BIASING_POTENTIALS':
            if '!' not in line:
                itype += 1
                mol = 'mol%i' % itype
                if mol not in input_data[section].keys():
                    input_data[section][mol] = {}
                for i in range(len(line.split())):
                    my_box = 'box%i' % (i + 1)
                    if my_box not in input_data[section][mol].keys():
                        input_data[section][mol][my_box] = {}
                    new_var = line.split()[i]
                    if 'd0' in new_var: new_var = new_var.rstrip('d0')
                    input_data[section][mol][my_box] = float(new_var)
        else:
            print(section, 'is missing formatting')
    return input_data


def PDB(file):
    data = {'coords': [], 'atoms': []}
    for line in open(file):
        if line.startswith('CRYST1'):
            a, b, c, alpha, beta, gamma = map(float, line.split()[1:7])
            data['box info'] = {'a': a, 'b': b, 'c': c,
                                'alpha': alpha, 'beta': beta, 'gamma': gamma}
        elif line.startswith('ATOM'):
            atom = line.split()[2]
            x, y, z = map(float, line.split()[5:8])
            data['atoms'].append(atom)
            data['coords'].append([x, y, z])
    return data


def xyz(file):
    data = {'atoms': [], 'coords': []}
    for line in open(file):
        if len(line.split()) == 4:
            x, y, z = map(float, line.split()[1:])
            atom = line.split()[0]
            data['atoms'].append(atom)
            data['coords'].append([x, y, z])
        elif len(line.split()) == 9:
            data['hmatrix'] = list(map(float, line.split()))
    return data


import os
import copy
import numpy as np
import properties
from mcflow.calc_tools import weighted_avg_and_std
import math
