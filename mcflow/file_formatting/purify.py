def removeExtraInfo(boxesToKeep,restartData,inputData):
    '''
    remove extra input information from an input file.
    Specify certain boxes to keep, keep only molecules in those boxes
    :param boxesToKeep: list of boxes to keep in as str
    :param restartData: restart data as returned from reader.read_restart
    :param inputData: input data as returned from reader.read_fort4
    :return: newRestartData, newInputData
    '''
    import copy
    newRestartData = restartData.copy() # don't want to change old stuff, could be bad
    # from restart data, find molecules to keep
    moleculesToKeep = []
    nchain_old = int(restartData['nchain'].split()[0])
    nchain_new = 0
    nbox_old = len(restartData['max displacement']['translation'].keys())
    nmolty_old = int(restartData['nmolty'].split()[0])
    nmolty_by_box = {}
    for i, box in enumerate(restartData['box types']):
        if box in boxesToKeep:
            mol = restartData['mol types'][i]
            old_box = 'box%s'%box
            if old_box not in nmolty_by_box.keys():
                nmolty_by_box[old_box] = {}
            if mol not in moleculesToKeep:
                moleculesToKeep.append(mol)
                new_mol = 'mol%i'%len(moleculesToKeep)
                nmolty_by_box[old_box][new_mol] = 0
            elif new_mol not in nmolty_by_box[old_box].keys():
                nmolty_by_box[old_box][new_mol] = 0
            nmolty_by_box[old_box][new_mol] += 1
            nchain_new += 1
    for box in ['box%s'%i for i in boxesToKeep]:
        if box not in nmolty_by_box.keys():
            nmolty_by_box[box] = {}
        for mol in ['mol%i'%i for i in range(1,len(moleculesToKeep)+1)]:
            if mol not in nmolty_by_box[box].keys():
                nmolty_by_box[box][mol] = 0
    moleculesToKeep = list(map(str,sorted(map(int,moleculesToKeep))))
    if not moleculesToKeep:
        print('no molecules predicted to be kept')
        print('purify: moleculesToKeep {} for boxes {}'.format(moleculesToKeep, boxesToKeep))
        quit()
    #   # remove extra boxes and molecules from restart data----------------------------------------------
    # only keep relevant data from maximum diplacements
    for disp_type, values in restartData['max displacement'].items():
        if disp_type != 'atom translation':
            newRestartData['max displacement'][disp_type] = {}
            for box in sort_keys(values.keys()):
                if box.strip('box') in boxesToKeep:
                    new_box = 'box%i'%(boxesToKeep.index(box.strip('box')) + 1)
                    if type(values[box]) != type({}):
                        newRestartData['max displacement'][disp_type][new_box] = values[box]
                    else:
                        newRestartData['max displacement'][disp_type][new_box] = {}
                        for mol in sort_keys(values[box].keys()):
                            if mol.strip('mol') in moleculesToKeep:
                                new_mol = 'mol%i'%(moleculesToKeep.index(mol.strip('mol')) + 1)
                                newRestartData['max displacement'][disp_type][new_box][new_mol] = values[box][mol]
    newRestartData['box dimensions'] = {}
    for box in restartData['box dimensions'].keys():
        if box.strip('box') in boxesToKeep:
            new_box = 'box%i'%(boxesToKeep.index(box.strip('box')) + 1)
            newRestartData['box dimensions'][new_box] = restartData['box dimensions'][box]
    #  update number of chains
    newRestartData['nchain'] = restartData['nchain'].replace('%i'%nchain_old, '%i'%nchain_new)
    # update number of molecule types
    newRestartData['nmolty'] = restartData['nmolty'].replace('%i'%nmolty_old, '%i'%len(moleculesToKeep))
    # remove extra mols from nunit
    newRestartData['nunit'] = {}
    new_mol_number = 0
    for mol in sort_keys(restartData['nunit'].keys()):
        if mol.strip('mol') in moleculesToKeep:
            new_mol_number += 1
            newRestartData['nunit']['mol%i'%new_mol_number] = restartData['nunit'][mol]
    # remove coordinates from restart file
    for dataType in ['box types','mol types','coords']:
        newRestartData[dataType] = []
    for i, box in enumerate(restartData['box types']):
        if box in boxesToKeep:
            newRestartData['coords'].append(restartData['coords'][i])
            # change mol index
            mol = restartData['mol types'][i]
            newRestartData['mol types'].append('%i'%(moleculesToKeep.index(mol)+1))
            # change box index
            newRestartData['box types'].append('%i'%(boxesToKeep.index(box)+1))
            
            
    if len(newRestartData['mol types']) != nchain_new:
        print('problem removing molecules')
        quit()
    #   # remove extra boxes and molecules from input data----------------------------------------------
    newInputData = copy.deepcopy(inputData)
    newInputData['&mc_shared']['nchain'] = inputData['&mc_shared']['nchain'].replace('%i'%nchain_old,'%i'%nchain_new)
    newInputData['&mc_shared']['nbox'] = inputData['&mc_shared']['nbox'].replace('%i'%nbox_old,'%i'%len(boxesToKeep))
    newInputData['&mc_shared']['nmolty'] = inputData['&mc_shared']['nmolty'].replace('%i'%nmolty_old,'%i'%len(moleculesToKeep))
    namelists = [i for i in inputData.keys() if i.startswith('&')]
    for namelist in namelists:
        for variable,value in inputData[namelist].items():
            if (type(value) == type({})) and (variable != 'pmsatc'):
                if len(value.keys()) == nmolty_old:
                    newInputData[namelist][variable] = {}
                    # initialize with no molecules and only put in those needed
                    for mol in sort_keys(value.keys()):
                        if mol.strip('mol') in moleculesToKeep:
                            new_mol_number = moleculesToKeep.index(mol.strip('mol')) + 1
                            newInputData[namelist][variable]['mol%i'%new_mol_number] = value[mol]
                    if variable == 'pmswmt':
                        idata = inputData[namelist][variable].copy()
                        probs = list(float(idata[mol].rstrip('d0')) for mol in sort_keys(idata.keys()))
                        if probs[-1] != 1.:
                            print('problem with swap probabilities from input file')
                            i = 0
                            for mol in sort_keys(idata.keys()):
                                i += 1
                                newInputData[namelist][variable][mol] = '%7.4fd0'%(i/len(idata.keys())) 
                elif len(value.keys()) == nbox_old:
                    newInputData[namelist][variable] = {}
                    # initialize with no boxes and only put in those needed
                    for box in sort_keys(value.keys()):
                        if box.strip('box') in boxesToKeep:
                            new_box_number = boxesToKeep.index(box.strip('box')) + 1
                            newInputData[namelist][variable]['box%i'%new_box_number] = value[box]
    newInputData['SIMULATION_BOX'] = {}
    for box in sort_keys(inputData['SIMULATION_BOX'].keys()):
        if box.strip('box') in boxesToKeep:
            new_box = 'box%i'%(boxesToKeep.index(box.strip('box')) + 1)
            newInputData['SIMULATION_BOX'][new_box] = {}
            for keyToCopy in ['nghost','defaults','temperature','pressure','initialization data']:
                newInputData['SIMULATION_BOX'][new_box][keyToCopy] = inputData['SIMULATION_BOX'][box][keyToCopy]
            newInputData['SIMULATION_BOX'][new_box]['dimensions'] = restartData['box dimensions'][box].rstrip('\n')
            boxlx = float(restartData['box dimensions'][box].split()[0])
            newInputData['SIMULATION_BOX'][new_box]['rcut'] = inputData['SIMULATION_BOX'][box]['rcut']
            for mol in ['mol%i'%i for i in range(1,len(moleculesToKeep)+1)]:
                newInputData['SIMULATION_BOX'][new_box][mol] = '%i'%nmolty_by_box[box][mol]
    for section in [i for i in ['MOLECULE_TYPE','MC_SWAP','INTRAMOLECULAR_OH15','UNIFORM_BIASING_POTENTIALS','INTERMOLECULAR_EXCLUSION']
                        if i in inputData.keys()]:
        newInputData[section] = {}
        for mol in sort_keys(inputData[section].keys()):
            if mol.strip('mol') in moleculesToKeep:
                new_mol_number = moleculesToKeep.index(mol.strip('mol')) + 1
                if section == 'MC_SWAP':
                    newInputData[section]['mol%i'%new_mol_number] = {}
                    box_pairs = []
                    for box_pair in inputData[section][mol]['box1 box2']:
                        # if keeping both boxes in box pair, keep pair in swap info
                        if ('%i'%box_pair[0] in boxesToKeep) and ('%i'%box_pair[1] in boxesToKeep):
                            box_pairs.append(box_pair)
                    if len(box_pairs) == 0:
                        # give it a box pair, make pair = a box molecules are currently in and box 1
                        box_pairs = []
                        for box in nmolty_by_box.keys():
                            if nmolty_by_box[box]['mol%i'%new_mol_number] > 0:
                                box_pairs.append( [1,int(box.strip('box'))] )
                    newInputData[section]['mol%i'%new_mol_number]['nswapb'] = len(box_pairs)
                    newInputData[section]['mol%i'%new_mol_number]['box1 box2'] = box_pairs
                    if len(inputData[section][mol]['pmswapb']) != len(box_pairs):
                        newInputData[section]['mol%i'%new_mol_number]['pmswapb'] = [(i+1)/len(box_pairs) for i in range(len(box_pairs))]
                    else:
                        newInputData[section]['mol%i'%new_mol_number]['pmswapb'] = inputData[section][mol]['pmswapb']
                elif section == 'UNIFORM_BIASING_POTENTIALS':
                    newInputData[section]['mol%i'%new_mol_number] = {}
                    for box in inputData[section][mol].keys():
                        if box.strip('box') in boxesToKeep:
                            new_box_number = boxesToKeep.index(box.strip('box')) + 1
                            newInputData[section]['mol%i'%new_mol_number]['box%i'%new_box_number] = inputData[section][mol][box]
                elif section == 'INTERMOLECULAR_EXCLUSION':
                    if 'mol%i'%new_mol_number in inputData[section][mol].keys():
                        mol1 = 'mol%i'%new_mol_number
                        newInputData[section][mol1] = {}
                        for unit1, mols in inputData[section][mol].items():
                            for mol2 in [i for i in mols.keys() if i.strip('mol') in moleculesToKeep]:
                                if mol1 not in input_data[section].keys():
                                    input_data[section][mol1] = {}
                                if unit1 not in input_data[section][mol1].keys():
                                    input_data[section][mol1][unit1] = {}
                                input_data[section][mol1][unit1][mol2] = mols[mol2]
                else:
                    newInputData[section]['mol%i'%new_mol_number] = inputData[section][mol]
    return newRestartData, newInputData

from MCFlow.file_formatting.writer import sort_keys
from MCFlow.file_formatting import reader, writer

if __name__ == '__main__':
    from MCFlow.parser import ChangeInput
    my_parser = ChangeInput()
    my_parser.molecules()
    my_parser.parser.add_argument('-b','--boxes',help='boxes to keep',type=str,nargs='+')
    args = vars(my_parser.parse_args())

    for feed in args['feeds']:
        for sim in args['indep']:
            path = feed + '/' + str(sim) + '/'
            input_data = reader.read_fort4(path + args['input'])
            nmolty, nbox = int(input_data['&mc_shared']['nmolty']), int(input_data['&mc_shared']['nbox'])
            restart_data = reader.read_restart(path + args['restart'], nmolty, nbox)
            newRes, newIn = removeExtraInfo(args['boxes'], restart_data, input_data)
            writer.write_fort4( newIn, path + 'fort.4.purified')
            writer.write_restart( newRes, path + 'fort.77.purified')
