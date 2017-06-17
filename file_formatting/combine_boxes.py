def combine(inRestartData,inInputData):
    '''
    combine restart and input files. no shared boxes currently
    :param listRestartData: restart Data for all simulations.
                            boxes will be numbered in this order.
    :param listInputData: input data (fort.4) for all simulations
    :return:
    '''
    listRestartData = inRestartData[:]
    listInputData = inInputData[:]
    if len(listRestartData) != len(listInputData):
        print('Amount of restart data does not equal amount of input data')
        quit()
    newRestartData = {'max displacement' : {},
                        'number of cycles' : listRestartData[0]['number of cycles'],
                      'box dimensions':{},
                      'nunit' : {},
                       'mol types':[],
                      'box types':[],
                      'coords':[]}
    nbox_new = 0
    nmolty_new = 0
    nchain_new = 0
    nbox_disp = 0
    moleculeTypesToKeep = []
    num_molty_by_box = {}
    nbox_old = 0
    for iRestart, fileRestart in enumerate(listRestartData):
        ibox = nbox_new
        if 'translation' in newRestartData['max displacement'].keys():
            nbox_old = len(newRestartData['max displacement']['translation'].keys())
        for disp_type, values in sorted(fileRestart['max displacement'].items()):
            nbox_disp = nbox_new
            if disp_type == 'atom translation':
                if disp_type not in newRestartData['max displacement'].keys():
                    newRestartData['max displacement'][disp_type] = values
            else:
                if disp_type not in newRestartData['max displacement'].keys():
                    newRestartData['max displacement'][disp_type] = {}
                for box in sort_keys(values.keys()):
                    nbox_disp += 1
                    new_box = 'box%i'%nbox_disp
                    if type(values[box]) != type({}):
                        newRestartData['max displacement'][disp_type][new_box] = values[box]
                    else:
                        newRestartData['max displacement'][disp_type][new_box] = {}
                        for mol in sort_keys(values[box].keys()):
                            mol_input = listInputData[iRestart]['MOLECULE_TYPE'][mol]
                            if mol_input not in moleculeTypesToKeep:
                                moleculeTypesToKeep.append(mol_input)
                            mol_num = moleculeTypesToKeep.index(mol_input) + 1
                            new_mol = 'mol%i'%mol_num
                            newRestartData['max displacement'][disp_type][new_box][new_mol] = values[box][mol]
        nbox_new = nbox_disp
        nmolty_new = len(moleculeTypesToKeep)
        ibox_current = len(newRestartData['box dimensions'].keys())
        ibox = ibox_current
        for box in sort_keys(fileRestart['box dimensions'].keys()):
            ibox += 1
            newRestartData['box dimensions']['box%i'%ibox] = fileRestart['box dimensions'][box]
        if 'nchain' not in newRestartData.keys():
            newRestartData['nchain'] = fileRestart['nchain']
            nchain_new = int(fileRestart['nchain'].split()[0])
        else:
            ichain = int(fileRestart['nchain'].split()[0])
            newRestartData['nchain'] = newRestartData['nchain'].replace('%i'%nchain_new,'%i'%(nchain_new + ichain))
            nchain_new += ichain
        newRestartData['nmolty'] = fileRestart['nmolty'].replace(fileRestart['nmolty'].split()[0],
                                                                 '%i'%len(moleculeTypesToKeep))
        for mol in sort_keys(fileRestart['nunit'].keys()):
            mol_input = listInputData[iRestart]['MOLECULE_TYPE'][mol]
            new_mol = 'mol%i'%(moleculeTypesToKeep.index(mol_input) + 1)
            newRestartData['nunit'][new_mol] = fileRestart['nunit'][mol]
        num_box_prev = nbox_old
        for i, boxNum in enumerate(fileRestart['box types']):
            new_box = '%i'%(int(boxNum) + num_box_prev)
            if ('box' + new_box) not in num_molty_by_box:
                num_molty_by_box['box' + new_box] = {}
            mol_new = '%i'%(moleculeTypesToKeep.index(listInputData[iRestart]['MOLECULE_TYPE']['mol' +
                                                fileRestart['mol types'][i]]) + 1)
            if ('mol' + mol_new) not in num_molty_by_box['box' + new_box]:
                num_molty_by_box['box' + new_box]['mol' + mol_new] = 1
            else:
                num_molty_by_box['box' + new_box]['mol' + mol_new] += 1
            newRestartData['mol types'].append(mol_new)
            newRestartData['box types'].append(new_box)
        newRestartData['coords'] += fileRestart['coords']
    # add in max displacements for all molecules in boxes they weren't in before
    for box, value in newRestartData['max displacement']['translation'].items():
        for mol in ['mol%i'%i for i in range(1,len(moleculeTypesToKeep)+1)]:
            if mol not in value.keys():
                rmtra = listInputData[0]['&mc_simple']['rmtra']
                newRestartData['max displacement']['translation'][box][mol] = ' '.join([rmtra, rmtra, rmtra]) + '\n'
    for box, value in newRestartData['max displacement']['rotation'].items():
        for mol in ['mol%i'%i for i in range(1,len(moleculeTypesToKeep)+1)]:
            if mol not in value.keys():
                rmrot = listInputData[0]['&mc_simple']['rmrot']
                newRestartData['max displacement']['rotation'][box][mol] = ' '.join([rmrot, rmrot, rmrot]) + '\n'
    for box, value in newRestartData['max displacement']['fluctuating charge'].items():
        for mol in ['mol%i'%i for i in range(1,len(moleculeTypesToKeep)+1)]:
            if mol not in value.keys():
                rmfluc = '0.100000000000000'
                newRestartData['max displacement']['fluctuating charge'][box][mol] = rmfluc
    for box in ['box%i'%i for i in range(1,nbox_new+1)]:
        if box not in num_molty_by_box.keys():
            num_molty_by_box[box] = {}
        for mol in ['mol%i'%i for i in range(1,len(moleculeTypesToKeep)+1)]:
            if mol not in num_molty_by_box[box].keys():
                num_molty_by_box[box][mol] = 0
    # combine input data
    newInputData = {}
    nbox_combined = 0
    for iInput, fileInput in enumerate(listInputData):
        namelists = [i for i in fileInput.keys() if i.startswith('&')]
        sections = [i for i in fileInput.keys() if i.isupper()]
        # namelist stuffs--------------------------------------
        if iInput == 0:
            valuesToUpdate = {}
            for namelist in namelists:
                if namelist not in newInputData.keys():
                    newInputData[namelist] = {}
                for variable, values in fileInput[namelist].items():
                    if variable == 'pmsatc':
                        newInputData[namelist][variable] = values
                    elif (variable == 'nchain') or (variable == 'nbox') or (variable == 'nmolty'):
                        old_value = fileInput[namelist][variable].split()[0]
                        if (variable == 'nchain'):
                            updated_value = '%i'%nchain_new
                        elif (variable == 'nbox'):
                            updated_value = '%i'%nbox_new
                        elif (variable == 'nmolty'):
                            updated_value = '%i'%nmolty_new
                        newInputData[namelist][variable] = fileInput[namelist][variable].replace(old_value,updated_value)
                    else:
                        # not a special case
                        if type(values) == type(''):
                            newInputData[namelist][variable] = values
                        elif 'mol1' in values.keys():
                            newInputData[namelist][variable] = {}
                            if namelist not in valuesToUpdate.keys():
                                valuesToUpdate[namelist] = [variable]
                            else:
                                valuesToUpdate[namelist].append(variable)
                            for mol in sort_keys(values.keys()):
                                mol_input = fileInput['MOLECULE_TYPE'][mol]
                                if mol_input in moleculeTypesToKeep:
                                    mol_new = moleculeTypesToKeep.index(mol_input) + 1
                                    newInputData[namelist][variable]['mol%i'%mol_new] = values[mol]
                        elif 'box1' in values.keys():
                            newInputData[namelist][variable] = {}
                            if namelist not in valuesToUpdate.keys():
                                valuesToUpdate[namelist] = [variable]
                                nbox_current = 0
                            else:
                                valuesToUpdate[namelist].append(variable)
                                nbox_current = len(newInputData[namelist][variable].keys())
                            ibox_current = 0
                            for box in sort_keys(values.keys()):
                                ibox_current += 1
                                new_box = 'box%i'%(nbox_current+ibox_current)
                                newInputData[namelist][variable][new_box] = values[box]
        else:
            for namelist, variables in valuesToUpdate.items():
                for var in variables:
                    values = fileInput[namelist][var].copy()
                    if 'mol1' in values.keys():
                        for mol in sort_keys(values.keys()):
                            mol_input = fileInput['MOLECULE_TYPE'][mol]
                            if (mol_input in moleculeTypesToKeep):
                                mol_new = 'mol%i'%(moleculeTypesToKeep.index(mol_input) + 1)
                                if (mol_new not in newInputData[namelist][var].keys()):
                                    newInputData[namelist][var][mol_new] = values[mol]
                    elif 'box1' in values.keys():
                        nbox_current = len(newInputData[namelist][var].keys())
                        ibox_current = 0
                        for box in sort_keys(values.keys()):
                            ibox_current += 1
                            new_box = 'box%i'%(nbox_current+ibox_current)
                            newInputData[namelist][var][new_box] = values[box]


        # section stuffs-------------------------------------------------------------
        for section in sections:
            if ((iInput == 0) and (section not in ['SIMULATION_BOX','MOLECULE_TYPE','MC_SWAP','INTRAMOLECULAR_OH15',
                                                  'UNIFORM_BIASING_POTENTIALS'])):
                # copy the section
                newInputData[section] = fileInput[section].copy()
            else:
                if section not in newInputData.keys():
                    newInputData[section] = {}
                if section == 'SIMULATION_BOX':
                    nbox_current = len(newInputData[section].keys())
                    ibox_current = 0
                    for box in sort_keys(fileInput[section].keys()):
                        ibox_current +=1
                        new_box = 'box%i'%(nbox_current + ibox_current)
                        newInputData[section][new_box] = fileInput[section][box].copy()
                elif section in ['MOLECULE_TYPE','INTRAMOLECULAR_OH15','MC_SWAP']:
                    for mol in sort_keys(fileInput[section].keys()):
                        mol_input = fileInput['MOLECULE_TYPE'][mol]
                        if (mol_input in moleculeTypesToKeep):
                            mol_new = 'mol%i'%(moleculeTypesToKeep.index(mol_input) + 1)
                            if (mol_new not in newInputData[section].keys()):
                                newInputData[section][mol_new] = fileInput[section][mol]
                                if section == 'MC_SWAP':
                                    for i, box_pair in enumerate(newInputData[section][mol_new]['box1 box2']):
                                        if (box_pair[0] == box_pair[1]) and (box_pair[0] != 1):
                                            # find a box that it's in and put it there
                                            j = newRestartData['mol types'].index(mol_new.strip('mol'))
                                            my_new_box = int(newRestartData['box types'][j])
                                            newInputData[section][mol_new]['box1 box2'][i] = [1,my_new_box]
                                        elif max(box_pair) <= nbox_combined:
                                            my_pair = box_pair[:]
                                            my_pair[my_pair.index(max(my_pair))] = nbox_combined + max(my_pair)
                                            newInputData[section][mol_new]['box1 box2'][i] = my_pair
                elif section == 'UNIFORM_BIASING_POTENTIALS':
                    for mol in sort_keys(fileInput[section].keys()):
                        mol_input = fileInput['MOLECULE_TYPE'][mol]
                        if (mol_input in moleculeTypesToKeep):
                            mol_new = 'mol%i'%(moleculeTypesToKeep.index(mol_input) + 1)
                            newInputData[section][mol_new] = {}
                            if (mol_new not in newInputData[section].keys()):
                                # copy all
                                for box, value in enumerate(fileInput[section][mol_new].keys()):
                                    newInputData[section][mol_new][box] = value.copy()
                            else:
                                # add boxes
                                nbox_current = len(newInputData[section][mol_new].keys())
                                ibox_current = 0
                                for box in sort_keys(fileInput[section][mol].keys()):
                                    ibox_current += 1
                                    new_box = 'box%i'%(nbox_current+ibox_current)
                                    newInputData[section][mol_new][new_box] = fileInput[section][mol][box]
        nbox_combined += int(fileInput['&mc_shared']['nbox'].split()[0])
    # add all molecules to correct boxes
    for box in newInputData['SIMULATION_BOX'].keys():
        for mol in num_molty_by_box[box]:
            newInputData['SIMULATION_BOX'][box][mol] = '%i'%num_molty_by_box[box][mol]
            if box not in newInputData['UNIFORM_BIASING_POTENTIALS'][mol].keys():
                newInputData['UNIFORM_BIASING_POTENTIALS'][mol][box] = '0.0d0'
    # correct pmswmt
    # if doing swap moves before, continue doing them. Split up swap moves evenly
    probabilities_swap = [float(newInputData['&mc_swap']['pmswmt'][i].rstrip('d0'))
                          for i in sort_keys(newInputData['&mc_swap']['pmswmt'].keys())]
    nswapty = len([i for i in probabilities_swap if i > 0.])
    total_prob = 0.
    for i, prob in enumerate(probabilities_swap):
        if prob > 0.:
            total_prob += 1/nswapty
            newInputData['&mc_swap']['pmswmt']['mol%i'%(i+1)] = '%7.4fd0'%total_prob
    # if doing volume moves before, continue doing them. Split up moves evenly
    probabilities_volume = [float(newInputData['&mc_volume']['pmvlmt'][i].rstrip('d0'))
                            for i in sort_keys(newInputData['&mc_volume']['pmvlmt'].keys())]
    nvolty = len([i for i in probabilities_volume if i > 0.])
    total_prob = 0.
    for i, prob in enumerate(probabilities_volume):
        if prob > 0.:
            total_prob += 1/nvolty
            newInputData['&mc_volume']['pmvlmt']['box%i'%(i+1)] = '%7.4fd0'%total_prob
    n_total = int(newInputData['&mc_shared']['nchain'])
    tavol = float(newInputData['&mc_volume']['tavol'])
    newInputData['&mc_volume']['pmvol'] = '%7.4fd0'%(1/(tavol*n_total))
    # correct translation, rotation, and cbmc stuffs
    nMol = []
    nBeads = []
    for molNum in ['%i'%i for i in range(1,len(moleculeTypesToKeep)+1)]:
        nMol.append( newRestartData['mol types'].count(molNum))
        if 'T T' in newInputData['MOLECULE_TYPE']['mol%s'%molNum]:
            nBeads.append( 0 )
        else:
            nBeads.append( int(newRestartData['nunit']['mol%s'%molNum]))
    cbmc_fraction, translation_fraction = calculateProbs(nMol, nBeads)
    for moveType in ['pmtrmt','pmromt']:
        nMol_trans = 0
        i = -1
        for mol in sort_keys(newInputData['&mc_simple'][moveType].keys()):
            i += 1
            if not newInputData['&mc_simple'][moveType][mol] == '0.0d0':
                nMol_trans += nMol[i]
        total_prob = 0.
        i = -1
        for mol in sort_keys(newInputData['&mc_simple'][moveType].keys()):
            i+=1
            if not newInputData['&mc_simple'][moveType][mol] == '0.0d0':
                total_prob += nMol[i]/nMol_trans
                newInputData['&mc_simple'][moveType][mol] = '%7.3fd0'%total_prob
    newInputData['&mc_swap']['pmswap'] = '0.600d0'
    newInputData['&mc_cbmc']['pmcb'] = '%7.4fd0'%(0.6 + 0.4*cbmc_fraction)
    newInputData['&mc_simple']['pmtra'] = '%7.4fd0'%(0.6 + 0.4*translation_fraction)
    cbmcProbs = calcCBMCmolTy(nMol, nBeads)
    i = -1
    for mol in sort_keys(newInputData['&mc_cbmc']['pmcbmt'].keys()):
        i+=1
        newInputData['&mc_cbmc']['pmcbmt'][mol] = '%7.4fd0'%cbmcProbs[i]
    return newRestartData, newInputData

def main(inputPaths, boxesToKeep):
    restarts = []
    inputs = []
    for my_path, my_boxes in zip(inputPaths, boxesToKeep):
        input_data = read_fort4(my_path + '/fort.4')
        nmolty, nbox = input_data['&mc_shared']['nmolty'], input_data['&mc_shared']['nbox']
        restart_data = read_restart(my_path + '/fort.77', nmolty, nbox)
        restart_new, fort4_new = purify.removeExtraInfo(my_boxes,restart_data,input_data)
        restarts.append(restart_new)
        inputs.append(fort4_new)
    restart_combined, fort4_combined = combine(restarts, inputs)
    return restart_combined, fort4_combined

def testing():
    path1 = '/Volumes/dejac002/Seagate Dashboard//BETO-simulations/C6-screening/setups/mol12-6.2e-08bar/1'
    path2 = '/Volumes/dejac002/Seagate Dashboard//BETO-simulations/C6-screening/setups/mol1-1.8e-01bar/1'
    paths = [path1, path2]
    boxesToKeep = [['1','2'],['2']]
    restart_all, fort4_all = main(paths, boxesToKeep, path1)
    write_restart(restart_combined, path1 + '/fort.77.new')
    write_fort4(fort4_combined, path1 + '/fort.4.new')

#import argparse
from file_formatting.reader import read_restart, read_fort4
from file_formatting.writer import write_restart, write_fort4, sort_keys
from specialized.probchanger import calculateProbs, calcCBMCmolTy
from file_formatting import purify
#
# parser = argparse.ArgumentParser(description = 'Add boxes together for single simulation')
# parser.add_argument('-p','--paths',help='paths to different input and restart files', type=str,nargs='+')

# args = vars(parser.parse_args())
# if len(args['paths']) < 2:
#     print('you cant combine only one file')
#     quit()
#
# print(read_fort4(args['paths'][0] + 'fort.4'), read_restart(args['paths'][0] + '/fort.77'))
