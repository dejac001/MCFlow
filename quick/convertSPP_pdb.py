'''
Take a pdb file that's been replicated with silanol groups.
Take out silanol groups for tabulation.
'''

import argparse, copy
import numpy as np
from file_formatting import reader, writer
from file_formatting.addMolecules import calculate_distance

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Separate SPP pdb file into portions to be'
                                                   ' tabulated and to be provided explicitly')
    parser.add_argument('-f','--file',help='input file', type=str,
                        default='/Users/dejacor/Documents/'
                                'Lulea-University-Butanol/'
                                'structures/SPP_film_8by4_4sheet.pdb')
    args = vars(parser.parse_args())

    # read input pdb
    data = reader.PDB(args['file'])
    main_dir = args['file'][:args['file'].rfind('/')] + '/'

    # silanol groups
    silanol_groups = []
    all_data = copy.deepcopy(data['coords'])
    abc = data['box info']['a'], data['box info']['b'], data['box info']['c']
    all_sioh_data = {'279':[],'180':[],'181':[],'182':[],'186':[]}
    fort4_topmon = {'279':'Oh','180':'Se','181':'Oe','186':'Ot','182':'H'}
    charges_fort4 = {'279': -0.375,'180':1.429,'181':-0.739, '182':0.435,'186':-0.45}
    charges_topmon = {'Si':1.5,'O':-0.75,'Oh':-0.375, 'Se':0.,'Oe':0.,'Ot':-0.3}
    for Hxyz in all_data['H']:
        my_group_data = {'182':[Hxyz]}
        data['coords']['H'].remove(Hxyz)
        for Oxyz in all_data['O']:
            # calculate distance
            distance = calculate_distance(Hxyz, Oxyz, abc)
            if distance < 2.:
                # this is a bond
                my_group_data['181'] = [Oxyz]
                data['coords']['O'].remove(Oxyz)
                for Sixyz in all_data['Si']:
                    distance = calculate_distance(Sixyz, Oxyz, abc)
                    if distance < 2.:
                        # this is a bond
                        my_group_data['180'] = [Sixyz]
                        data['coords']['Si'].remove(Sixyz)

                        my_key = '279'
                        # 186 is oxygen on dimer
                        # 279 is oxygen on monomer
                        for Oxyz2 in [ i for i in all_data['O'] if i != Oxyz]:
                            distance = calculate_distance(Sixyz, Oxyz2, abc)
                            if (distance < 2):
                                # this is a bond
                                if my_key not in my_group_data.keys():
                                    my_group_data[my_key] = []
                                if (Oxyz2 in data['coords']['O']):
                                    # this oxygen has not been added to silanols yet
                                    my_group_data[my_key].append( Oxyz2 )
                                    data['coords']['O'].remove(Oxyz2)
                                else:
                                    # oxygen has been marked in a previous silanol group
                                    # O2SiOH-O-SiOHO2
                                    # need to find a way to make dimer here
                                    my_key = '186'
                                    if '279' in my_group_data.keys():
                                        my_group_data['186'] = my_group_data.pop('279')
                                    else:
                                        my_group_data['186'] = []
                                    my_group_data['186'].append( Oxyz2 )
        silanol_groups.append(my_group_data.copy())

    # need to find a way to make dimer here
    initial_monomer_list = [silanol_groups[k] for k in range(len(silanol_groups)) if '279' in silanol_groups[k].keys()]
    monomer_list = [silanol_groups[k] for k in range(len(silanol_groups)) if '279' in silanol_groups[k].keys()]
    final_sioh_groups = []
    for i in [k for k in range(len(silanol_groups)) if '186' in silanol_groups[k].keys()]: # dimers
        for monomer in initial_monomer_list:
            # determine if they share an oxygen
            for iO, O1 in enumerate(silanol_groups[i]['186']):
                if O1 in monomer['279']:
                    dimer =  copy.deepcopy(silanol_groups[i])
                    dimer['186'].append( dimer['186'].pop(iO) ) # put shared silanol at end
                    for atom in monomer.keys():
                        if atom == '279':
                            atom_num = '186'
                        else:
                            atom_num = atom
                        for xyz in monomer[atom]:
                            for dxyz in dimer[atom_num]:
                                if np.allclose(xyz,dxyz):
                                    break
                            else:
                                # loop fell w/o finding any that are the same
                                dimer[atom_num].append( xyz )
                    final_sioh_groups.append(dimer)
                    monomer_list.remove( monomer )
                    break
    final_sioh_groups += monomer_list
    atom_order = ['182','181','180','279','186']
    for i, group_data in enumerate(final_sioh_groups):
        my_charge = 0
        my_data = {'atoms':[],'coords':[]}
        for atom_num in [i for i in  atom_order if i in group_data]:
            for xyz in group_data[atom_num]:
                my_charge += charges_fort4[atom_num]
                all_sioh_data[atom_num].append(xyz)
                my_data['coords'].append( xyz )
                my_data['atoms'].append( fort4_topmon[atom_num] )
                if atom_num in fort4_topmon.keys():
                    # put back in tabulation if needed
                    atom = fort4_topmon[atom_num]
                    if atom != 'H':
                        if atom not in data['coords'].keys():
                            data['coords'][atom] = []
                        data['coords'][atom].append(xyz)
        if (abs(my_charge) > 1e-05):
            print('charge of silanol group is %8.4f'%my_charge)
            print(group_data)
        writer.xyz(main_dir + 'silanol%i.xyz'%i, my_data)
    pdb_data = {'atoms':[],'coords':[]}
    for atom in data['coords']:
        for xyz in data['coords'][atom]:
            pdb_data['atoms'].append( atom )
            pdb_data['coords'].append( xyz )
    writer.xyz( args['file'].rstrip('.pdb') + '_to_tabulate.xyz', pdb_data)
    
    total_charge_fort4 = 0
    for atom in charges_fort4.keys():
        total_charge_fort4 += charges_fort4[atom]*len(all_sioh_data[atom])
    total_charge_tab = 0
    for atom in data['coords'].keys():
        if len(data['coords'][atom]) > 0:
            total_charge_tab += charges_topmon[atom]*len(data['coords'][atom])
    print('total charge fort4: ',total_charge_fort4)
    print('total charge tabulation: ',total_charge_tab)
    total_O = len(all_sioh_data['279']) + len(all_sioh_data['186']) + len(all_sioh_data['181']) + len(data['coords']['O'])
    total_Si = len(all_sioh_data['180']) + len(data['coords']['Si'])
    total_H = len(all_sioh_data['182'])
    if total_H%2 != 0:
        print('error in data manipulation')
    else:
        nWater = total_H/2
        nSiO2 = (total_O - nWater)/2
        print('unit cell is %i SiO2 + %i H2O'%(nSiO2, nWater))
        print('unit cell vectors were ',abc)
        print('-----total H:',nWater*2, total_H)
        print('-----total Si:',nSiO2, total_Si)
        print('-----total O:',nWater + nSiO2*2, total_O)
