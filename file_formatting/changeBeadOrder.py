def changeBeadOrder(mols, newOrders,  input_data, restart_data):
    '''
    new order is lists for each mol as string
    Ex: newOrders[0] might be  [4,6,1,2,3,5]
            for mol# = mols[0]
            changes 4 --> 1 (old bead 4 to new bead 1)
            changes 6 --> 2
            changes 1 --> 3
            etc...
    '''
    new_input = input_data.copy()
    new_rest = restart_data.copy()
    new_input['MOLECULE_TYPE'] = input_data['MOLECULE_TYPE'].copy()
    for imol, mol in enumerate(mols):
        new_input['MOLECULE_TYPE']['mol%s'%mol] = ''
        mol_lines = input_data['MOLECULE_TYPE']['mol%s'%mol].split('\n')
        mol_units = [0 for i in range(int(restart_data['nunit']['mol%s'%mol]))]
        iunit = -1
        for iline, my_line_str in enumerate(mol_lines):
            if '!' in my_line_str:
                line_str = my_line_str[:my_line_str.find('!')] # strip comments 
            else:
                line_str = my_line_str
            if 'growpoint' in mol_lines[iline-1]:
                new_input['MOLECULE_TYPE']['mol%s'%mol] += my_line_str + '\n'
            elif (line_str.replace(' ','').isdigit()):
                # we need to change bead info
                if (('unit' in mol_lines[iline-1]) and 
                        ('ntype' in mol_lines[iline-1])
                         and ('leaderq' in mol_lines[iline-1])):
                    ibead_old = int(line_str.split()[0])
                    iunit = newOrders[imol].index(ibead_old)
                    mol_units[iunit] = ''
                    new_bead = '%i'%(newOrders[imol].index(ibead_old) + 1)
                    new_line = [new_bead, line_str.split()[1], new_bead]
                else:
                    new_line = []
                    for bead in line_str.split()[:-1]:
                        ibead_old = int(bead)
                        new_line.append('%i'%(newOrders[imol].index(ibead_old) + 1))
                    new_line.append(line_str.split()[-1])
                mol_units[iunit] += ' '.join(new_line) + '\n'
            elif iunit == -1:
                new_input['MOLECULE_TYPE']['mol%s'%mol] += my_line_str + '\n'
            else:
                mol_units[iunit] += my_line_str + '\n'
        for i in mol_units:
            new_input['MOLECULE_TYPE']['mol%s'%mol] += i
        if 'first_bead_to_swap' in new_input['&mc_cbmc'].keys():
            try:
                new_input['&mc_cbmc']['first_bead_to_swap']['mol%s'%mol] = '1'
            except TypeError:
                new_input['&mc_cbmc']['first_bead_to_swap'] = {'mol1':'1'}
    for jmol, mol in enumerate(restart_data['mol types']):
        if mol in mols:
            molNum = mols.index(mol)
            # we need to change coordinate orders
            mol_coords = restart_data['coords'][jmol][:]
            new_rest['coords'][jmol] = [0 for i in range(len(mol_coords))]
            for ibead, bead in enumerate(mol_coords):
                new_rest['coords'][jmol][newOrders[molNum].index(ibead+1)] = bead.copy()
            if 0 in new_rest['coords'][jmol]:
                print('error making new coords')
                quit()
    return new_input, new_rest
                

from MCFlow.file_formatting.reader import read_fort4, read_restart
from MCFlow.file_formatting.writer import write_fort4, write_restart

if __name__ == '__main__':
    import argparse, os

    parser = argparse.ArgumentParser(description = 'Analyze and collect data for vapor adsorption')
    parser.add_argument('-f','--in_file',help ='fort.4 file to change',type=str,
                            default=os.getcwd() + '/fort.4')
    parser.add_argument('-r','--in_rest',help ='fort.77 file to change',type=str,
                            default=os.getcwd() + '/fort.77')
    parser.add_argument('-b','--new_beads',help ='''new beads for each molecule type:
    New order is lists for each mol.    Ex: you might pass in `4,6,1,2,3,5' to get -->  [4,6,1,2,3,5]. This changes
    4 --> 1 (old bead 4 to new bead 1), 6 --> 2, 1 --> 3, etc...N.B: mols are read in string as list''',type=str,nargs='+')
    parser.add_argument('-m','--mols',help='mols to change coresponding to those passed into new beads', type=str,nargs='+')
    args = vars(parser.parse_args())
    if len(args['new_beads']) != len(args['mols']):
        print('beads and mols not specified correctly')
        for key in ['new_beads','mols']:
            print('{} {} len: {}'.format(key, args[key], len(args[key])))
        quit()
    new_beads = [list(map(int,i.split(','))) for i in args['new_beads']]
    input_data = read_fort4(args['in_file'])
    nmolty, nbox = int(input_data['&mc_shared']['nmolty']), int(input_data['&mc_shared']['nbox'])
    restart_data = read_restart(args['in_rest'],nmolty, nbox)
    new_in, new_res = changeBeadOrder(args['mols'], new_beads, input_data.copy(), restart_data.copy())
    if args['in_file'].rfind('/') != -1:
        new_path = args['in_file'][:args['in_file'].rfind('/')] + '/'
    else:
        new_path = ''
    write_fort4(new_in, new_path + 'fort.4.new')
    write_restart(new_res, new_path + 'fort.77.new')
