def add2fort77(old_file, coords, box):
    '''
    coords: {'mol#':[],...
            where mol# is an integer. The indexes of list
            contain printable coordinates of all bead of this molecule in 
            order of fort.4 file.
    '''
    f = open('fort.77.new','w')
    molec_types = []
    box_types = []
    FinishDispl = False
    molecType = False
    boxType = False
    Start = True
    for line in open(old_file):
        if Start and (len(line.split()) == 1) and (int(line.split()[0]) < 10000):
            nchain_old = int(line.split()[0])
            FinishDispl = True
            Start = False
            f.write(line.replace('%i'%nchain_old, '%i'%(sum([nchain_old] + 
                                            [len(coords[i]) for i in coords.keys()]))))
        elif FinishDispl:
            if line.split()[0] == '1': 
                molecType = True
                FinishDispl = False
                molec_types += line.split()
            else:
                f.write(line)
        elif molecType:
            molec_types += line.split()
            if len(molec_types) == nchain_old:
                # write everything at once
                my_line = ''
                for molec in molec_types:
                    my_line += '    %s'%molec
                for molNum in sorted(coords.keys()):
                    for i in range(len(coords[molNum])):
                        my_line += '    %s'%molNum
                f.write(my_line + '\n')
                boxType = True
                molecType = False
        elif boxType:
            box_types += line.split()
            if len(box_types) == nchain_old:
                my_line = ''
                for ibox in box_types:
                    my_line += '    %s'%ibox
                for molNum in sorted(coords.keys()):
                    for i in range(len(coords[molNum])):
                        my_line += '    %i'%box
                f.write(my_line + '\n')
                boxType = False
        else:
            f.write(line)
    for molNum in sorted(coords.keys()):
        for myPos in coords[molNum]:
            f.write(myPos)
    f.close()


import os

if __name__ == '__main__':
    molecules = []
    path_to_struc = '../../../structures/'
    for file in [i for i in os.listdir(path_to_struc) if 'fort77_mol' in i]:
        mols = {}
        print(file)
        for line in open(path_to_struc + file):
            if len(line.split()) == 4:
                bead, x, y, z = line.split()
                if bead in mols.keys(): bead = bead + str(len([i for i in mols.keys() if bead in i]))
                mols[bead] = '%s    %s      %s\n'%(x,y,z)
        molecules.append(mols)
    str_coordinates = {'2':[],'3':[]} #mols 2 and 3
    for mol in molecules:
        my_mol_str = ''
        for bead in [i for i in mols.keys() if 'Cr' in i]:
           my_mol_str += mol[bead] + '-0.3750000\n'
        if 'C' in mol.keys(): # stands for Si
            my_mol_str += mol['C'] + '1.4290000\n'
            qO, qH = -0.739, 0.435
            nmoltype = '2'
        elif 'Sn' in mol.keys():
            my_mol_str += mol['Sn'] + '1.5550000\n'
            qO, qH = -0.887, 0.457
            nmoltype = '3'
        my_mol_str += mol['Os'] + '%e\n'%qO
        my_mol_str += mol['H'] + '%e\n'%qH
        str_coordinates[nmoltype].append(my_mol_str)
    add2fort77('fort.77',str_coordinates, 1)
