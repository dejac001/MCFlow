def keep_zeo(path, fort4, restart_file, temp, boxlx_vapor):
    '''
    Takes fort4 and restart file from previous simulation. Only
    keeps structure in zeolite box for a subsequent simulation.
    N.B.:   give full path to all file names
    '''
    def writeNewBoxInfo(mlclNum, boxNum):
        if len(mlclNum) != len(boxNum):
            print('number of molecules not consistent')
            print(len(mlclNum), len(boxNum))
            quit()
        nchain_new = 0
        for i in range(len(mlclNum)):
            if boxNum[i] == '1':
                nchain_new += 1
                f77new.write('{} '.format(mlclNum[i]))
        f77new.write('\n')
        for i in range(len(mlclNum)):
            if boxNum[i] == '1':
                f77new.write('{} '.format(boxNum[i]))
        f77new.write('\n')
        return nchain_new
    def findCoordToKeep(mlclNum, boxNum, beads, nmolty_keep):
        total_coordinates = 0
        coordinates_to_keep = []
        for i in range(len(mlclNum)):
            moltype = int(mlclNum[i])
            mol_beads = beads[moltype-1]
            if boxNum[i] == '1':
                coordinates_to_keep += [total_coordinates + j for j in range(1,mol_beads+1)]
                nmolty_keep[moltype-1] += 1
            total_coordinates += mol_beads
        return coordinates_to_keep, nmolty_keep
    
    for line in open(fort4):
        if ('=' in line):
            if 'nbox' in line:
                nbox_old = int(line.split()[-1])
            elif 'nchain' in line:
                nchain_old = int(line.split()[-1])
            elif 'nmolty' in line:
                nmolty_old = int(line.split()[-1])
                nmolty_present = [0 for i in range(nmolty_old)]


    f77new = open(path + '/fort.77','w')
    nline = 0
    mol_num_in_beads = 0 
    fluc_charge = False
    trans_rotation_info = 2+2*nbox_old*nmolty_old
    fluc_charge_start = trans_rotation_info + 1
    displ_start = 3
    for line in open(restart_file):
        nline += 1
        if nline < displ_start:
            f77new.write(line)
        elif nline <= trans_rotation_info:
            my_box = math.floor( (nline - displ_start)/(nmolty_old*2) ) + 1
            if my_box == 1:
                f77new.write(line)
            elif my_box == 2:
                if (nline - displ_start) % 2 == 2:
                    # line is odd, transl_displ is line
                    f77new.write(' %f %f %f\n'%(boxlx_vapor/2, boxlx_vapor/2, boxlx_vapor/2))
                else:
                    f77new.write(line)
        elif nline == fluc_charge_start:
            fluc_charge = True
            num_fluq_charge = len(line.split())
            f77new.write(line)
        elif fluc_charge:
            num_fluq_charge += len(line.split())
            if num_fluq_charge <= 2*nmolty_old:
                f77new.write(line)
            if num_fluq_charge == nbox_old*nmolty_old:
                fluc_charge = False
                volume_displ = nline + 1
                nchain_start = nline + 1 + nbox_old
        elif nline == volume_displ:
            f77new.write(' 1000.0   1000.0\n')
            ibox = 0
        elif nline <= nchain_start:
            ibox += 1
            if ibox == 1:
                f77new.write(line)
            elif ibox == 2:
                f77new.write(' {:f} {:f} {:f}\n'.format(boxlx_vapor,
                                                   boxlx_vapor,
                                                   boxlx_vapor))
        elif nline <= nchain_start + 2:
            if nline == nchain_start + 1:
                f77new.write('     NCHAIN\n')
            else:
                f77new.write(line)
            beads_section = True
        elif beads_section:
            mol_num_in_beads += len(line.split())
            if mol_num_in_beads == len(line.split()):
                beads_by_type = []
            beads_by_type += list(map(int, line.split()))
            f77new.write(line)
            if mol_num_in_beads == nmolty_old:
                beads_section = False
                moleculeNumbers = True
                molNum = []
            elif mol_num_in_beads > nmolty_old:
                print('something went wrong in reading beads by molty section')
                quit()
        elif moleculeNumbers:
            for i in line.split():
                molNum.append(i)
            if len(molNum) == nchain_old:
                moleculeNumbers = False
                boxNumbers = True
                boxNum = []
        elif boxNumbers:
            for i in line.split():
                boxNum.append(i)
            if len(boxNum) == nchain_old:
                boxNumbers = False
                coordinates = True
                nchain_new = writeNewBoxInfo(molNum, boxNum)
                coordinate_numbers, nchain_molty  = findCoordToKeep(molNum, boxNum, beads_by_type, nmolty_present)
                coordinates_start = nline + 1
                number_of_coordinates = 0
        elif coordinates:
            if ((nline - coordinates_start) %2) == 0:
                number_of_coordinates += 1
            if number_of_coordinates in coordinate_numbers:
                f77new.write(line)
    f77new.close()
    print(nchain_molty)
    with open(path + '/fort.77') as f:
        file_str = f.read()
        file_str = file_str.replace('NCHAIN', '%i'%sum(nchain_molty))
    with open(path + '/fort.77','w') as f:
        f.write(file_str)



    beads_for_calc = beads_by_type
    beads_for_calc[0] = 0
    cbmc_frac, transl_frac = calculateProbs(nchain_molty, beads_for_calc)
    x_tot = 0.
    mole_frac = []
    for x in nchain_molty:
        x_tot += x/sum(nchain_molty)
        mole_frac.append( x_tot)
    cbmc_dof = 0
    for k in range(len(beads_for_calc)):
        cbmc_dof += beads_for_calc[k]*nchain_molty[k]
    cbmc_molty_frac = [l*u/cbmc_dof for (l,u) in zip(beads_for_calc, nchain_molty)]
    f4new = open(path + '/fort.4','w')
    ibox = 0
    sim_box_end = False
    bp_section, swap_section = False, False
    for line in open(fort4):
        if ('=' in line):
            if 'nbox' in line:
                f4new.write(makeNewLine(indent, 'nbox','2'))
            elif 'nchain' in line:
                f4new.write(makeNewLine(indent,'nchain','%i'%nchain_new))
            elif 'time_limit' in line:
                f4new.write(makeNewLine(indent,'time_limit','6000'))
            elif 'pmvlmt' in line:
                f4new.write(makeNewLine(indent, 'pmvlmt','0.0d0 0.0d0'))
            elif ('pmvol' in line) or (' pmswat ' in line):
                f4new.write(makeNewLine(indent, line.split()[0], '0.0d0'))
            elif ('rmvolume' in line) or ('allow_cutoff_failure' in line):
                continue
            elif 'pmswap' in line:
                pswap = float(line.split()[-1].rstrip('d0'))
                f4new.write(line)
                pmcb = (1-pswap)*cbmc_frac + pswap
                pmtra = (1-pswap)*transl_frac + pswap
            elif ' pmcb ' in line:
                f4new.write(makeNewLine(indent, 'pmcb', '%5.4fd0'%pmcb))
            elif ' pmtra ' in line:
                f4new.write(makeNewLine(indent, 'pmtra', '%5.4fd0'%pmtra))
            elif (' pmtrmt' in line) or (' pmromt ' in line) or (' pmswmt ' in line):
                L = ' '.join(['%6.3fd0'%i for i in mole_frac])
                f4new.write(makeNewLine(indent, line.split()[0], L))
            elif ' pmcbmt ' in line:
                L = ' '.join(['%6.3fd0'%i for i in cbmc_molty_frac])
                f4new.write(makeNewLine(indent, line.split()[0], L))
            else:
                f4new.write(line)
        elif len(line.split()) == 13:
            ibox += 1
            if ibox <= 2:
                L = line.split()
                L[-2] = '%f'%temp
                if ibox == 2:
                    L[0], L[1], L[2], L[3] = ['%f'%boxlx_vapor, 
                                                '%f'%boxlx_vapor,
                                                '%f'%boxlx_vapor,
                                                '%f'%(boxlx_vapor/2)]
                f4new.write(' '.join(L) + '\n')
        elif (ibox == 1) or (ibox == 2):
            if (len(line.split()) == nmolty_old +1) and (line.split()[0][0] != '!'):
                if ibox == 1:
                    f4new.write(' '.join(['%i'%k for k in nchain_molty] + ['0']) + '\n')
                else:
                    L = ['0' for i in range(nmolty_old+2)]
                    f4new.write(' '.join(L) + '\n')
            else:
                f4new.write(line)
        elif ibox < 1:
            f4new.write(line)
        elif (ibox > 2) and ('SIMULATION_BOX' in line):
            f4new.write(line)
            sim_box_end = True
        elif line.startswith('MC_SWAP'):
            swap_section = True
            f4new.write(line)
            for i in range(nmolty_old):
                f4new.write('! nswapb pmswapb\n1 1.0d0\n! box1 box2\n1 2\n')
        elif line.startswith('END MC_SWAP'):
            swap_section = False
            f4new.write(line)
        elif line.startswith('UNIFORM_BIASING_POTENTIALS'):
            bp_section = True
            f4new.write(line)
            for i in range(nmolty_old):
                f4new.write('0.0d0 0.0d0\n')
        elif line.startswith('END UNIFORM_BIASING_POTENTIALS'):
            f4new.write(line)
            bp_section = False
        elif sim_box_end and (not swap_section) and (not bp_section):
            f4new.write(line)
    f4new.close()

indent = ' '*4
from writer import makeNewLine
from MCFlow.changeProbs import calculateProbs
import math


if __name__ == '__main__':
    import os
    main_dir = os.getcwd()
    if not os.path.exists(main_dir + '/jobs'):
        os.makedirs(main_dir + '/jobs')
    njobs = 0
    nfiles = 1
    f = open('jobs/run1.txt','w')
    for mol in ['192mol','128mol','96mol','64mol']:
        for T in [343,353,363,373,383,393,403,413,423,433,443,453,463,473,483,493]:
            for L in [24, 48, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,
                        32768,65536,131072,262144,524288]:
                for i in range(1, 33):
                    njobs += 1
                    if njobs%512 == 0:
                        f.close()
                        nfiles += 1
                        f = open('jobs/run%i.txt'%nfiles,'w')
                    path = '%s/%iK/%i/%i'%(mol,T, int(L), i)
                    if not os.path.exists(path):
                        os.makedirs(path)
                    f.write(path + '\n')
                    keep_zeo(main_dir + '/' + path,
                            'fort.4.%s'%mol,
                            'fort.77.%s'%mol, T, L)
