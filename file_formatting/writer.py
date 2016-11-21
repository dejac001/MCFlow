def makeNewLine(pad,name,val):
    #:: val is a string
    newLine = '%s%-15s = %s\n'%(pad,name,val)
    return newLine

def xyz(file_name, DATA):
    '''
    for arbitrary number of atoms
    DATA = type(dict)
    DATA.keys() == ['atoms','coords']
    DATA['coords'] = [[x0,y0,z0],...[xN,yN,zN]]

    writes coordinates/atoms in same order as was read
    '''
    if len(DATA['atoms']) != len(DATA['coords']):
        print('number of atoms doesnt equal number of coordinates')
        for key in ['atoms','coords']:
            print('%s:  '%key, len(DATA[key]))
        quit()
    total_atoms = len(DATA['atoms'])
    f = open(file_name, 'w')
    f.write('%i\n\n' % total_atoms)
    for el, coords in zip(DATA['atoms'],DATA['coords']):
        f.write('%-7s%-12.4f%-12.4f%-12.4f\n' % (el, coords[0], coords[1], coords[2]))
    f.close()

def write_new_bias(PATH, old_begin, nfiles, boxlx_vapor, biasPotNew, nGhost,
                   nindep, rcut_vapor, nstepnew = 25000, tag='equil-'):
    nMolTy = 100000
    for i in nindep:
        os.chdir(PATH + '/%d' % i)
        f = open('fort.4')
        fnew = open('fort.4.' + tag + str(old_begin + nfiles),'w')
        box = 0
        flag = 0
        for line in f:
            if 'nstep' in line:
                fnew.write(makeNewLine(indent, 'nstep','%i'%nstepnew))
            elif 'iprint' in line:
                fnew.write(makeNewLine(indent, 'iprint','%i'%int(nstepnew/10)))
            elif 'iblock' in line:
                iblock = int(nstepnew/10)
                if iblock > 1000: iblock = 1000
                fnew.write(makeNewLine(indent, 'iblock','%i'%iblock))
            elif ('nmolty' in line) and ('nmolty' == line.split()[0]):
                nMolTy = int(line.split()[-1])
                fnew.write(line)
            elif (len(line.split()) == 13) and (' F ' in line):
                box += 1
                if box == 3: # in vapor box
                    if line.split()[-4] == 'T':
                        print('vapor box was ideal gas',os.getcwd(),'-- keeping ideal')
                        my_rcut = 14.0
                    else:
                        my_rcut = boxlx_vapor*rcut_vapor
                    L = line.split()
                    L[0] = '%7.3f'%boxlx_vapor
                    L[1] = '%7.3f'%boxlx_vapor
                    L[2] = '%7.3f'%boxlx_vapor
                    L[3] = '%7.2fd0'% my_rcut
                    L[-4] = 'T' # ideal gas (for now)
                    newline = ' '.join( i for i in L) + '\n'
                    fnew.write(newline)
                else:
                    fnew.write(line)
            elif ((len(line.split()) == (nMolTy+1)) and (line.split()[-1].isnumeric())
                and (box == 3)):
                nMol = list(map(int,line.split()))
                nMol[-1] = nGhost
                newline = ' '.join('%i'%i for i in nMol) + '\n'
                fnew.write(newline)
            elif line.startswith('UNIFORM_BIASING_POTENTIALS'):
                fnew.write('UNIFORM_BIASING_POTENTIALS\n')
                mols = list(map(int, biasPotNew.keys()))
                if len(mols) > nMolTy: mols = mols[:nMolTy]
                ordered_mols = list(map(str, sorted(mols)))
                print('ordered mols for biasing potentials is {}'.format(ordered_mols))
                for mol in ordered_mols:
                    for box in sorted(biasPotNew[mol].keys()):
                        fnew.write('%.2fd0 '%(biasPotNew[mol][box]))
                    fnew.write('\n')
                fnew.write('END UNIFORM_BIASING_POTENTIALS\n\n\n')
                fnew.write('SPECIFIC_ATOM_TRANSL\n')
                fnew.write('! How many atoms should we be performing translations on?\n')
                fnew.write('! What is the index of those atoms in their respective molecules?\n')
                fnew.write('! What is the molecule type for each atom?\n')
                fnew.write('END SPECIFIC_ATOM_TRANSL\n')
                fnew.close()
                flag = 1
            elif flag != 1:
                fnew.write(line)
        f.close()
        fnew.close()

def write_new_prob(PATH, pswapMolTy, pswatchMolTy, pmvol, pmswap, pswatch,
                   composition, nindep, boxLengths, rcut_vapor, time, nstepnew=25000, tag='equil-',Pdir1=False):
    import runAnalyzer
    def writePMT(probs):
        val = ''
        for pmt in probs:
            if pmt < 1e-08:
                val += '0.0d0  '
            elif pmt < 1e-04:
                val += '%6.4e ' % pmt
            else:
                val += '%6.4fd0  ' %pmt
        return val
    for i in nindep:
        os.chdir(PATH + '/%i' % i)
        f = open('fort.4')
        newRun = runAnalyzer.findNextRun(PATH + '/%i'%i, tag)
        fnew = open('fort.4.' + tag + str(newRun),'w')
        box = 0
        numSwatchTy = -1
        val_pmsatc = ''
        for line in f:
            # write customary variables
            if 'nstep' in line:
                fnew.write(makeNewLine(indent,'nstep','%-7i' % nstepnew))
                fnew.write(makeNewLine(indent,'time_limit','%-7i' % time))
            elif 'time_limit' in line:
                continue
            elif 'rmin ' in line:
                fnew.write(makeNewLine(indent,'rmin','1.0d0'))
            elif 'seed' in line:
                fnew.write(makeNewLine(indent,'seed','%i'%i))
            elif 'lstop' in line:
                fnew.write(makeNewLine(indent,'lstop','F'))
            elif ' nmolty' in line:
                molReadInFort4 = int(line.split()[-1])
                fnew.write(line)
            elif 'iratio' in line:
                fnew.write(makeNewLine(indent,'iratio',500))
            elif line.find('iprint') != -1:
                fnew.write(makeNewLine(indent,'iprint','%-7i' % int(nstepnew/10)))
            elif line.find('iblock') != -1:
                iblock = int(nstepnew/10)
                if iblock > 1000: iblock = 1000
                fnew.write(makeNewLine(indent,'iblock','%-7i' % iblock))
            elif line.find('pmswmt') != -1:
                value = writePMT(pswapMolTy)
                if len(pswapMolTy) < molReadInFort4:
                    print('extraneous molecules in fort.4... take lots of care to check'
                          'that probability changer is putting correct values for pmswmt...'
                          '(it\'s likely wrong)')
                fnew.write(makeNewLine(indent,'pmswmt',value))
            elif 'nswaty' in line:
                numSwatchTy = int(line.split()[-1])
                fnew.write(line)
                value = writePMT(pswatchMolTy)
                for extraMols in range(numSwatchTy-len(pswatchMolTy)):
                    value += '0.00d0  '
                fnew.write(makeNewLine(indent,'pmsatc',value))
            elif ' pmsatc ' in line:
                continue
            elif 'pmvol' in line and ('pmvolx' not in line)\
                    and ('pmvoly' not in line) and ('pmvolb' not in line):
                if pmvol <= 0.:
                    print('-----NO VOLUME MOVES WILL BE DONE!!!-----')
                value = line.split()[-1]
                if 'd' in value:
                    number = float(value[:value.find('d')])
                else:
                    number = float(value)
                if number <=0.:
                    print(' no volume moves done previously. You doing NVT with gas and zeolite?')
                    print(' keeping pmvol = 0')
                    pmvol = 0.
                fnew.write(makeNewLine(indent, 'pmvol','%0.5fd0'%pmvol))
            elif 'tavol' in line:
                fnew.write(makeNewLine(indent, 'tavol', '0.40d0'))
            elif 'rmvolume' in line:
                fnew.write(makeNewLine(indent, 'rmvolume', line.split()[-1]))
                if 'equil' in tag: fnew.write(makeNewLine(indent, 'allow_cutoff_failure', '1'))
            elif 'allow_cutoff_failure' in line:
                continue
            elif 'pmswat' in line:
                value = line.split()[-1]
                if 'd' in value:
                    number = float(value[:value.find('d')])
                else:
                    number = float(value)
                if number <=0.:
                    print(' no swatch moves done previously')
                    print(' keeping pmswat = 0')
                    pswatch = 0.
                fnew.write(makeNewLine(indent, 'pmswat', '%0.5fd0'%pswatch))
            elif (' pmswap ' or 'pmswap=') in line:
                swapold = line.split()[-1]
                if 'd' in swapold:
                    number = swapold[:swapold.find('d')]
                elif 'e' in swapold:
                    number = swapold
                pmswap_old = float(number)
                if pmswap_old == 0.:
                    pmswap_old = pmvol
                fnew.write(makeNewLine(indent,'pmswap','%5.4fd0'%pmswap))
            elif ('pmcb ' or 'pmcb=') in line:
                # we want to update pmcb here so that the correct ratio of pmcb/pmtra is maintained
                value = line.split()[-1]
                if 'd' in value:
                    number = value[:value.find('d')]
                else:
                    number = value
                pmcb_old = float(number) - pmswap_old
                if pmcb_old <= 0.:
                    print('pmcb was zero in default file, keeping zero')
                    pmcb_old = 0
                    pmtra_old = (1-pmswap_old)/2
                    pmcb = 0
                else:
                    pmcb = (1-pmswap)*pmcb_old/(1-pmswap_old) + pmswap
                    pmtra_old = (1- pmswap_old - pmcb_old)/2
                if pmcb <= pmswap:
                    print('pmcb = {} is less than pmswap = {}'.format(pmcb, pmswap))
                fnew.write(makeNewLine(indent,'pmcb','%5.4fd0'%pmcb))
            elif 'pmtra' in line:
                if pmcb != 0.:
                   pmtra = (1-pmswap)*pmtra_old/(1-pmswap_old) + pmcb
                else:
                   pmtra = (1-pmswap)*pmtra_old/(1-pmswap_old) + pmswap
                if pmtra > 1:
                    print('pmtra greater than 1!!!')
                    print('pmswap, pmtra_old, pmswap_old, pmcb, pmcb_old ='+
                            '{}, {}, {}, {}, {}'.format(pmswap, pmtra_old,
                                                    pmswap_old, pmcb, pmcb_old))
                    quit()
                fnew.write(makeNewLine(indent,'pmtra','%5.4fd0'%pmtra))
            elif (len(line.split()) == 13) and ('=' not in line): # simulation box
                box += 1
                L = line.split()
                boxlx = boxLengths['box%i'%box]['mean']
                L[0] = '%7.3f'%boxlx
                L[1] = '%7.3f'%boxlx
                L[2] = '%7.3f'%boxlx
                if L[-4] == 'T':
                    lideal= True
                elif L[-4] == 'F':
                    lideal = False
                else:
                    print('error reading simulation box info')
                    print(line)
                if (boxlx > 100) and (not lideal):
                    rcut = boxlx*rcut_vapor
                else:
                    rcut = 14.0000
                if rcut_vapor < 1.0: L[3] = '%7.2fd0'% (rcut)
                newline = ' '.join( i for i in L)
                newline += '\n'
                fnew.write(newline)
            elif Pdir1 and line.find('MC_SWAP') != -1 and line.find('END') == -1:
                # haven't figured out how to put swap directions in yet
                fnew.write(line)
                fnew.write(f.readline())
                fnew.write('2 %0.4fd0 1.0d0\n' % (Pdir1))
                f.readline() # go over line above in iteration
            else:
                fnew.write(line)
        f.close()
        fnew.close()

def sort_keys(keys):
    my_sort = []
    nums = []
    for key in keys:
        my_str = ''
        my_name = ''
        for letter in key:
            if letter.isdigit():
                my_str += letter
            else:
                my_name += letter
        nums.append(int(my_str))
    for num in sorted(nums):
        for key in keys:
            if key == my_name + '%i'%num:
                my_sort.append(key)
    return my_sort

def write_restart(data, newfile):
    f = open(newfile,'w')
    f.write(data['number of cycles'])
    f.write(data['max displacement']['atom translation'])
    for box in sort_keys(data['max displacement']['translation'].keys()):
        for molty in sort_keys(data['max displacement']['translation'][box].keys()):
            f.write(data['max displacement']['translation'][box][molty])
            f.write(data['max displacement']['rotation'][box][molty])
    for box in sort_keys(data['max displacement']['fluctuating charge'].keys()):
        for molty in data['max displacement']['fluctuating charge'][box].keys():
            f.write(' ' + data['max displacement']['fluctuating charge'][box][molty] )
        f.write('\n')
    max_displ_volume = ''
    for box in sort_keys(data['max displacement']['volume'].keys()):
        max_displ_volume += ' ' + data['max displacement']['volume'][box]
    f.write(max_displ_volume + '\n')
    for box in sort_keys(data['box dimensions'].keys()):
        if len(data['box dimensions'][box].split()) == 9:
            f.write('0.000000000000000000E+00 '*9 + '\n')
        f.write(data['box dimensions'][box])
    nchain = int(data['nchain'].rstrip('\n'))
    f.write(data['nchain'])
    f.write(data['nmolty'])
    index = 0
    for molty in sort_keys(data['nunit'].keys()):
        index += 1
        f.write('%12s'%data['nunit'][molty])
        if (index > 0) and ((index % 6) == 0):
            f.write('\n')
    if ((index % 6) != 0): f.write('\n')
    if ((len(data['mol types']) != nchain) or (len(data['mol types'])
                                              != len(data['box types']))):
        print('error in writing restart info')
        print('number of molecules not consistent')
        print('error: nchain: ',nchain)
        print('error: mol types: ', len(data['mol types']))
        print('error: box types: ', len(data['box types']))
        print('error: coords: ', len(data['coords']))
        print('quitting...')
        quit()
    for identity in [data['mol types'], data['box types']]:
        index = 0
        for chain in identity:
            index += 1
            f.write('%12s'%chain)
            if (index > 0) and ((index % 6) == 0):
                f.write('\n')
        if ((index % 6) != 0): f.write('\n')
    for mol in data['coords']:
        for bead in mol:
            f.write(bead['xyz'])
            f.write(bead['q'])
    f.close()

def write_fort4(data, newfile):
    indent = ' '*4
    f = open(newfile,'w')
    namelists = [i for i in data.keys() if i.startswith('&')]
    sections = [i for i in data.keys() if i not in namelists]
    variable_order = {'&mc_shared':('seed','nbox','nmolty','nchain','nstep', 'time_limit',
                                    'iratio','rmin','softcut','linit','lreadq'),
                      '&analysis':('iprint','imv','iblock','iratp','idele',
                                    'iheatcapacity','ianalyze'),
                      '&mc_volume':('tavol','iratv','pmvlmt','nvolb','pmvolb',
                                    'pmvol','pmvolx','pmvoly','rmvol',
                                    'allow_cutoff_failure'),
                      '&mc_swatch':('pmswat','nswaty','pmsatc'),
                      '&mc_swap':('pmswap','pmswmt'),
                      '&mc_cbmc':('rcutin','pmcb','pmcbmt','pmall','nchoi1','nchoi',
                                  'nchoir','nchoih','nchtor','nchbna','nchbnb','icbdir','icbsta',
                                  'rbsmax','rbsmin','avbmc_version','first_bead_to_swap',
                                  'pmbias','pmbsmt','pmbias2','pmfix','lrig','lpresim',
                                  'iupdatefix'),
                      '&mc_simple':('armtra','rmtra','rmrot','tatra','tarot','pmtra',
                                    'pmtrmt','pmromt')}
    namelist_order = ['&mc_shared','&analysis','&mc_volume','&mc_swatch','&mc_swap',
                     '&mc_cbmc','&mc_simple']
    section_order = ['SIMULATION_BOX','MOLECULE_TYPE','SAFE_CBMC','MC_SWAP','MC_SWATCH',
                     'INTERMOLECULAR_EXCLUSION','INTRAMOLECULAR_OH15','UNIFORM_BIASING_POTENTIALS',
                     'SPECIFIC_ATOM_TRANSL']
    if len(namelists) != len(namelist_order):
        print('amount of namelists not provided correctly')
        quit()
    elif len(sections) != len(section_order):
        print('Not enough sections provided')
        print('- missing sections: ',[i for i in section_order if i not in sections])
        print('- proceeding without this section')
        section_order = [i for i in section_order if i in sections]
    for NL in namelist_order:
        f.write(NL + '\n')
        if len(variable_order[NL]) > len(data[NL].keys()):
            print('Expected more variables for namelist: {} '
                  'including {}'.format(NL, [i for i in variable_order[NL] if i not in data[NL].keys()]))
        elif len(variable_order[NL]) < len(data[NL].keys()):
            print('More variables from previous input than expected for namelist {} '
                  'including {}'.format(NL, [i for i in data[NL].keys() if i not in variable_order[NL]]))

        for variable in [i for i in variable_order[NL] if i in data[NL].keys()]:
            value = data[NL][variable]
            if type(value) == type(''):
                f.write(makeNewLine(indent,variable,value))
            elif type(value) == type({}):
                val = ''
                for typ in sort_keys(value):
                    val += value[typ] + ' '
                f.write(makeNewLine(indent,variable,val.rstrip(' ')))
        f.write('/\n\n')
    for SEC in section_order:
        f.write(SEC + '\n')
        if SEC == 'SIMULATION_BOX':
            for box in sort_keys(data[SEC].keys()):
                f.write('! boxlx   boxly   boxlz   rcut  kalp   rcutnn numDimensionIsIstropic lsolid lrect lideal ltwice temperature pressure(MPa)\n')
                for var in ['dimensions', 'rcut', 'defaults', 'temperature', 'pressure']:
                    f.write(data[SEC][box][var] + ' ')
                f.write('\n')
                f.write('! nchain_1 ... nchain_nmolty ghost_particles\n')
                if len([i for i in data[SEC][box].keys() if 'mol' in i]) != int(data['&mc_shared']['nmolty']):
                    print('error in box molecule specs')
                    print(box)
                    print([i for i in data[SEC][box].keys() if 'mol' in i],  int(data['&mc_shared']['nmolty']))
                    print('newfile going to be: ',newfile)
                    quit()
                for molnum in sort_keys([i for i in data[SEC][box].keys() if 'mol' in i]):
                    f.write(data[SEC][box][molnum] + ' ')
                f.write(data[SEC][box]['nghost'] + '\n')
                f.write('! inix iniy iniz inirot inimix zshift dshift use_linkcell rintramax\n')
                f.write(data[SEC][box]['initialization data'])
                f.write('\n')
        elif SEC in ['MOLECULE_TYPE','MC_SWATCH']:
            for itype in sort_keys(data[SEC].keys()):
                f.write(data[SEC][itype] + '\n')
        elif SEC == 'MC_SWAP':
            for itype in sort_keys(data[SEC].keys()):
                f.write('! nswapb pmswapb\n')
                f.write('%i '%data[SEC][itype]['nswapb'] +
                        ' '.join(['%6.4fd0'%i for i in data[SEC][itype]['pmswapb']]) + '\n')
                f.write('! box1 box2\n')
                for boxpair in data[SEC][itype]['box1 box2']:
                    f.write(' '.join(['%i'%i for i in boxpair]) + '\n')
        elif SEC == 'INTERMOLECULAR_EXCLUSION':
            for mol1, vals1 in data[SEC].items():
                for unit1, mols2 in vals1.items():
                    for mol2, units2 in mols2.items():
                        for unit2 in units2:
                            f.write('%s %s %s %s\n'%(mol1.strip('mol'), unit1.strip('unit'),
                                                mol2.strip('mol'), unit2.strip('unit')))
        elif SEC == 'INTRAMOLECULAR_OH15':
            for itype in sort_keys(data[SEC].keys()):
                for intra in data[SEC][itype]:
                    f.write(itype.strip('mol') + intra)
        elif SEC == 'UNIFORM_BIASING_POTENTIALS':
            for molnum in sort_keys(data[SEC].keys()):
                for box in sort_keys(data[SEC][molnum].keys()):
                    f.write(data[SEC][molnum][box] + ' ')
                f.write('\n')
        f.write('END ' + SEC + '\n\n')
    f.close()


indent = ' '*4

import os
