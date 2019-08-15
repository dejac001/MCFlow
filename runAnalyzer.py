class NoFilesAnalyzed(Exception):
    pass
import MCFlow.file_organization as fo

def findFinalSimRun(path, tag):
    nfiles = 0
    if ('prod' in tag) and ('butanol-water' in path) and ('323K' in path):
        while os.path.isfile('%s/%s%i/run.%s%i'%(path, tag, nfiles+1, tag, nfiles+1))\
                and os.path.isfile('%s/%s%i/movie_Hbond_moltype.txt'%(path, tag, nfiles+1)):
            nfiles += 1
    else:
        while os.path.isfile(fo.read(path,'run.',tag,nfiles+1)):
            nfiles += 1
    return nfiles

def findNextRun(path, tag):
    Runs = [i.split(tag)[-1] for i in os.listdir(path) if (i.startswith(tag) and len(i.split(tag))>1)]
    assert len(Runs) > 0, 'No runs found' + path
    Runs = [i for i in Runs if len(i) > 0]
    Runs = set(Runs)
    sorted_numbers = sorted([int(i) for i in Runs])
    sorted_numbers.reverse()
    for runNum in sorted_numbers:
        if os.path.isfile( fo.read(path, 'run.', tag, runNum) ):
            return runNum+1


def findFinalSimIntvl(path, fileNumStart, numIntvl, tag='equil-'):
    '''
    nIntvl is integer interval of number of files at end of finished runs
    '''
    fileNum = fileNumStart
    while os.path.isfile(fo.read(path,'run.',tag, fileNum)):
        fileNum += 1
    if fileNum <= numIntvl:
        print('interval for analysis too long, just going to analyze all runs')
        runAnalStart = fileNumStart
    elif fileNum > numIntvl:
        runAnalStart = fileNum - numIntvl
    return runAnalStart

def what2Analyze(path, oldTagName, fileNumStart, numIntvl):
    if '-' in oldTagName:
        TAG = oldTagName[:(oldTagName.find('-')+1)]
    else:
        TAG = oldTagName
    if 'prod' in oldTagName:
        if numIntvl == 0:
            # use all files for production
            old_begin = fileNumStart
            nfiles = findFinalSimRun(path, TAG)
        else:
            old_begin = findFinalSimIntvl(path, fileNumStart, numIntvl, TAG)
            nfiles = numIntvl
    elif 'equil' in oldTagName:
        if numIntvl == 0: numIntvl = 1
        old_begin = findFinalSimIntvl(path, fileNumStart, numIntvl, TAG)
        nfiles = numIntvl
    else:
        raise NoFilesAnalyzed
    return old_begin, nfiles

def getRealRho(rhoBias, bias, T):
    '''
    get real number densities: note: units in chain/nm**3
    '''
    rhoReal = {}
    for mol in rhoBias.keys():
        for box in rhoBias[mol].keys():
            if type(rhoBias[mol][box]) == type(list()):
                if mol not in rhoReal.keys(): rhoReal[mol] = {}
                rhoReal[mol][box] = [i*math.exp(bias[mol][box]/T) for i in rhoBias[mol][box]]
            elif rhoBias[mol][box] > 0.:
                if mol not in rhoReal.keys(): rhoReal[mol] = {}
                rhoReal[mol][box] = rhoBias[mol][box] * math.exp(  bias[mol][box]/T )
    assert rhoReal, ('All number densities = 0.0. It is likely that ' +
                'the previous simulation did not run to iblock cycles')
    return rhoReal



def calc95conf(stdev, numIndep):
    T_values = {'8':2.365,'4':3.182, '16':2.131, '32':2.04,'2':4.303,
                    '6':2.571,'7':2.447,'1':1.0, '5':2.776}
    assert '%i'%numIndep in T_values.keys(), 'No T-value stored for %i indep'%numIndep
    return stdev/math.pow(numIndep, 0.5)*T_values['%i'%numIndep]

def checkRun(tag, listOfDBs, feed):
    '''
    check to see if all DBs have data corresponding
    to the same run
    '''
    assert tag != None, 'No tag provided to find run'
    run_names = []
    for db in listOfDBs:
        run_names += list(db[feed].keys())
    data_run = ''
    for run in set(run_names):
        if tag in run:
            data_run += run
    assert data_run, 'No run analyzed of run type prompted for feed %s'%feed
    assert data_run in listOfDBs[0][feed].keys(), 'More than one run found to analyze {}'.format(listOfDBs[0][feed].keys())
    return data_run

def calcDGfromNumDens(rho, N_i, T):
    dG = {}
    for mlcl_name in N_i.keys():
        mol = mlcl_name.split()[0]
        if mol not in rho.keys(): continue
        if (N_i[mlcl_name] > 0):
            assert mol in rho.keys(), '{} not in rho.keys()!'.format(mlcl_name)
            dG[mol] = {}
            if mol not in rho.keys():
                continue
            for boxFrom in sorted(rho[mol].keys()):
                for boxTo in [i for i in list(rho[mol].keys()) if i != boxFrom]:
                    key = boxFrom + '--' + boxTo
                    try:
                        DELTAG =( -R['kJ/(mol*K)'] *T*math.log(rho[mol][boxTo])
                        - -R['kJ/(mol*K)'] *T*math.log(rho[mol][boxFrom]))
                        if key not in dG[mol].keys(): dG[mol][key] = []
                        dG[mol][key].append(DELTAG )
                    except ValueError:
                        print('math domain error for DG for mol {} from {} to {}'.format(mlcl_name, boxFrom, boxTo))
    return dG

def getRelMols(N, box):
    mols = []
    for mol in N.keys():
        means = {box:np.mean(value) for box,value in N[mol].items()}
        if (means[box] > 0.):# and (means[box] < sum(means.values())):
            mols.append(mol)
    return sorted(mols)

def getTotalMols(N):
    Ntotal = {}
    for box in N['1'].keys():
        mols = getRelMols(N, box)
        mols.reverse()
        N_total = np.zeros(len(N['1'][box]),dtype=int)
        for imol, mol in enumerate(mols):
            N_total = N_total + np.array(N[mol][box], dtype=int)
        Ntotal[box] = N_total
    return Ntotal

def getKX(N, Ntot):
    K = {}
    X = {}
    for box in N['1'].keys():
        K[box] = {}
        X[box] = {}
        mols = getRelMols(N,box)
        mols.reverse()
        if len(mols) <= 1:
            K[box] = { i:[0.] for i in mols }
            X[box] = {i:[1.] for i in mols}
            continue
        for imol, mol1 in enumerate(mols):
            for mol2 in mols:
                if mol1 == mol2: continue
                # calculate K
                K[box][mol1 +'/'+ mol2] = []
                for (i,j) in zip(N[mol1][box], N[mol2][box]):
                    if j > 0.:
                        K[box][mol1 +'/'+ mol2].append( i/j )
        for imol, mol1 in enumerate(mols):
            x_np = np.divide( N[mol1][box], Ntot[box])
            x_np = x_np[np.where(x_np < np.inf)]
            X[box][mol1] = x_np.tolist()
    return K, X

def getRhoTotal(dens):
    rho = {}
    for mol, vals in dens.items():
        for box, value in vals.items():
            if box not in rho.keys(): rho[box] = 0.
            if type(value) == type(dict()):
                rho[box] += value['mean']
            else:
                rho[box] += value
    return rho

def calcdHig_mixt(U, T, rhomol, Ntot):
    if 'box3' not in U.keys():
        p_box = 'box%i'%(len(U.keys()))
    else:
        p_box = 'box3'
    rho_total = getRhoTotal(rhomol)
    p = rho_total[p_box]/N_av*R['nm**3*kPa/(mol*K)']*T
    dH = {}
    boxes = sorted(list(U.keys()))
    boxes.reverse()
    for i, boxFrom in enumerate(boxes):
        for boxTo in [j for j in boxes[i:] if j != boxFrom]:
            key = boxFrom + '-->' + boxTo
            dH[key] = []
            for u_from, u_to, N_from, N_to in zip(
                U[boxFrom], U[boxTo],
                Ntot[boxFrom], Ntot[boxTo]
            ):
                if (N_to > 0) and (N_from > 0) and boxFrom in rho_total.keys() and boxTo in rho_total.keys():
                    v_from = 1000./rho_total[boxFrom]*N_av
                    v_to = 1000./rho_total[boxTo]*N_av
                    delta_H =( (u_to/N_to - u_from/N_from)*R['kJ/(mol*K)'] +
                           p*(v_from - v_to)/R['\AA**3*kPa/(mol*K)']*R['kJ/(mol*K)'])
                    dH[key].append(delta_H)
    return dH

def calcdH_mixt(U, P, rhomol, Ntot):
    if 'box3' not in U.keys():
        p_box = 'box%i'%(len(U.keys()))
    else:
        p_box = 'box3'
    rho_total = getRhoTotal(rhomol)
    dH = {}
    boxes = sorted(list(U.keys()))
    boxes.reverse()
    for i, boxFrom in enumerate(boxes):
        for boxTo in [j for j in boxes[i:] if j != boxFrom]:
            key = boxFrom + '-->' + boxTo
            dH[key] = []
            for p, u_from, u_to, N_from, N_to in zip(P[p_box],
                U[boxFrom], U[boxTo],
                Ntot[boxFrom], Ntot[boxTo]
            ):
                if (N_to > 0) and (N_from > 0) and boxFrom in rho_total.keys() and boxTo in rho_total.keys():
                    v_from = 1000./rho_total[boxFrom]*N_av
                    v_to = 1000./rho_total[boxTo]*N_av
                    delta_H =( (u_to/N_to - u_from/N_from)*R['kJ/(mol*K)'] +
                           p*(v_from - v_to)/R['\AA**3*kPa/(mol*K)']*R['kJ/(mol*K)'])
                    dH[key].append(delta_H)
    return dH

def rhoMCC(nMol, boxLx):
    number_densities = {}
    for mol in nMol.keys():
        number_densities[mol] = {}
        for box, val in nMol[mol].items():
            V = [i**3/1000 for i in boxLx[box]]
            number_densities[mol][box] = [j/k for j,k in zip(val, V)]
    return number_densities

def getRhoByMCC(feed, seed, path, type, guessStart, interval,
                verbosity, **kwargs):
    from MCFlow.file_formatting import reader
    my_dir = '%s/%s/%i' % (path, feed, seed)
    (old_begin, nfiles) = what2Analyze(my_dir, type,
                                       guessStart, interval)
    if verbosity > 0:
        print('old_begin = {}, nfiles = {}'.format(old_begin, nfiles))
    # get data from old files
    (N, P, boxLengths,
     ncycle_old, molWeights, E) = reader.read_fort12(my_dir, old_begin,
                                                     nfiles, tag=type)

    (number_densities, chemical_potential,
     swap_info, biasPot, volumes,
     totalComposition, cbmc_info,
     T, zeolite) = reader.go_through_runs(my_dir, ncycle_old,
                                          old_begin, nfiles,
                                          tag=type)
    number_densities = rhoMCC(N, boxLengths)
    # do calculations for other data that may be needed
    number_dens_real = getRealRho(number_densities, biasPot, T)
    os.chdir(path)
    return number_dens_real


def getFileData(feeds, indep, path, type, guessStart, interval,
                verbosity, liq=False, mol='-1', **kwargs):
    from MCFlow.file_formatting import reader
    general_data = {key:{} for key in feeds}
    # TODO: add molec name (i.e. 15PDO or WATER) into general_data
    for feed in feeds:
        if verbosity > 0: print('-'*12 + 'Dir is %s'%feed + '-'*12)
        nNotFound = 0
        for seed in indep:
            try:
                my_dir = '%s/%s/%i'%(path,feed,seed)
                (old_begin, nfiles) = what2Analyze(my_dir, type,
                                                           guessStart,interval)
                if verbosity > 0:
                    print('old_begin = {}, nfiles = {}'.format(old_begin, nfiles))
                # get data from old files
                (N, P, boxLengths,
                 ncycle_old, molWeights, E) = reader.read_fort12(my_dir, old_begin,
                                                                 nfiles, tag=type)
                Ntotal = getTotalMols(N)

                K, mole_frac = getKX(N, Ntotal)
                (number_densities, chemical_potential,
                 swap_info, biasPot, volumes,
                 totalComposition, cbmc_info,
                 T, zeolite) = reader.go_through_runs(my_dir, ncycle_old,
                                                                          old_begin, nfiles,
                                                                          tag=type)
                # do calculations for other data that may be needed
                number_dens_real = getRealRho(number_densities, biasPot, T)
                deltaG = calcDGfromNumDens(number_dens_real, totalComposition, T)
                if liq:
                    concentrations = {}
                    assert kwargs['box'], 'box needed for liquid phase'
                    if 'box' not in kwargs['box']: kwargs['box'] = 'box' + kwargs['box']
                    l_box = kwargs['box']
                    c = calc_tools.g_mL(N[mol][l_box], boxLengths[l_box],
                                            MW=molWeights[mol])
                    concentrations[mol] = {l_box:c}
                # initialize vars
                if (seed == indep[0]) and (feed == feeds[0]):
                    k_ratio = properties.AnyProperty(K)
                    X = properties.AnyProperty(mole_frac)
                    boxlx = properties.AnyProperty(boxLengths)
                    CBMC = properties.AnyProperty(cbmc_info)
                    Press = properties.AnyProperty(P)
                    Nmlcl = properties.AnyProperty(N)
                    SWAP = properties.AnyProperty(swap_info)
                    U = properties.AnyProperty(E)
                    rho = properties.AnyProperty(number_dens_real)
                    dG = properties.AnyProperty( deltaG )
                    data = {'CBMC':CBMC, 'P':Press, 'N':Nmlcl,
                            'SWAP':SWAP, 'U':U, 'rho':rho,
                            'boxlx':boxlx, 'dG':dG, 'K':k_ratio,
                             'X':X}
                    if liq:
                        if (verbosity > 1):
                            print('Doing analysis for C [ g/mL ] for mol {}'.format(mol))
                            print('Box 2 should be liquid phase')
                        C = properties.AnyProperty( concentrations )
                        data['Conc'] = C
                CBMC.addVals(cbmc_info)
                Press.addVals(P)
                Nmlcl.addVals(N)
                SWAP.addVals(swap_info)
                U.addVals(E)
                rho.addVals(number_dens_real)
                boxlx.addVals(boxLengths)
                dG.addVals(deltaG)
                k_ratio.addVals(K)
                X.addVals(mole_frac)
                if liq:
                    C.addVals(concentrations)
            except FileNotFoundError:
                nNotFound += 1
                print('File not found for dir {}'.format(my_dir))
                print('Not averaging for this directory')
                continue
        if nNotFound == len(indep): raise NoFilesAnalyzed
        for cls in data.values():
            cls.avgVals(feed)
        general_data[feed]['numIndep'] = len(indep)
        general_data[feed]['indepSims'] = list(indep)

        # get general data
        general_data[feed]['compositions'] = totalComposition
        general_data[feed]['ncycle'] = ncycle_old
        general_data[feed]['molecular weight'] = molWeights
        general_data[feed]['temperature'] = T
        if zeolite:
            general_data[feed]['zeolite'] = zeolite
    os.chdir(path)
    return data, general_data

import math, os, sys
import numpy as np
#np.seterr(all='raise')
try:
    from MCFlow import properties, calc_tools, file_formatting
    from MCFlow.chem_constants import R, N_av
except ImportError:
    print('import error')
    print(sys.path)
