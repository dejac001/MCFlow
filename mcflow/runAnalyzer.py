import file_organization as fo
from mcflow.file_formatting import reader
import math
import os
import numpy as np
import properties
import calc_tools
from chem_constants import R, N_av


class NoFilesAnalyzed(Exception):
    pass


# ----finding simulation paths/names--------------------------------------------------------
def findFinalSimRun(path, tag):
    nfiles = 0
    if ('prod' in tag) and ('butanol-water' in path) and ('323K' in path):
        while os.path.isfile('%s/%s%i/run.%s%i' % (path, tag, nfiles + 1, tag, nfiles + 1)) \
                and os.path.isfile('%s/%s%i/movie_Hbond_moltype.txt' % (path, tag, nfiles + 1)):
            nfiles += 1
    else:
        while os.path.isfile(fo.read(path, 'run.', tag, nfiles + 1)):
            nfiles += 1
    assert nfiles > 0, 'No files found!'
    return nfiles


def findNextRun(path, tag):
    """

    :param path: path to files
    :type path: str
    :param tag: tag associated with run file (prod- or equil- associated with production or equilibration files
    :type tag: str
    :return:
    """
    Runs = [i.split('-')[-1] for i in os.listdir(path) if tag in i]
    assert len(Runs) > 0, 'No runs found' + path
    Runs = [i for i in Runs if len(i) > 0]
    Runs = set(Runs)
    sorted_numbers = sorted([int(i) for i in Runs])
    sorted_numbers.reverse()
    for runNum in sorted_numbers:
        if os.path.isfile(fo.read(path, 'run.', tag, runNum)):
            return runNum + 1


def findFinalSimIntvl(path, fileNumStart, numIntvl, tag='equil-'):
    '''
    nIntvl is integer interval of number of files at end of finished runs
    '''
    fileNum = fileNumStart
    while os.path.isfile(fo.read(path, 'run.', tag, fileNum)):
        fileNum += 1
    if fileNum <= numIntvl:
        print('interval for analysis too long, just going to analyze all runs')
        runAnalStart = fileNumStart
    elif fileNum > numIntvl:
        runAnalStart = fileNum - numIntvl
    return runAnalStart


def what2Analyze(path, oldTagName, fileNumStart, numIntvl):
    """figure out which files to analyze

    :param oldTagName: tag name for file
    :type oldTagName: tag name for file
    :param fileNumStart: guess start for file num
    :param numIntvl: integer number of intervale, defaults to 0 in parsers
    """
    if '-' in oldTagName:
        # FIXME: is there a reason for this?
        TAG = oldTagName[:(oldTagName.find('-') + 1)]
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
    assert data_run, 'No run analyzed of run type prompted for feed %s' % feed
    assert data_run in listOfDBs[0][feed].keys(), 'More than one run found to analyze {}'.format(
        listOfDBs[0][feed].keys())
    return data_run


# ----finding number densities / performing calculations with them-----------------------------------------------------
def getRealRho(rhoBias, bias, T):
    """get real number densities: note: units in chain/nm**3"""
    rhoReal = {}
    for mol in rhoBias.keys():
        for box in rhoBias[mol].keys():
            if type(rhoBias[mol][box]) == type(list()):
                if mol not in rhoReal.keys(): rhoReal[mol] = {}
                rhoReal[mol][box] = [i * math.exp(bias[mol][box] / T) for i in rhoBias[mol][box]]
            elif rhoBias[mol][box] > 0.:
                if mol not in rhoReal.keys(): rhoReal[mol] = {}
                rhoReal[mol][box] = rhoBias[mol][box] * math.exp(bias[mol][box] / T)
    assert rhoReal, ('All number densities = 0.0. It is likely that ' +
                     'the previous simulation did not run to iblock cycles')
    return rhoReal


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


# ---------------calculation of numbers of molecules-----------------------------------
def getRelMols(N, box):
    mols = []
    for mol in N.keys():
        means = {box: np.mean(value) for box, value in N[mol].items()}
        if (means[box] > 0.):  # and (means[box] < sum(means.values())):
            mols.append(mol)
    return sorted(mols)


def getTotalMols(N):
    Ntotal = {}
    for box in N['1'].keys():
        mols = getRelMols(N, box)
        mols.reverse()
        N_total = np.zeros(len(N['1'][box]), dtype=int)
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
        mols = getRelMols(N, box)
        mols.reverse()
        if len(mols) <= 1:
            K[box] = {i: [0.] for i in mols}
            X[box] = {i: [1.] for i in mols}
            continue
        for imol, mol1 in enumerate(mols):
            for mol2 in mols:
                if mol1 == mol2: continue
                # calculate K
                K[box][mol1 + '/' + mol2] = []
                for (i, j) in zip(N[mol1][box], N[mol2][box]):
                    if j > 0.:
                        K[box][mol1 + '/' + mol2].append(i / j)
        for imol, mol1 in enumerate(mols):
            x_np = np.divide(N[mol1][box], Ntot[box])
            x_np = x_np[np.where(x_np < np.inf)]
            X[box][mol1] = x_np.tolist()
    return K, X


# ---------------calculation of transfer free energies or delta energies-----------------------------------
def calcDGfromNumDens(rho, N_i, T):
    """Calculation of transfer free energies; only do one direction to save storage space"""
    dG = {}
    for mlcl_name in N_i.keys():
        mol = mlcl_name.split()[0]
        if mol not in rho.keys(): continue
        if (N_i[mlcl_name] > 0):
            assert mol in rho.keys(), '{} not in rho.keys()!'.format(mlcl_name)
            dG[mol] = {}
            if mol not in rho.keys():
                continue
            boxes = sorted(rho[mol].keys(), reverse=True)
            for i, boxFrom in enumerate(boxes):
                for boxTo in boxes[(i + 1):]:
                    key = boxFrom + '--' + boxTo
                    try:
                        DELTAG = (-R['kJ/(mol*K)'] * T * math.log(rho[mol][boxTo])
                                  - -R['kJ/(mol*K)'] * T * math.log(rho[mol][boxFrom]))
                        if key not in dG[mol].keys(): dG[mol][key] = []
                        dG[mol][key].append(DELTAG)
                    except ValueError:
                        print('math domain error for DG for mol {} from {} to {}'.format(mlcl_name, boxFrom, boxTo))
    return dG


def calc_dU_dH(U, P, Ntot, box_length_x, box_length_y, box_length_z):
    """

    :param U: Internal energy from fort12 file (units of K*nchain)
    :param P: Pressures from fort12 file (units of kPa)
    :param Ntot: total number of chains in each box
    :param box_length_x: length of box in x-direction in Angstrom
    :param box_length_y: length of box in z-direction in Angstrom
    :param box_length_z: length of box in z-direction in Angstrom
    :return: dH, dU
    """
    if 'box3' not in U.keys():
        p_box = 'box%i' % (len(U.keys()))
    else:
        p_box = 'box3'
    dH = {}
    dU = {}
    boxes = sorted(list(U.keys()), reverse=True)
    for i, boxFrom in enumerate(boxes):
        for boxTo in boxes[(i + 1):]:
            key = boxFrom + '-->' + boxTo
            dH[key] = []
            dU[key] = []
            for i in range(len(box_length_x[boxTo])):
                N_from, N_to = Ntot[boxFrom][i], Ntot[boxTo][i]
                if N_to > 0 and N_from > 0:
                    u_from, u_to = U[boxFrom][i], U[boxTo][i]
                    delta_U = (u_to / N_to - u_from / N_from) * R['kJ/(mol*K)']
                    dU[key].append(delta_U)
                    p = P[p_box][i]
                    if np.isnan(p):
                        # pressure was not calculated, no need to calculate dH
                        continue
                    V_from = box_length_x[boxFrom][i] * box_length_y[boxFrom][i] * box_length_z[boxFrom][i]
                    V_to = box_length_x[boxTo][i] * box_length_y[boxTo][i] * box_length_z[boxTo][i]
                    molar_volume_from, molar_volume_to = V_from / N_from * N_av, V_to / N_to * N_av
                    delta_H = delta_U + p * (molar_volume_to - molar_volume_from) / R['\AA**3*kPa/(mol*K)'] * R[
                        'kJ/(mol*K)']
                    dH[key].append(delta_H)
    return dU, dH


def getFileData(feeds, indep, path, type, guessStart, interval,
                verbosity, liq=False, mol='-1', **kwargs):
    general_data = {key: {} for key in feeds}
    for feed in feeds:
        if verbosity > 0: print('-' * 12 + 'Dir is %s' % feed + '-' * 12)
        nNotFound = 0
        for seed in indep:
            try:
                my_dir = '%s/%s/%i' % (path, feed, seed)
                (old_begin, nfiles) = what2Analyze(my_dir, type,
                                                   guessStart, interval)
                if verbosity > 0:
                    print('old_begin = {}, nfiles = {}'.format(old_begin, nfiles))
                # get data from old files
                (N, P, boxlx, boxly, boxlz,
                 ncycle_old, molWeights, U) = reader.read_fort12(my_dir, old_begin,
                                                                 nfiles, tag=type)
                Ntotal = getTotalMols(N)

                K, mole_frac = getKX(N, Ntotal)
                (number_densities, chemical_potential,
                 swap_info, biasPot, average_volumes,
                 totalComposition, cbmc_info,
                 T, zeolite) = reader.go_through_runs(my_dir, ncycle_old,
                                                      old_begin, nfiles,
                                                      tag=type)

                # do calculations for other data that may be needed
                number_dens_real = getRealRho(number_densities, biasPot, T)
                deltaG = calcDGfromNumDens(number_dens_real, totalComposition, T)
                if 'energies' in kwargs.keys() and kwargs['energies'] == 'Yes':
                    deltaU, deltaH = calc_dU_dH(U, P, Ntotal, boxlx, boxly, boxlz)
                if liq:
                    concentrations = {}
                    assert kwargs['box'], 'box needed for liquid phase'
                    if 'box' not in kwargs['box']: kwargs['box'] = 'box' + kwargs['box']
                    l_box = kwargs['box']
                    c = calc_tools.g_mL(N[mol][l_box], boxlx[l_box],
                                        MW=molWeights[mol])
                    concentrations[mol] = {l_box: c}
                # initialize vars
                if (seed == indep[0]) and (feed == feeds[0]):
                    k_ratio = properties.AnyProperty(K)
                    X = properties.AnyProperty(mole_frac)
                    CBMC = properties.AnyProperty(cbmc_info)
                    Press = properties.AnyProperty(P)
                    Nmlcl = properties.AnyProperty(N)
                    SWAP = properties.AnyProperty(swap_info)
                    rho = properties.AnyProperty(number_dens_real)
                    dG = properties.AnyProperty(deltaG)
                    data = {'CBMC': CBMC, 'P': Press, 'N': Nmlcl,
                            'SWAP': SWAP, 'rho': rho,
                            'dG': dG, 'K': k_ratio,
                            'X': X}
                    if liq:
                        if (verbosity > 1):
                            print('Doing analysis for C [ g/mL ] for mol {}'.format(mol))
                            print('Box 2 should be liquid phase')
                        C = properties.AnyProperty(concentrations)
                        data['Conc'] = C
                    if 'energies' in kwargs.keys() and kwargs['energies'] == 'Yes':
                        dH = properties.AnyProperty(deltaH)
                        dU = properties.AnyProperty(deltaU)
                        data['dH'] = dH
                        data['dU'] = dU
                CBMC.addVals(cbmc_info)
                Press.addVals(P)
                Nmlcl.addVals(N)
                SWAP.addVals(swap_info)
                rho.addVals(number_dens_real)
                dG.addVals(deltaG)
                k_ratio.addVals(K)
                X.addVals(mole_frac)
                if liq:
                    C.addVals(concentrations)
                if 'energies' in kwargs.keys() and kwargs['energies'] == 'Yes':
                    dU.addVals(deltaU)
                    dH.addVals(deltaH)
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
