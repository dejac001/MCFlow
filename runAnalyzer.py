def findFinalSimRun(path, tag):
    nfiles = 0
    if ('prod' in tag) and ('butanol-water' in path) and ('323K' in path):
        while os.path.isfile('%s/%s%i/run.%s%i'%(path, tag, nfiles+1, tag, nfiles+1))\
                and os.path.isfile('%s/%s%i/movie_Hbond_moltype.txt'%(path, tag, nfiles+1)):
            nfiles += 1
    else:
        while os.path.isfile('%s/%s%i/run.%s%i'%(path, tag, nfiles+1, tag, nfiles+1)):
            nfiles += 1
    return nfiles

def findNextRun(path, tag):
    runNum = 1
    while os.path.isfile(path +'/%s%i/run.%s%i'%(tag, runNum, tag, runNum)):
        runNum += 1
    return runNum

def findFinalSimIntvl(path, fileNumStart, numIntvl, tag='equil-'):
    '''
    nIntvl is integer interval of number of files at end of finished runs
    '''
    fileNum = fileNumStart
    while os.path.isfile('%s/%s%i/run.%s%i'%(path, tag, fileNum, tag, fileNum)):
        fileNum += 1
    if fileNum <= numIntvl:
        print('interval for analysis too long, just going to analyze all runs')
        runAnalStart = fileNumStart
    elif fileNum > numIntvl:
        runAnalStart = fileNum - numIntvl
    return runAnalStart

def what2Analyze(path, oldTagName, fileNumStart, numIntvl):
    TAG = oldTagName[:(oldTagName.find('-')+1)]
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
    return old_begin, nfiles

def getRealRho(rhoBias, bias, T):
    '''
    get real number densities: note: units in chain/nm**3
    '''
    rhoReal = {}
    for mol in rhoBias.keys():
        for box in rhoBias[mol].keys():
            if rhoBias[mol][box] > 0.:
                if mol not in rhoReal.keys(): rhoReal[mol] = {}
                rhoReal[mol][box] = rhoBias[mol][box] * math.exp(  bias[mol][box]/T )
    return rhoReal

def getConc(c_data, mol='', box=''):
    for key, value in sorted(c_data.items()):
        if len(key) == 1:
            mol = key
        elif 'box' in key:
            box = key
        if isinstance(value, dict):
            if 'mean' in value.keys():
                return mol, box, value
            else:
                return getConc(value, mol, box)

def calc95conf(stdev, numIndep):
    T_values = {'8':2.365,'4':3.182}
    assert '%i'%numIndep in T_values.keys(), 'No T-value stored for %i indep'%numIndep
    return stdev/math.pow(numIndep, 0.5)*T_values['%i'%numIndep]

def checkRun(tag, feeds1, feeds2):
    data_run = ''
    for feed in set(feeds1 + feeds2):
        if tag in feed: data_run += feed
    assert feed in feeds1, 'More than one run found to analyze'
    return data_run

def calcDGfromNumDens(rho, N_i, T):
    dG = {}
    for mlcl_name in N_i.keys():
        mol = mlcl_name.split()[0]
        if (N_i[mlcl_name] > 0):
            assert mol in rho.keys(), '{} not in rho.keys()!'.format(mlcl_name)
            dG[mol] = {}
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

def getFileData(feeds, indep, path, type, guessStart, interval,
                verbosity, liq=False, mol=['-1'], **kargs):
    general_data = {key:{} for key in feeds}
    for feed in feeds:
        if verbosity > 0: print('-'*12 + 'Dir is %s'%feed + '-'*12)
        for seed in indep:
            my_dir = '%s/%s/%i'%(path,feed,seed)
            (old_begin, nfiles) = what2Analyze(my_dir, type,
                                                       guessStart,interval)
            if verbosity > 0:
                print('old_begin = {}, nfiles = {}'.format(old_begin, nfiles))
            try:
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
            except FileNotFoundError:
                print('File not found for dir {}'.format(my_dir))
                print('Not averaging for this directory')
                continue
            # do calculations for other data that may be needed
            number_dens_real = getRealRho(number_densities, biasPot, T)
            if liq:
                concentrations = {}
                for mlcl in mol:
                    c = calc_tools.g_mL(N[mlcl]['box2'], boxLengths['box2'],
                                        MW=molWeights[mlcl])
                    concentrations[mlcl] = {'box2':c}
                deltaG = calcDGfromNumDens(number_dens_real, totalComposition, T)
            # initialize vars
            if (seed == indep[0]) and (feed == feeds[0]):
                boxlx = properties.AnyProperty(boxLengths)
                CBMC = properties.AnyProperty(cbmc_info)
                Press = properties.AnyProperty(P)
                Nmlcl = properties.AnyProperty(N)
                SWAP = properties.AnyProperty(swap_info)
                U = properties.AnyProperty(E)
                rho = properties.AnyProperty(number_dens_real)
                data = {'CBMC':CBMC, 'P':Press, 'N':Nmlcl,
                        'SWAP':SWAP, 'U':U, 'rho':rho,
                        'boxlx':boxlx}
                if liq:
                    if (verbosity > 1):
                        print('Doing analysis for C [ g/mL ] for mol {}'.format(mol))
                        print('Box 2 should be liquid phase')
                    C = properties.AnyProperty( concentrations )
                    dG = properties.AnyProperty( deltaG )
                    data['Conc'] = C
                    data['dG'] = dG
            CBMC.addVals(cbmc_info)
            Press.addVals(P)
            Nmlcl.addVals(N)
            SWAP.addVals(swap_info)
            U.addVals(E)
            rho.addVals(number_dens_real)
            boxlx.addVals(boxLengths)
            if liq:
                C.addVals(concentrations)
                dG.addVals(deltaG)
        for cls in data.values():
            cls.avgVals(feed)
        # get general data
        general_data[feed]['compositions'] = totalComposition
        general_data[feed]['ncycle'] = ncycle_old
        general_data[feed]['molecular weight'] = molWeights
        general_data[feed]['numIndep'] = len(indep)
        if zeolite:
            general_data[feed]['zeolite'] = zeolite
    return data, general_data

import math, calc_tools
import properties, os
from file_formatting import reader
from chem_constants import R
