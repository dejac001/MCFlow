class NoSwapsAccepted(BaseException):
    pass


def calculateProbs(nMol, nBeadMolty):
    '''
    currently only works for binary butanol/water. Assumes that molecules are not perfectly
    linear (like CO2) and have 3 translational and 3 rotation degrees of freedom

    Does not work for one atom molecules which have no rotational DOF

    note that translational_fraction returned is cbmc_fraction + translation_fraction

    :param nMol: number of molecules of each type in tuple or list
    :param nBeadMolty: list or tuple. number of cbmc DOF for each molecule type.
    :return: ratio of probabilities for cbmc and translation, and ratio of type of cbmc for each
    '''
    if len(nBeadMolty) != len(nMol):
        print('----------number of beads for each molecule not provided correctly')
        print('----------',nMol, nBeadMolty)
        raise TypeError
    dof_cbmc_molty = []
    dof_trans_molty = []
    dof_rot_molty = []
    for molNum in range(len(nMol)):
        if nBeadMolty[molNum] == 0:
            # rigid molecule
            dof_cbmc_molty.append(0)
            dof_trans_molty.append(3*nMol[molNum])
            dof_rot_molty.append(3*nMol[molNum])
        else:
            dof_cbmc_molty.append(nBeadMolty[molNum]*nMol[molNum])
            dof_trans_molty.append(3*nMol[molNum])
            dof_rot_molty.append(3*nMol[molNum])
    translation_fraction = sum(dof_trans_molty)/ (sum(dof_cbmc_molty) + sum(dof_rot_molty)
                                                  + sum(dof_trans_molty))
    rotation_fraction = sum(dof_rot_molty)/ (sum(dof_cbmc_molty) + sum(dof_rot_molty)
                                                  + sum(dof_trans_molty))
    cbmc_fraction = sum(dof_cbmc_molty)/ (sum(dof_trans_molty) + sum(dof_rot_molty)+
                                          sum(dof_cbmc_molty))
    # now format for fort.4 file
    translation_fraction = cbmc_fraction + translation_fraction
    return cbmc_fraction, translation_fraction

def calcCBMCmolTy(nMlcl, nBeads):
    probabilityCBMCmolTy = []
    totalProb = 0.
    for mol in range(len(nMlcl)):
        if sum([i*j for i,j in zip(nMlcl, nBeads)]) > 0.:
            cbmcFrac = nMlcl[mol]*nBeads[mol]/sum([i*j for i,j in zip(nMlcl, nBeads)])
        else:
            cbmcFrac = 0
        totalProb += cbmcFrac
        probabilityCBMCmolTy.append( totalProb )
    return probabilityCBMCmolTy


def analyzeTransfers(my_transferInfo, ncycle, numberMoleculeTypes,
                         nActPerCycle=1, tavol = 0.4, maxProb=0.6):
    '''
    Works for swaps and swatches of arbitrary numbers of molecules.
    Does NOT account for number of molecules.
    Works for any amount of boxes.

    :param my_transferInfo: swap info in dictionary form read from run files. Includes swatch info
    :param nindep: number of independent sims
    :param ncycle: number of total cycles
    :param numberMoleculeTypes: number of molecules of each type
    :param nActPerCycle: number of swap moves in total that you want to accept per cycle
    :param tavol:
    :param maxProb: maximum probability of swap
    :return:
    '''
    numMolType = numberMoleculeTypes.copy() # don't want to change in main() too
    transferInfo = my_transferInfo.copy()
    nchain = sum([ i for i in numMolType.values()])

    swapInfo = {}
    swatchInfo = {}
    pmvol = nActPerCycle / (nchain*tavol) # desired pmvol
    newSwaps = {}
    newSwatches = {}
    pctAct = {}
    nActCycle = 0
    numberSwatchtype = 0
    numberSwapType = 0

    def addToMatrix(myMatrix, myIndex, acceptance):
        if myIndex == 0:
            myMatrix[0,0] = acceptance
        elif myIndex != len(myMatrix[:,0]):
            myMatrix[myIndex,myIndex] = acceptance
            myMatrix[myIndex-1,myIndex] = -1*acceptance
        else:
            myMatrix[myIndex-1,myIndex] = -1*acceptance

    # get acceptance info and distinguish between swaps and swatches
    moves = sorted(transferInfo.keys())
    for moveType in moves:
        if ' and ' in moveType: numberSwatchtype += 1
        pctAct[moveType]={}
        boxes = sorted(transferInfo[moveType].keys())
        for boxpair in boxes:
            accepted = transferInfo[moveType][boxpair]['accepted']['mean']
            attempted = transferInfo[moveType][boxpair]['attempted']['mean']
            nActCycle += accepted/ncycle
            if attempted > 0.:
                pctAct[moveType][boxpair] = accepted/attempted
            else:
                print('no moves attempted for move: {} b/t boxes {}'.format(moveType, boxpair))
            if (accepted == 0.):
                if ' and ' in moveType:
                    molecules = moveType.split(' and ')
                    if (numMolType[molecules[0]] > 0) and (numMolType[molecules[1]] > 0):
                        print('no swatch moves accepted for type {} between {}'.format(moveType,boxpair))
                        print('     - You need to take this move out of your fort.4 file!')
                        transferInfo[moveType].pop( boxpair , None)
                        pctAct[moveType].pop( boxpair , None)
#                       quit()
                    else:
                        # we aren't trying to do this swatch in the fort.4
                        numberSwatchtype -= 1
                elif (numMolType[moveType] != 0) and (attempted > 0):
                    print('no swap moves accepted for type {} between {}; path is {}'.format(moveType,boxpair,
                                                                                                os.getcwd()))
                    raise NoSwapsAccepted
#                   quit()
            elif boxpair[0] == boxpair[1]:
                print('Take out same box swatch from fort.4 file for {} between {}'.format(moveType,boxpair))
                transferInfo[moveType].pop( boxpair , None)
                pctAct[moveType].pop( boxpair , None)
        if ' and ' in moveType:
            swatchInfo[moveType] = transferInfo.pop(moveType)
        elif numMolType[moveType] != 0:
            swapInfo[moveType] = transferInfo.pop(moveType)
            if accepted > 0.: numberSwapType +=1
            # ^ keep all moveTypes for swatch, but not for swap
    print('number of accepted transfer moves was {}'.format(nActCycle))
    if numberSwatchtype == 0 and swatchInfo.keys():
        print('-we accidentally did a swatch')
        swatchInfo = {}

    # determine if there are extraneous moves which can be taken out
    # (1)
    # if the molecule is not an impurity, and you are using two different types of swatches to
    # transfer it in the same direction, decide which direction is best
#    original_swatches = list(swatchInfo.keys())
#    for swatch1 in original_swatches:
#        mol1, mol2 = swatch1.split(' and ')
#        if swatch1 not in swatchInfo.keys(): continue
#        for swatch2 in [i for i in original_swatches if i != swatch1]:
#            if swatch2 not in swatchInfo.keys(): continue
#            if ((mol1 in swatch2) and (len(pctAct[mol1].keys()) > 0) or
#                (mol2 in swatch2) and (len(pctAct[mol2].keys()) > 0)):
#                dirns = [i for i in swatchInfo[swatch1].keys() if i in swatchInfo[swatch2].keys()]
#                for dirn in dirns:
#                    if pctAct[swatch1][dirn] >  pctAct[swatch2][dirn]:
#                        # remove swatch2
#                        swatchInfo[swatch2].pop( dirn , None)
#                        pctAct[swatch2].pop( dirn , None)
#                        print('     - remove swatches for {} between {} in fort.4'.format(swatch2,dirn))
#                        if len(swatchInfo[swatch2].keys())  == 0:
#                            swatchInfo.pop( swatch2 , None)
#                            pctAct.pop( swatch2, None)
#                            numberSwatchtype -= 1
#                    else:
#                        # remove swatch1
#                        swatchInfo[swatch1].pop( dirn , None)
#                        pctAct[swatch1].pop( dirn , None)
#                        print('     - remove swatches for {} between {} in fort.4'.format(swatch1,dirn))
#                        if len(swatchInfo[swatch1].keys())  == 0:
#                            swatchInfo.pop( swatch1 , None)
#                            pctAct.pop( swatch1, None)
#                            numberSwatchtype -= 1
#    # (2)
#    # for a given swatch type, if doing swaps for both molecules in same direction, determine
#    # whether you should keep the swatch move and replace the swap move(s) and vice versa
#    original_swatches = list(swatchInfo.keys())
#    for swatchType in original_swatches:
#        mol1, mol2 = swatchType.split(' and ')
#        for dirn in pctAct[mol1].keys(): # for directions in swaps for mol1
#            if (dirn in swatchInfo[swatchType].keys()) and (dirn in pctAct[mol2].keys()):
#                # doing swaps for both molecules in swatch move in same direction
#                if ((pctAct[swatchType][dirn] < pctAct[mol1][dirn]) and
#                        (pctAct[swatchType][dirn] < pctAct[mol2][dirn])):
#                    # remove the swatch move
#                    swatchInfo[swatchType].pop( dirn , None)
#                    pctAct[swatchType].pop( dirn , None)
#                    print('     - remove swatches for {} between {} in fort.4'.format(swatchType,dirn))
#                    if len(swatchInfo[swatchType].keys())  == 0:
#                        swatchInfo.pop( swatchType , None)
#                        pctAct.pop( swatchType, None)
#                        numberSwatchtype -= 1
#                else:
#                    if pctAct[swatchType][dirn] > pctAct[mol1][dirn]:
#                        print('doing {} instead of {} between {}'.format(swatchType, mol1, dirn))
#                        pctAct[mol1] = {}
#                        if swapInfo[mol1][dirn]['accepted']['mean'] > 0.:
#                            swapInfo[mol1][dirn]['accepted']['mean'] = 0
#                            numberSwapType -= 1
#                    if pctAct[swatchType][dirn] > pctAct[mol2][dirn]:
#                        print('doing {} instead of {} between {}'.format(swatchType, mol2, dirn))
#                        pctAct[mol2] = {}
#                        if swapInfo[mol2][dirn]['accepted']['mean'] > 0.:
#                            swapInfo[mol2][dirn]['accepted']['mean'] = 0
#                            numberSwapType -= 1
#
    # accept the same number for each swap type
    swapMatrix = np.zeros((numberSwapType, numberSwapType))
    swap_b = np.array([[0] for i in range(numberSwapType)])
    swap_b[-1] = 1
    index = -1
    for moveType in sorted(swapInfo.keys()):
        newSwaps[moveType] = {}
        if len(swapInfo[moveType].keys()) == 2:
            # molecule is swapping in two directions
            dir1, dir2 = swapInfo[moveType].keys()
            if (swapInfo[moveType][dir1]['attempted']['mean'] > 0) and (swapInfo[moveType][dir2]['attempted']['mean'] > 0):
                A = np.matrix([[pctAct[moveType][dir1], -1*pctAct[moveType][dir2]],
                               [1                      , 1                       ]])
                b = np.matrix([[0],
                               [1]])
                pDir1, pDir2 = np.linalg.solve(A, b)
                newSwaps[moveType][dir1] = pDir1[0]
                newSwaps[moveType][dir2] = pDir2[0]
                effective_acceptance = pDir1[0]*pctAct[moveType][dir1]
                index += 1
                addToMatrix(swapMatrix, index, effective_acceptance)
        elif len(swapInfo[moveType].keys()) == 1:
            direction = list(swapInfo[moveType].keys())[0]
            if (swapInfo[moveType][direction]['attempted']['mean'] > 0) and (direction in pctAct[moveType].keys()):
                newSwaps[moveType][direction] = 1.0
                index += 1
                addToMatrix(swapMatrix, index, pctAct[moveType][direction])
        elif len(swapInfo[moveType].keys()) == 3:
            dir1, dir2, dir3 = swapInfo[moveType].keys()
            A = np.matrix([[pctAct[moveType][dir1], -1*pctAct[moveType][dir2], 0.],
                            [0.,    pctAct[moveType][dir2], -1*pctAct[moveType][dir3]],
                            [1.,                1.,             1.]])
            b = np.matrix([[0],[0],[1]])
            pDir1, pDir2, pDir3 = np.linalg.solve(A, b)
            newSwaps[moveType][dir1] = pDir1[0]
            newSwaps[moveType][dir2] = pDir2[0]
            newSwaps[moveType][dir3] = pDir3[0]
            effective_acceptance = pDir1[0]*pctAct[moveType][dir1]
            index += 1
            addToMatrix(swapMatrix, index, effective_acceptance)
        elif len(swapInfo[moveType].keys()) == 4:
            dir1, dir2, dir3, dir4 = swapInfo[moveType].keys()
            A = np.matrix([[pctAct[moveType][dir1], -1*pctAct[moveType][dir2], 0., 0.],
                            [0.,    pctAct[moveType][dir2], -1*pctAct[moveType][dir3], 0.],
                            [0.,    0., pctAct[moveType][dir3], -1*pctAct[moveType][dir4]],
                            [1.,                1.,             1., 1.]])
            b = np.matrix([[0],[0],[0],[1]])
            pDir1, pDir2, pDir3, pDir4 = np.linalg.solve(A, b)
            newSwaps[moveType][dir1] = pDir1[0]
            newSwaps[moveType][dir2] = pDir2[0]
            newSwaps[moveType][dir3] = pDir3[0]
            newSwaps[moveType][dir4] = pDir4[0]
            effective_acceptance = pDir1[0]*pctAct[moveType][dir1]
            index += 1
            addToMatrix(swapMatrix, index, effective_acceptance)
        else:
            print('!!!probability script not ready for %i different types of directions'%(len(swapInfo[moveType].keys())))
            quit()
    swapMatrix[-1,:] = 1

    # accept same number for each swatch type
    if swatchInfo.keys():
        swatchMatrix = np.zeros((numberSwatchtype, numberSwatchtype))
        swatch_b = np.array([[0] for i in range(numberSwatchtype)])
        swatch_b[-1] = 1
        index = -1
        for moveType in sorted(swatchInfo.keys()):
            molecules = moveType.split(' and ')
            if (numMolType[molecules[0]] > 0) and (numMolType[molecules[1]] > 0):
                index += 1
                newSwatches[moveType] = {}
                if len(swatchInfo[moveType].keys()) == 2:
                    dir1, dir2 = swatchInfo[moveType].keys()
                    A = np.matrix([[pctAct[moveType][dir1], -1*pctAct[moveType][dir2]],
                                 [1                      , 1                       ]])
                    b = np.matrix([[0],
                                   [1]])
                    pDir1, pDir2 = np.linalg.solve(A, b)
                    newSwatches[moveType][dir1] = pDir1[0]
                    newSwatches[moveType][dir2] = pDir2[0]
                    effective_acceptance = pDir1[0]*pctAct[moveType][dir1]
                    addToMatrix(swatchMatrix, index, effective_acceptance)
                elif len(swatchInfo[moveType].keys()) == 1:
                    direction = list(swatchInfo[moveType].keys())[0]
                    newSwatches[moveType][direction] = 1.0
                    addToMatrix(swatchMatrix, index, pctAct[moveType][direction])
                elif len(swatchInfo[moveType].keys()) == 3:
                    dir1, dir2, dir3 = swatchInfo[moveType].keys()
                    A = np.matrix([[pctAct[moveType][dir1], -1*pctAct[moveType][dir2], 0.],
                                    [0.,    pctAct[moveType][dir2], -1*pctAct[moveType][dir3]],
                                    [1.,                1.,             1.]])
                    b = np.matrix([[0],[0],[1]])
                    pDir1, pDir2, pDir3 = np.linalg.solve(A, b)
                    newSwatches[moveType][dir1] = pDir1[0]
                    newSwatches[moveType][dir2] = pDir2[0]
                    newSwatches[moveType][dir3] = pDir3[0]
                    effective_acceptance = pDir1[0]*pctAct[moveType][dir1]
                    addToMatrix(swatchMatrix, index, effective_acceptance)
                elif len(swatchInfo[moveType].keys()) == 4:
                    dir1, dir2, dir3, dir4 = swatchInfo[moveType].keys()
                    A = np.matrix([[pctAct[moveType][dir1], -1*pctAct[moveType][dir2], 0., 0.],
                                    [0.,    pctAct[moveType][dir2], -1*pctAct[moveType][dir3], 0.],
                                    [0.,    0., pctAct[moveType][dir3], -1*pctAct[moveType][dir4]],
                                    [1.,                1.,             1., 1.]])
                    b = np.matrix([[0],[0],[0],[1]])
                    pDir1, pDir2, pDir3, pDir4 = np.linalg.solve(A, b)
                    newSwatches[moveType][dir1] = pDir1[0]
                    newSwatches[moveType][dir2] = pDir2[0]
                    newSwatches[moveType][dir3] = pDir3[0]
                    newSwatches[moveType][dir4] = pDir4[0]
                    effective_acceptance = pDir1[0]*pctAct[moveType][dir1]
                    addToMatrix(swatchMatrix, index, effective_acceptance)
                else:
                    print('!!!probability script not ready for 3 different types of directions')
                    quit()
        swatchMatrix[-1,:] = 1

    # solve
    swap_probabilities = np.linalg.solve(swapMatrix, swap_b)
    pSwap = [i[0] for i in swap_probabilities]
    if swatchInfo.keys():
        swatch_probabilities = np.linalg.solve(swatchMatrix, swatch_b)
        pSwatch = [i[0] for i in swatch_probabilities]
    # add back into data
    molNum = -1
    for moveType in sorted(numMolType.keys()): # sorted here is key!
        if (numMolType[moveType] > 0) and (moveType in swapInfo.keys()):
            for dirn in swapInfo[moveType].keys():
                accepted = swapInfo[moveType][dirn]['accepted']['mean']
                if (accepted > 0) and ('total' not in newSwaps[moveType].keys()):
                    molNum += 1
                    newSwaps[moveType]['total'] = pSwap[molNum]
                elif (accepted == 0.):
                    newSwaps[moveType][dirn] = 0.
                    newSwaps[moveType]['total'] = 0.
                    pctAct[moveType][dirn] = 0.
        else:
            if moveType not in newSwaps.keys(): newSwaps[moveType] = {}
            newSwaps[moveType]['total'] = 0.
            if moveType in swapInfo.keys():
                for dirn in swapInfo[moveType].keys():
                    pctAct[moveType][dirn] = 0.
    molNum = -1
    for swatchType in sorted(swatchInfo.keys()):
        molecules = swatchType.split(' and ')
        if (numMolType[molecules[0]] > 0) and (numMolType[molecules[1]] > 0):
            molNum += 1
            newSwatches[swatchType]['total'] = pSwatch[molNum]
        else:
            newSwatches[swatchType]['total'] = 0.
    # make lists for moves by molecule type for fort.4
    orderedMols = sortMolKeys(newSwaps)
    print('ordered mols for calculating swaps is {}'.format(orderedMols))
    normSwaps = []
    totalProb = 0
    for molNum in orderedMols:
        totalProb += newSwaps[molNum]['total']
        if newSwaps[molNum]['total'] == 0.:
            normSwaps.append( 0. )
        else:
            normSwaps.append(totalProb)
    if swatchInfo.keys():
        orderedMoves = sortMolKeys(newSwatches)
        print('ordered mols for calculating swatches is {}'.format(orderedMoves))
        normSwatches = []
        totalProb = 0
        for moveType in orderedMoves:
            totalProb += newSwatches[moveType]['total']
            dirn = list(swatchInfo[moveType].keys())[0]
            swatchNum = int(swatchInfo[moveType][dirn]['swatchNum']['mean'])
            if swatchNum > len(normSwatches):
                for j in range(swatchNum - (len(normSwatches)+1)):
                    normSwatches.append( 0. )
                normSwatches.append(totalProb)
            else:
                normSwatches.append(totalProb)

        # determine alpha
        for mol in swapInfo.keys():
#           if len(swapInfo[mol].keys()) == 1: # if only 1 direction
            dirn = list(swapInfo[mol].keys())[0]
            if swapInfo[mol][dirn]['accepted']['mean'] > 0.:
                swapDir = dirn
                swapMol = mol
        for mol in swatchInfo.keys():
            if ((numMolType[molecules[0]] > 0) and (numMolType[molecules[1]] > 0)):
#               and (len(swatchInfo[mol].keys())==1)):
                swatchDir = list(swatchInfo[mol].keys())[0]
                swatchMol = mol

        # alpha: pswatch/pswap
        alpha = newSwaps[swapMol]['total']*pctAct[swapMol][swapDir]/(newSwatches[swatchMol]['total']
                                                                     *pctAct[swatchMol][swatchDir])
    else:
        alpha = 0.

    # now calculate pswaptot
    # nchain*num_cycles*pmswap*\sum(pmswmt_normalized*pctActMt) = number_of_accepted_moves
    # pmswap = nActPerCycle/(nchain*\sum(pmswmt_normalized*pctActMt))
    denominator = 0
    for mol in swapInfo.keys():
        if len(swapInfo[mol].keys()) > 1:
            for dirn in swapInfo[mol].keys():
                denominator += newSwaps[mol]['total']*newSwaps[mol][dirn]*pctAct[mol][dirn]
        else:
            dirn = list(swapInfo[mol].keys())[0]
            denominator += newSwaps[mol]['total']*pctAct[mol][dirn]
    denominator *= nchain

    denominatorSwatch = 0
    for mol in swatchInfo.keys():
        if len(swatchInfo[mol].keys()) > 1:
            for dirn in swatchInfo[mol].keys():
                denominatorSwatch += newSwatches[mol]['total']*newSwatches[mol][dirn]*pctAct[mol][dirn]
        else:
            dirn = list(swatchInfo[mol].keys())[0]
            denominatorSwatch += newSwatches[mol]['total']*pctAct[mol][dirn]
    denominatorSwatch *= nchain

    pSwapTot = nActPerCycle/ (denominator + alpha*denominatorSwatch)


    if swatchInfo.keys():
        if (pmvol + pSwapTot*(alpha + 1)) > maxProb:
            print('predicted probability too high, setting maxProb = %.4f'%maxProb)
            pSwapTot = maxProb / (    pmvol  +  alpha  + 1    )
            actualAcceptedPerCycle = pSwapTot * (denominator + alpha*denominatorSwatch)
            pmvol = actualAcceptedPerCycle/(nchain*tavol)
            print('actual accepted swaps and swatches per cycle will be %6.4f'%actualAcceptedPerCycle)
    else: #not doing swatches
        normSwatches = []
        if (pmvol + pSwapTot) > maxProb:
            print('predicted probability too high, setting maxProb = %.4f'%maxProb)
            pSwapTot = maxProb - pmvol
            actualAcceptedPerCycle = pSwapTot * (denominator)
            pmvol = actualAcceptedPerCycle/(nchain*tavol)
            print('actual accepted swaps per cycle will be %6.4f (we arent doing swatches)'%actualAcceptedPerCycle)
    if pmvol < 0.1/(nchain*tavol):
        # if we accept less than 1 volume move per 10 cycles, do 1 volume move per 10 cycles
        pmvol  = 0.1/(nchain*tavol)
    return (newSwaps, newSwatches, pctAct, nActCycle, normSwaps, normSwatches,
            pmvol, pmvol+pSwapTot*alpha, pmvol+pSwapTot*alpha+pSwapTot)

from MCFlow.runAnalyzer import getFileData, findNextRun
from MCFlow.file_formatting import writer, reader
from MCFlow.getData import outputDB, outputGenDB
import numpy as np
import os, shutil
from MCFlow.dataUtil import sortMolKeys

if __name__ == '__main__':
    from analysis_parsers import Change

    args = vars(Change().parse_args())
    feeds = args.pop('feeds')

    for feed in feeds:
        args['feeds'] = [feed]
        data, gen_data = getFileData(**args)
        outputDB(args['path'], args['feeds'],args['type'], data )
        outputGenDB(args['path'], args['feeds'],args['type'], gen_data )
        nbox = len(data['rho'].averages[feed].keys())
        # TODO: format return from analyzeTransfers to fit well with fort4 data dict
        try:
            (newSwaps, newSwatches, pctAct,
            nActCycle, normSwaps, normSwatches,
            pmvol, pswatch_norm, pswap_norm) = analyzeTransfers(data['SWAP'].averages[feed],
                                                            gen_data[feed]['ncycle'],
#                                                       nActPerCycle=3656/200, tavol=0.3,
                                                            gen_data[feed]['compositions'])
        except NoSwapsAccepted:
            print('no swaps accepted for {}'.format(feed))
            continue
        print(newSwaps)
        # read and write
        for sim in args['indep']:
            my_path = '%s/%s/%i/'%(args['path'],feed,sim)
            try:
                input_data = reader.read_fort4(my_path + 'fort.4')
            except FileNotFoundError:
                input_data = reader.read_fort4(my_path + 'old-fort4')
            # add in swaps
            if input_data['&mc_swap']['pmswap'][-2:] == 'd0':
                input_data['&mc_swap']['pmswap'] = input_data['&mc_swap']['pmswap'][:-2]
            pswap_old = float(input_data['&mc_swap']['pmswap'])
            input_data['&mc_swap']['pmswap'] = '%e'%pswap_norm
            input_data['&mc_shared']['seed'] = '%i'%sim
            nmolty = int(input_data['&mc_shared']['nmolty'])
            assert len(newSwaps) == nmolty, 'Incorrect swaps: too many'
            total_pswap = 0.
            for molNum in range(1,nmolty+1):
                for mol in newSwaps.keys():
                    if mol.split()[0] == '%i'%molNum:
                        key = mol
                        val = newSwaps[key]
                mol = 'mol%s'%(key.split()[0])
                # total value for molecule type
                total = val['total']
                if (total < 1e-12):
                    my_tot = '-1.0'
                else:
                    total_pswap = total_pswap + total
                    my_tot = '%4e'%total_pswap
                input_data['&mc_swap']['pmswmt'][mol] = my_tot
                # add in swap directions
                dirs = [i for i in val.keys() if i != 'total']
                input_data['MC_SWAP'][mol]['nswapb'] = len(dirs)
                input_data['MC_SWAP'][mol]['pmswapb'] = []
                input_data['MC_SWAP'][mol]['box1 box2'] = []
                my_sum = 0.
                for drn in dirs:
                    my_sum = my_sum + val[drn] # numpy object, can't add to itself
                    input_data['MC_SWAP'][mol]['pmswapb'].append( my_sum )
                    input_data['MC_SWAP'][mol]['box1 box2'].append( list(map(int,drn)) )
            # take out extraneous swatches
            unique_swatches = []
            for iswatch, info in input_data['MC_SWATCH'].items():
                in_line_1 = []; in_line_2 = []; splist = []
                for i, line in enumerate(info.split('\n')):
                    if not line.startswith('!'):
                        if len(in_line_1) == 0:
                            in_line_1 = list(line.split())
                        elif len(in_line_2) == 0:
                            in_line_2 = list(line.split())
                        else:
                            splist.append( list(map(int,line.split())) )
                            if len(splist) == int(in_line_1[2]): break
                swatch_label = ('! moltyp1<->moltyp2 nsampos 2xncut\n' + ' '.join(in_line_1) + '\n'
                                + '! gswatc 2x(ifrom, iprev)\n' + ' '.join( in_line_2) + '\n' +
                                '! splist\n')
                for sameBead in splist:
                    swatch_label += '%i %i\n'%(sameBead[0], sameBead[1])
                if swatch_label not in unique_swatches:
                    unique_swatches.append(swatch_label)
            # add in swatches
            input_data['MC_SWATCH'] = {}
            input_data['&mc_swatch'] = {'pmsatc':{},'nswaty':'0','pmswat':'0'}
            if pswatch_norm - pmvol < 1e-8:
                input_data['&mc_swatch']['pmswat'] = '-1.0'
            else:
                input_data['&mc_swatch']['pmswat'] = '%e'%pswatch_norm
            nSwatch = 0
            pmswatch_total = 0.
            for key, value in newSwatches.items():
                mol_names = key.split(' and ')
                mol1_num, mol2_num = [i.split()[0] for i in mol_names]
                for swatch in unique_swatches:
                    mol1, mol2 = swatch.split('\n')[1].split()[:2]
                    if mol1 == mol1_num and mol2 == mol2_num:
                        nSwatch += 1
                        s_num = 'swatch%i'%nSwatch
                        info = swatch + '! nswtcb pmswtcb\n'
                        # add in direstions
                        dirs = [i for i in value.keys() if i != 'total']
                        info += '%i'%len(dirs)
                        my_sum = 0.
                        for drn in sorted(dirs):
                            my_sum += value[drn]
                            info += ' %e'%(my_sum)
                        info += '\n! box numbers\n'
                        for drn in sorted(dirs):
                            info += drn[0] + ' ' + drn[1] + '\n'
                        input_data['MC_SWATCH'][s_num] = info
                        pmswatch_total += value['total']
                        input_data['&mc_swatch']['pmsatc'][s_num] = '%e'%pmswatch_total
            input_data['&mc_swatch']['nswaty'] = '%i'%nSwatch
            # add in other stuff
            #   general info
            input_data['&mc_shared']['time_limit'] = '%i'%args['time']
            input_data['&mc_shared']['nstep'] = '%i'%args['nstep']
            input_data['&mc_shared']['iratio'] = '500'
            input_data['&mc_shared']['rmin'] = '1.0'
            iblock = int(args['nstep']/10)
            if iblock > 1000: iblock = 1000
            input_data['&analysis']['iblock'] = '%i'%iblock
            input_data['&mc_volume']['pmvol'] = '%e'%pmvol
            #   other probabilities
            if input_data['&mc_cbmc']['pmcb'][-2:] == 'd0':
                input_data['&mc_cbmc']['pmcb'] = input_data['&mc_cbmc']['pmcb'][:-2]
            pmcb_old = float(input_data['&mc_cbmc']['pmcb']) - pswap_old
            pswap_new = pswap_norm
            if pmcb_old <= 0.:
                print('pmcb was zero previously, keeping it that way')
                pmcb_old = 0.
                pmtra_old = (1-pswap_old)/2.
                pmcb = -1.0
                pmtra = (1-pswap_new)*pmtra_old/(1-pswap_old) + pswap_new
            else:
                pmcb = (1-pswap_new)*pmcb_old/(1-pswap_old) + pswap_new
                pmtra_old = (1-pswap_old - pmcb_old)/2
                pmtra = (1-pswap_new)*pmtra_old/(1-pswap_old) +  pmcb
            assert pmtra < 1., 'Pmtra %e too large!'%pmtra
            input_data['&mc_cbmc']['pmcb'] = '%e'%pmcb
            input_data['&mc_simple']['pmtra'] = '%e'%pmtra
            #   box lengths
            for box in range(1,int(input_data['&mc_shared']['nbox'])+1):
                boxLengths = data['boxlx'].averages[feed]['box%i'%box]['mean']
                input_data['SIMULATION_BOX']['box%i'%box]['dimensions'] = (
                                    '%8f %8f %8f'%(boxLengths, boxLengths, boxLengths)
                                                                        )
                if (boxLengths > 100.) and (args['rcut'] < 1.):
                    input_data['SIMULATION_BOX']['box%i'%box]['rcut'] = '%5e'%(boxLengths*args['rcut'])
                if (input_data['SIMULATION_BOX']['box%i'%box]['defaults'].split()[-2] == 'T'):
                    # if ideal gas
                    input_data['SIMULATION_BOX']['box%i'%box]['rcut'] = '14.0'
            # make new run
            old_fort4 = [i for i in os.listdir(my_path) if 'fort.4.%s'%args['type'] in i]
            for fort4 in old_fort4:
                os.remove(my_path + fort4)
            try:
                shutil.move(my_path + 'fort.4',my_path + 'old-fort4')
            except FileNotFoundError:
                pass
            nextRun = 'fort.4.%s%i'%(args['type'],findNextRun(my_path,args['type']))
            writer.write_fort4(input_data, my_path + nextRun)
            input_data = None
