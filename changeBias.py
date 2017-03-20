def newBias(number_densities, boxLengths, N, biasOld, T, pressure, vaporBox='box3'):
    '''
    changing the biasing potential here is based on the following:
    G(desired) = G(real) + R*E(bias)
    where G = -RT*ln(rho)
    :param number_densities: already REAL numberdensities !!!
    :param boxLengths:
    :param N:
    :param biasOld:
    :param T:
    :param pressure: pressure in bar
    :param vaporBox:
    :return:
    '''
    # get all box volumes
    # determine new vapor box size
    # calculate real dG's for impurities/water
    # calculate new biasing potentials for impurities
    # calculate new biasing potential for solvent in vapor
    def getBiasMol(molNum, nVapor, volVapor, rhoReal, vaporBox='box3'):
        molNum = molNum.strip('mol')
        volVapor = volVapor/1000 # convert to nm**3
        rho_vapor = rhoReal[molNum][vaporBox]['mean']
        rho_desired = nVapor/volVapor
        bias = -T*math.log( rho_desired/rho_vapor )
        return bias
    def getVaporVolume(sorbates, rhoReal, vaporBox='box3', liquidBox='box2', nVapor=[2]):
        '''
        calculate so that 2 sorbate molecules in vapor box
        '''
        bias_sorbates = {}
        if len(sorbates) == 1:
            mol = sorbates[0]
            volume_AA3 = nVapor[0]/rhoReal[mol.strip('mol')][vaporBox]['mean']*1000
            boxlx_AA3 = math.pow( volume_AA3, 1/3 )
            least_volatile = mol
            print('- To have {} sorbate molecules of type {} in vapor phase, boxlength should be {}'.format(
                nVapor[0], mol, boxlx_AA3))
        else:
            # we have multiple sorbates
            Kmin = 10**5
            for mol in sorbates:
                K = rhoReal[mol][vaporBox]['mean']/rhoReal[mol][liquidBox]['mean']
                if K < Kmin:
                    Kmin = K
                    least_volatile = mol
                    volume_AA3 = nVapor[sorbates.index(mol)]/rhoReal[mol][vaporBox]['mean']*1000
                    boxlx_AA3 = math.pow( volume_AA3, 1/3 )
            print('- To have {} molecules of least volatile mol#{} in vapor phase,'
                                    ' boxlength should be {}'.format(
                nVapor[sorbates.index(least_volatile)],least_volatile, boxlx_AA3))
        return boxlx_AA3, volume_AA3, least_volatile
    def getBiasImpurity(molNum, nLiquid, nVapor, volLiq, volVapor,
                        rhoReal, T, vaporBox='box3', liquidBox='box2'):
        if vaporBox != 'box3':
            print(' vapor box must be box 3, quitting')
            quit()
        molNum = molNum.strip('mol')
        num_dens_ratio_real = rhoReal[molNum][vaporBox]['mean']/rhoReal[molNum][liquidBox]['mean']
        num_dens_ratio_desired = (nVapor/volVapor)/(nLiquid/volLiq)
        bias = -T*math.log( num_dens_ratio_desired/num_dens_ratio_real )
        return bias
    def themGhosts(vapor_volume,T,P):
        '''
        assume vapor volume will take size of ghost molecules if ideal gas
        :param vapor_volume: in \AA**3
        :value 83.14 is cm**3*bar/mol/K
        :param P in bar
        :return: number of ghost molecules
        '''
        nGhost = int(   vapor_volume*10**(-24)*P/(83.14*T)*N_av    )
        return nGhost

    pressure = pressure[vapor_box]['mean']
    N = {'mol%s'%mol: {box: N[mol][box] for box in N[mol].keys()} for mol in N.keys()} # get same notation
    bias_new = {mol: {box: 0 for box in N[mol].keys()} for mol in N.keys()}
    nbox = len(N['mol1'].keys())
    rho = number_densities # already have taken out biasing potentials!!

    sorbates = []
    impurities = []
    solvent = -1
    for mol in N.keys():
        total_number_mol = sum([N[mol][i]['mean'] for i in N[mol].keys()])
        if total_number_mol > 400.:
            if solvent == -1:
                # molecule is solvent
                nchain_solvent = total_number_mol
                solvent = mol
                print('- solvent molecule is {}'.format(mol))
            else:
                # molecule is a sorbate
                nchain_sorbate = total_number_mol
                sorbates.append(mol)
                print('- sorbate molecule is {}'.format(mol))
        elif total_number_mol > 3.:
            # molecule is a sorbate
            nchain_sorbate = total_number_mol
            sorbates.append(mol)
            print('- sorbate molecule is {}'.format(mol))
        elif total_number_mol > 0.:
            # molecule is impurity
            nchain_impurity = total_number_mol
            impurities.append( mol )
            print('- impurity molecule is {}'.format(mol))

    if (N[sorbates[0]]['box1']['mean'] == 0.) and (nbox == 3):
        print(' no "sorbates" in zeolite')
        nSorbate_vapor = [0.50*nchain_sorbate]
#   elif nchain_sorbate > 250:
#       nSorbate_vapor = [0.25*nchain_sorbate]
    else:
        nSorbate_vapor = []
        for mol in sorbates:
#           print(N[mol]['box2']['mean'],N[mol]['box3']['mean'])
#           nSorbate_vapor.append( (N[mol]['box2']['mean'] + N[mol]['box3']['mean'])/2 )
            print('average number of sorbate mol#{} in liquid was {}'.format(mol, N[mol]['box2']['mean']))
            if (N[mol]['box2']['mean'] + N[mol]['box3']['mean']) < 4.:
                nSorbate_vapor.append(
                    (N[mol]['box2']['mean']
                     +N[mol]['box3']['mean'])/2) # make number in vapor the same as that in liquid
            elif ((N[mol]['box2']['mean'] > 10) and (N[mol]['box3']['mean'] > 2)
                    and (N[mol]['box3']['mean'] < 10)):
                # if more than 10 in liq phase and between 2 and 10 in vapor phase, keep old value
                print('-/-/ More than 10 molecules in liq phase and between 2 and '
                      '10 in vapor phase--keeping old value')
                nSorbate_vapor.append( N[mol]['box3']['mean'] )
            else:
                nSorbate_vapor.append( 2 )

    print('nSorbate_vapor is   ', nSorbate_vapor)
    print('vaporBox is ',vaporBox)
    boxlength_vapor_AA3, volume_vapor_AA3, least_volatile = getVaporVolume(sorbates, rho,
                                                                           vaporBox=vaporBox, nVapor=nSorbate_vapor)
    if boxlength_vapor_AA3 > 10000.:
        print(' ### boxlength too high, making 10,000 \AA and putting biasing potential on sorbate')
        boxlength_vapor_AA3 = 10000.
        volume_vapor_AA3 = 10000.**3
        least_volatile = '-1'
    nGhost = themGhosts(volume_vapor_AA3,T,pressure)
    liquid_volume_AA3 = math.pow( boxLengths['box2']['mean'], 3)
    for mol in impurities:
        if (N[mol]['box2']['mean'] > 0.) and (N[mol]['box3']['mean'] > 0.):
            nchain_liquid = 0.5*nchain_impurity
            nchain_vapor = 0.5*nchain_impurity
            bias_new[mol][vaporBox] = getBiasImpurity(mol, nchain_liquid, nchain_vapor, liquid_volume_AA3,
                                                  volume_vapor_AA3, rho, T)
        elif (N[mol]['box1']['mean'] > 0.) and (N[mol]['box3']['mean'] > 0.):
            # impurity between vapor and zeolite
            zeo_volume_AA3 = math.pow( boxLengths['box2']['mean'], 3)
            nchain_zeo = 0.5*( N[mol]['box1']['mean'] + N[mol]['box3']['mean'] )
            nchain_vapor = nchain_zeo
            bias_new[mol][vaporBox] = getBiasImpurity(mol, nchain_zeo, nchain_vapor, zeo_volume_AA3,
                                                  volume_vapor_AA3, rho, T, liquidBox='box1')
        else:
            bias_new[mol][vaporBox] = 0.
        print('Bias for mol {} (impurity) changed from {} to {}'.format(mol, biasOld[mol][vaporBox],
                                                                        bias_new[mol][vaporBox] ))
    for mol in sorbates:
        if mol != least_volatile:
            if N[mol]['box2']['mean'] == 0.:
                bias_new[mol][vaporBox] = 0.
            else:
                bias_new[mol][vaporBox] = getBiasMol(mol, nSorbate_vapor[sorbates.index(mol)],
                                                        volume_vapor_AA3, rho, vaporBox=vaporBox)
            print('Bias for mol {} (sorbate) changed from {} to {}'.format(mol,
                                                                           biasOld[mol][vaporBox],
                                                                           bias_new[mol][vaporBox] ))
    bias_new[solvent][vaporBox] = getBiasMol(solvent, 30, volume_vapor_AA3, rho, vaporBox=vaporBox)
    print(' using {} ghost molecules'.format(nGhost))
    return boxlength_vapor_AA3, bias_new, nGhost

import math
from chem_constants import N_av

if __name__ == '__main__':
    from parser import Change
    from runAnalyzer import getFileData, findNextRun
    from file_formatting import writer, reader
    from getData import outputGenDB, outputDB
    from setup_production import iaverage

    args = vars(Change().parse_args())
    feeds = args.pop('feeds')

    vapor_box = 'box3' # default vapor box to box3

    for feed in feeds:
        args['feeds'] = [feed]
        data, gen_data = getFileData(**args)
        nbox = len(data['rho'].averages[feed].keys())
        input_data = reader.read_fort4(args['path']+'/' + feed + '/1/fort.4')
        (vapor_boxlx_AA3, bias, nGhost) = newBias(data['rho'].averages[feed],
                                                            data['boxlx'].averages[feed],
                                                            data['N'].averages[feed],
                                                  input_data['UNIFORM_BIASING_POTENTIALS'],
                                                  gen_data[feed]['temperature'],
                                                  data['P'].averages[feed], vaporBox=vapor_box)
        # change data
        input_data['SIMULATION_BOX'][vapor_box]['dimensions'] = '%8.2f %8.2f %8.2f'%(vapor_boxlx_AA3,
                                                                            vapor_boxlx_AA3,vapor_boxlx_AA3)
        input_data['UNIFORM_BIASING_POTENTIALS'] = bias
        input_data['SIMULATION_BOX'][vapor_box]['rcut'] = '14.0d0'
        input_data['SIMULATION_BOX'][vapor_box]['nghost'] = '%i'%nGhost
        input_data['&mc_shared']['nstep'] = '%i'%args['nstep']
        defaults = input_data['SIMULATION_BOX'][vapor_box]['defaults']
        input_data['SIMULATION_BOX'][vapor_box]['defaults'] = defaults.replace(
                                                defaults.split()[-2], 'T' ) # make vapor box ideal
        iprint, iblock = iaverage(args['nstep'])
        input_data['&analysis']['iprint'] = '%i'%iprint
        input_data['&analysis']['iblock'] = '%i'%iblock

        # write new files
        for sim in gen_data[feed]['indepSims']:
            input_data['&mc_shared']['seed'] = '%i'%sim
            my_path = '%s/%s/%i'%(args['path'],feed,sim)
            nextRun = findNextRun(my_path, args['type'])
            new_file = '%s/fort.4.%s%i'%(my_path, args['type'], nextRun)
            writer.write_fort4(input_data, new_file)
        outputDB(args['path'], args['feeds'],args['type'], data )
        outputGenDB(args['path'], args['feeds'],args['type'], gen_data )
