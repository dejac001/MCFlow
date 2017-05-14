__author__ = 'dejacor'

def getResults(data, mol, MW):
    conv_factor = 1/N_av*10**21*MW
    P = data.P[data.feed][data.run]
    rho = data.rho[data.feed][data.run]
    nIndep = data.gen_data[data.feed][data.run]['numIndep']
    if rho[mol]['box2']['mean'] < rho[mol]['box1']['mean']:
        gas_box, liq_box = 'box2', 'box1'
    else:
        gas_box, liq_box = 'box1', 'box2'
    print('gas box is %s and liquid box is %s'%(gas_box,liq_box))
    return (
        rho[mol][gas_box]['mean'], calc95conf(rho[mol][gas_box]['stdev'],nIndep),
        rho[mol][liq_box]['mean'], calc95conf(rho[mol][liq_box]['stdev'],nIndep),
        P[gas_box]['mean'], calc95conf(P[gas_box]['stdev'],nIndep)
            )
    # return (
    #     rho[mol][gas_box]['mean'], rho[mol][gas_box]['stdev']/np.sqrt(nIndep),
    #     rho[mol][liq_box]['mean'], rho[mol][liq_box]['stdev']/np.sqrt(nIndep),
    #     P[gas_box]['mean'], P[gas_box]['stdev']/np.sqrt(nIndep)
    #         )

import numpy as np
from chem_constants import N_av

if __name__ == '__main__':
    import argparse, os
    from writeXvY import IdealGasAds as VLCC_data
    from MCFlow.runAnalyzer import checkRun, calc95conf
    parser = argparse.ArgumentParser(description='generate VLCC input for fortran files')
    parser.add_argument('-f','--feeds',help='feeds',nargs = '+',type=str)
    parser.add_argument('-T','--temperatures',help='temperatures corresponding to feeds',
                        nargs='+',type=float)
    parser.add_argument('-t','--type',help='type of run',type=str,
                        choices=['equil-','prod-'],default='prod-')
    parser.add_argument('-MW','--molWeight',help='molecular weight of molecule (g / mol)',
                        type=float,default=104.15)
    parser.add_argument('-m','--mol',help='molecule',type=str,default='1')
    args = vars(parser.parse_args())
    assert args['feeds'], 'No feeds given'
    assert (len(args['feeds']) == len(args['temperatures'])), 'Num feeds not equal to num temps'
    assert args['temperatures'][0] == min(args['temperatures']), 'First T is not smallest'
    my_data = VLCC_data()
    my_data.files.append( 'P-data.db')
    my_data.path = os.getcwd()
    my_data.feeds = args['feeds']
    my_data.P = {}
    my_data.variables.append( my_data.P )
    my_data.readDBs()

    # gather data
    temperatures = []; gas_densities = {'mean':[],'error':[]}; liquid_densities = {'mean':[],'error':[]}
    pressures = {'mean':[],'error':[]}
    for feed, T in zip(args['feeds'], args['temperatures']):
        my_data.run = checkRun(args['type'],my_data.variables,feed)
        my_data.feed = feed
        rhoG, drhoG, rhoL, drhoL, p, dp = getResults(my_data, args['mol'],
                                                     args['molWeight'])
        temperatures.append( T )
        gas_densities['mean'].append(rhoG)
        gas_densities['error'].append(drhoG)
        liquid_densities['mean'].append(rhoL)
        liquid_densities['error'].append(drhoL)
        pressures['mean'].append(p)
        pressures['error'].append(dp)

    # write to file
    with open('vlcc.inp','w') as f:
        f.write(str(args['feeds'])+str(args['temperatures']) + '\n')
        f.write('NDATA  BETA    NITER\n')
        f.write('%i      0.326   100\n'%(len(temperatures)))
        f.write(' T      RHOG      ERROR      RHOL      ERROR    '
                '  PRESSURE      ERROR (standard deviation or standard error)\n')
        cc_string = ''
        for TK, rG, drG, rL, drL, p, dp in zip(temperatures,gas_densities['mean'],gas_densities['error'],
                                               liquid_densities['mean'],liquid_densities['error'],
                                               pressures['mean'],pressures['error']):
            f.write('%e %e %e %e %e %e %e\n'%(TK, rG, drG, rL, drL, p, dp))
            cc_string += '%e %e %e\n'%(TK, p, dp)
        f.write('TMIN NPOINT for VLCC plots\n')
        f.write('%e  200\n'%(temperatures[0]*0.995))
        f.write('NDATA  STP\n')
        f.write('%i      101.325\n'%(len(temperatures)))
        f.write(' T    PRESSURE    ERROR\n')
        f.write(cc_string)
        f.write('TMIN TMAX NPOINT for Clausius-Clapeyron plots\n')
        f.write('%e %e 200\n'%(0.99*temperatures[0],1.01*temperatures[-1]))