__author__ = 'dejacor'

def getResults(data, mol, MW):
    conv_factor = 1/N_av*10**21*MW
    P = data.P[data.feed][data.run]
    rho = data.rho[data.feed][data.run]
    nIndep = data.gen_data[data.feed][data.run]['numIndep']
    if 'high dens' in rho[mol].keys() and 'low dens' in rho[mol].keys():
        gas_box, liq_box = 'low dens', 'high dens'
        for box in 'box1','box2':
            print(P[box]['mean'],P[box]['stdev'])
        p_box = 'box2'
        print('using dpd for %s'%data.feed)
    else:
        if  rho[mol]['box2']['mean'] < rho[mol]['box1']['mean']:
            gas_box, liq_box = 'box2', 'box1'
        else:
            gas_box, liq_box = 'box1', 'box2'
        for box in gas_box, liq_box:
            rho[mol][box]['95conf'] = calc95conf(rho[mol][box]['stdev'], nIndep)
        p_box = gas_box
    print('gas box is %s and liquid box is %s'%(gas_box,liq_box))
    return (
            rho[mol][gas_box]['mean']*conv_factor,
            rho[mol][gas_box]['95conf']*conv_factor,
            rho[mol][liq_box]['mean']*conv_factor,
            rho[mol][liq_box]['95conf']*conv_factor,
        P[p_box]['mean'], calc95conf(P[p_box]['stdev'],nIndep)
            )
    # return (
    #     rho[mol][gas_box]['mean'], rho[mol][gas_box]['stdev']/np.sqrt(nIndep),
    #     rho[mol][liq_box]['mean'], rho[mol][liq_box]['stdev']/np.sqrt(nIndep),
    #     P[gas_box]['mean'], P[gas_box]['stdev']/np.sqrt(nIndep)
    #         )

import numpy as np
from chem_constants import N_av
from MCFlow.file_organization import equilName, prodName

if __name__ == '__main__':
    import argparse, os
    from writeXvY import IdealGasAds as VLCC_data
    from MCFlow.runAnalyzer import checkRun, calc95conf
    parser = argparse.ArgumentParser(description='generate VLCC input for fortran files')
    parser.add_argument('-f','--feeds',help='feeds',nargs = '+',type=str)
    parser.add_argument('-T','--temperatures',help='temperatures corresponding to feeds',
                        nargs='+',type=float)
    parser.add_argument('-t','--type',help='type of run',type=str,
                        choices=[equilName,prodName],default=prodName)
    parser.add_argument('-MW','--molWeight',help='molecular weight of molecule (g / mol)',
                        type=float,default=104.15)
    parser.add_argument('-m','--mol',help='molecule',type=str,default='1')
    args = vars(parser.parse_args())
    assert args['feeds'], 'No feeds given'
    assert (len(args['feeds']) == len(args['temperatures'])), 'Num feeds %i not equal to num temps %i'%(len(args['feeds']), len(args['temperatures']))
    assert args['temperatures'][0] == min(args['temperatures']), 'First T is not smallest'
    my_data = VLCC_data()
    my_data.files.append( 'P-data.json')
    my_data.path = os.getcwd()
    my_data.feeds = args['feeds']
    my_data.P = {}
    my_data.variables.append( my_data.P )
    my_data.read_json()

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
            f.write('%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n'%(TK, rG, drG, rL, drL, p, dp))
            cc_string += '%12f %12f %12f\n'%(TK, p, dp)
        f.write('TMIN NPOINT for VLCC plots\n')
        f.write('%f  200\n'%(temperatures[0]*0.995))
        f.write('NDATA  STP\n')
        f.write('%i      101.325\n'%(len(temperatures)))
        f.write(' T    PRESSURE    ERROR\n')
        f.write(cc_string)
        f.write('TMIN TMAX NPOINT for Clausius-Clapeyron plots\n')
        f.write('%e %e 200\n'%(0.99*temperatures[0],1.01*temperatures[-1]))
