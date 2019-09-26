'''
Write isotherm from previously generated databank files
'''
def getX(c_data, mol='', box=''):
    for key, value in sorted(c_data.items()):
        if (len(key) == 1) or (len(key) == 2):
            mol = key
        elif 'box' in key:
            box = key
        if isinstance(value, dict):
            if 'mean' in value.keys():
                return value, mol, box
            else:
                return getX(value, mol, box)

from MCFlow.runAnalyzer import checkRun, calc95conf
from MCFlow.writeXvY import writeAGR
from MCFlow.file_formatting import writer
from MCFlow.chem_constants import R, N_av
import math

if __name__ == '__main__':
    from analysis_parsers import Main
    import shelve

    my_parser = Main()
    my_parser.parser.add_argument('-TK','--Temp',help='Temperature in Kelvin',
                                  type=float)

    args = vars(my_parser.parse_args())

    assert args['Temp'], 'Temperature needed for ideal gas law'

    rho = {}; C = {}; gen_data = {}
    for file, var in zip(['rho-data.db','Conc-data.db', 'general-data.db'],
                         [rho, C, gen_data]):
        with shelve.open('%s/%s'%(args['path'], file)) as db:
            for feed in args['feeds']:
                assert feed in db.keys(), 'Feed {} not in database'.format(feed)
                var[feed] = db[feed]

    mol_data = {}
    for feed in args['feeds']:
        # determine if run has been completed
        run = checkRun(args['type'], [C,rho], feed)
        # gen data
        numIndep = gen_data[feed][run]['numIndep']
        # concentrations
        conc, c_mol, box = getX(C[feed])
        # initialize variables if needed
        if c_mol not in mol_data.keys():
            mol_data[c_mol] = {}
            for key in ['kH', 'P']:
                mol_data[c_mol][key] = {'mean':[],'95conf':[],'feed':[]}
        henry_constant = mol_data[c_mol]['kH']
        pressure = mol_data[c_mol]['P']
        rho_mean, rho_stdev = (rho[feed][run][c_mol]['box3']['mean'],
                    rho[feed][run][c_mol]['box3']['stdev'])
        # convert number density to pressure [kPa]
        # (molec/nm**3)*(mol/molec)*(nm**3*kPa/(mol*K))*K = kPa
        p_mean = rho_mean/N_av*R['nm**3*kPa/(mol*K)']*args['Temp']
        p_stdev = rho_stdev/N_av*R['nm**3*kPa/(mol*K)']*args['Temp']

        henry_constant['mean'].append( conc['mean']/p_mean )
        kH_stdev = conc['mean']/p_mean*math.pow(
            math.pow(conc['stdev']/conc['mean'], 2) +
            math.pow(p_stdev/p_mean, 2), 0.5 )
        henry_constant['95conf'].append( calc95conf(kH_stdev, numIndep) )
        pressure['mean'].append( p_mean )
        pressure['95conf'].append( calc95conf(p_stdev, numIndep) )
        if args['verbosity'] > 0:
            print('liquid analysis for feed %s is with mol %s, box %s'%(feed, c_mol, box))
        henry_constant['feed'].append( feed )

    # write out results
    for mol in mol_data.keys():
        file_name = 'kH-mol%s.dat'%(mol)
        file_name = file_name.replace('/','_')
        x_data = mol_data[mol]['kH']
        y_data = mol_data[mol]['P']
        message = 'kH (g/(mL*kPa)      P (kPa)      d(kH)    d(P)'
        writeAGR(x_data['mean'], y_data['mean'],
                 x_data['95conf'], y_data['95conf'],
                 x_data['feed'], file_name, description=message)
