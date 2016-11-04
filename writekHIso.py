'''
Write isotherm from previously generated databank files
'''
from runAnalyzer import checkRun, calc95conf
from writeXvY import writeAGR
from chem_constants import R, N_av
import math

if __name__ == '__main__':
    from parser import Plot
    import shelve

    my_parser = Plot()
    my_parser.isotherm()
    my_parser.kH()

    args = vars(my_parser.parse_args())
    assert args['units'], 'Units must be defined for isotherm'
    assert args['box'], 'Box needed for number density in kH isotherm'
    assert args['Temp'], 'Temperature needed for ideal gas law'
    assert (args['henry'] and
            len(args['henry']) == 2), 'kH(mean) and kH(95% conf) needed'
    assert args['mol'], 'Mol needed for x axis'

    N = {}; rho = {}; gen_data = {}
    for file, var in zip(['N-data.db','rho-data.db', 'general-data.db'],
                         [N, rho, gen_data]):
        with shelve.open('%s/%s'%(args['path'], file)) as db:
            for feed in args['feeds']:
                assert feed in db.keys(), 'Feed {} not in database'.format(feed)
                var[feed] = db[feed]

    kH_mean, kH_95conf = args['henry']
    for feed in args['feeds']:
        if args['verbosity'] > 0:
            print('starting feed {}'.format(feed))
        # determine if run has been completed
        run = checkRun(args['type'], [rho, N], feed)
        # gen data
        numIndep = gen_data[feed][run]['numIndep']
        # initialize variables if needed
        mol_data = {}
        mols_adsorbed = sorted([i for i in N[feed][run].keys()
                                    if N[feed][run][i]['box1']['mean'] > 1e-06])
        # pressure info ----------------------
        rho_mean, rho_stdev = (rho[feed][run][args['mol']]['box%s'%args['box']]['mean'],
                                rho[feed][run][args['mol']]['box%s'%args['box']]['stdev'])
        for mol in mols_adsorbed:
            if mol not in mol_data.keys():
                mol_data[mol] = {}
                for key in ['loading', 'concentrations']:
                    mol_data[mol][key] = {'mean':[],'95conf':[],'feed':[]}
            loadings = mol_data[mol]['loading']
            concentrations = mol_data[mol]['concentrations']
            # loading info ----------------------
            if args['units'] == 'molec/uc':
                qfactor = 1/gen_data[feed][run]['zeolite']['unit cells']
            elif args['units'] == 'g/g':
                qfactor = (gen_data[feed][run]['molecular weight'][args['mol']]/
                                N_av)/gen_data[feed]['zeolite']['mass (g)']
            elif args['units'] == 'mol/kg':
                qfactor = gen_data[feed][run]['zeolite'][' mol/kg / 1 mlcl adsorbed']
            loadings['mean'].append( N[feed][run][mol]['box1']['mean']*qfactor )
            loadings['95conf'].append( calc95conf(
                N[feed][run][mol]['box1']['stdev']*qfactor, numIndep ) )
            loadings['feed'].append( feed )
            # pressure info ----------------------
            # convert number density to pressure [kPa]
            # (molec/nm**3)*(mol/molec)*(nm**3*kPa/(mol*K))*K = kPa
            p_mean = rho_mean/N_av*R['nm**3*kPa/(mol*K)']*args['Temp']
            p_stdev = rho_stdev/N_av*R['nm**3*kPa/(mol*K)']*args['Temp']
            p_95conf = calc95conf(p_stdev, numIndep)
            C_95conf = p_mean*kH_mean*math.pow(
                        math.pow(kH_95conf/kH_mean, 2) +
                        math.pow(p_95conf/p_mean, 2), 0.5 )
            concentrations['mean'].append( p_mean*kH_mean )
            concentrations['95conf'].append( C_95conf  )
            concentrations['feed'].append( feed )
        # write after each feed!
        for mol in mol_data.keys():
            iso_name = 'isotherm-mols-%s-mol%s-%s'%('_'.join(mols_adsorbed),mol, args['units'])
            iso_name = iso_name.replace('/','_')
            help = ''
            x_data = mol_data[mol]['concentrations']
            y_data = mol_data[mol]['loading']
            if feed == args['feeds'][0]: help = 'C(g/mL)    Q(%s)    dC     dQ'%args['units']
            writeAGR(x_data['mean'], y_data['mean'],
                     x_data['95conf'], y_data['95conf'],
                     x_data['feed'], iso_name, description=help)

