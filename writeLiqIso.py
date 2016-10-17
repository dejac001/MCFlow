'''
Write isotherm from previously generated databank files
'''
from runAnalyzer import checkRun, getConc, calc95conf
from file_formatting.writer import writeAGR

if __name__ == '__main__':
    from parser import Plot
    import shelve
    from chem_constants import N_av

    my_parser = Plot()
    my_parser.isotherm()

    args = vars(my_parser.parse_args())
    assert args['units'], 'Units must be defined for isotherm'

    N = {}; C = {}; gen_data = {}
    for file, var in zip(['N-data.db','Conc-data.db', 'general-data.db'],
                         [N, C, gen_data]):
        with shelve.open('%s/%s'%(args['path'], file)) as db:
            for feed in args['feeds']:
                assert feed in db.keys(), 'Feed {} not in database'.format(feed)
                var[feed] = db[feed]


    for feed in args['feeds']:
        # determine if run has been completed
        run = checkRun(args['type'], list(C[feed].keys()), list(N[feed].keys()))
        # gen data
        numIndep = gen_data[feed][run]['numIndep']
        # concentrations
        mol, box, conc = getConc(C[feed])
        if args['verbosity'] > 0:
            print('liquid analysis for feed %s is with mol %s, box %s'%(feed, mol, box))
        # initialize variables if needed
        mol_data = {}
        mols_adsorbed = sorted([i for i in N[feed][run].keys()
                                    if N[feed][run][i]['box1']['mean'] > 1e-06])
        for mol in mols_adsorbed:
            if mol not in mol_data.keys():
                mol_data[mol] = {}
                for key in ['loading', 'concentrations']:
                    mol_data[mol][key] = {'mean':[],'95conf':[],'feed':[]}
            loadings = mol_data[mol]['loading']
            concentrations = mol_data[mol]['concentrations']
            if args['units'] == 'molec/uc':
                qfactor = 1/gen_data[feed][run]['zeolite']['unit cells']
            elif args['units'] == 'g/g':
                qfactor = (gen_data[feed][run]['molecular weight'][args['mol']]/
                                N_av)/gen_data[feed]['zeolite']['mass (g)']
            elif args['units'] == 'mol/kg':
                qfactor = gen_data[feed][run]['zeolite'][' mol/kg / 1 mlcl adsorbed']
            loadings['mean'].append( N[feed][run][mol]['box1']['mean']*qfactor )
            loadings['95conf'].append( calc95conf( N[feed][run][mol]['box1']['stdev']*qfactor ) )
            loadings['feed'].append( feed )
            concentrations['mean'].append( conc['mean'] )
            concentrations['95conf'].append( calc95conf(conc['stdev'], numIndep)  )
            concentrations['feed'].append( feed )
        # write after each feed!
        for mol in mol_data.keys():
            iso_name = 'isotherm-mols-%s-mol%s-%s'%('_'.join(mols_adsorbed),mol, args['units'])
            iso_name = iso_name.replace('/','_')
            x_data = mol_data[mol]['concentrations']
            y_data = mol_data[mol]['loadings']
            writeAGR(x_data['mean'], y_data['mean'],
                     x_data['95conf'], y_data['95conf'],
                     x_data['feed'], iso_name)
