
from runAnalyzer import getConc, checkRun, calc95conf

if __name__ == '__main__':
    from file_formatting.writer import writeAGR
    from parser import Main
    import shelve

    my_parser = Main()

    args = vars(my_parser.parse_args())

    dG = {}; C = {}; N = {}; gen_data = {}
    for file, var in zip(['dG-data.db','Conc-data.db', 'N-data.db','general-data.db'],
                         [dG, C, N, gen_data]):
        with shelve.open('%s/%s'%(args['path'], file)) as db:
            for feed in args['feeds']:
                assert feed in db.keys(), 'Feed {} not in database'.format(feed)
                var[feed] = db[feed]


    for feed in args['feeds']:
        # determine if run has been completed
        run = checkRun(args['type'], list(C[feed].keys()), list(dG[feed].keys()))
        # gen data
        numIndep = gen_data[feed][run]['numIndep']
        # concentrations
        c_mol, box, conc = getConc(C[feed])
        if args['verbosity'] > 0:
            print('liquid analysis for feed %s is with mol %s, box %s'%(feed, mol, box))
        # initialize variables if needed
        mol_data = {}
        solutes = sorted([i for i in N[feed][run].keys()
                                    if ((N[feed][run][i]['box2']['mean'] > 1e-06)
                          and (N[feed][run][i]['box2']['mean'] < 300))      ])
        for mol in solutes:
            if mol not in mol_data.keys():
                mol_data[mol] = {}
                for key in ['deltaG', 'concentrations']:
                    mol_data[mol][key] = {'mean':[],'95conf':[],'feed':[]}
            deltaG = mol_data[mol]['deltaG']
            concentrations = mol_data[mol]['concentrations']
            concentrations['mean'].append( conc['mean'] )
            concentrations['95conf'].append( calc95conf(conc['stdev'], numIndep) )
            concentrations['feed'].append( feed )
            dG_mean, dG_stdev = (dG[feed][run][mol]['box3--box2']['mean'],
                                dG[feed][run][mol]['box3--box2']['stdev'])
            deltaG['mean'].append( dG_mean )
            deltaG['95conf'].append( calc95conf(dG_stdev, numIndep) )
            deltaG['feed'].append( feed )
        # write out results after each feed
        for mol in mol_data.keys():
            file_name = 'dG-mol%s_vs_C-mol%s.dat'%(mol, c_mol)
            file_name = file_name.replace('/','_')
            x_data = mol_data[mol]['concentrations']
            y_data = mol_data[mol]['deltaG']
            writeAGR(x_data['mean'], y_data['mean'],
                     x_data['95conf'], y_data['95conf'], x_data['feed'], file_name)

