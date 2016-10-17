'''
Write isotherm from previously generated databank files
'''
T_values = {'8':2.365,'4':3.182}

if __name__ == '__main__':
    from parser import Plot
    import shelve
    from chem_constants import N_av

    my_parser = Plot()
    my_parser.isotherm()

    args = vars(my_parser.parse_args())
    assert args['units'], 'Units must be defined for isotherm'

    with shelve.open(args['path'] + '/N-data.db') as db:
        N = {}
        for feed in args['feeds']:
            try:
                N[feed] = db[feed]
            except KeyError:
                print('Feed {} not in database'.format(feed))

    try:
        c_box = 'box2'
        with shelve.open(args['path'] + '/Conc-data.db') as db:
            C = {}
            assert db.keys(), 'No liquid concentrations found'
            for feed in args['feeds']:
                try:
                    C[feed] = db[feed]
                except KeyError:
                    print('Feed {} not in database'.format(feed))
    except AssertionError:
        c_box = 'box3'
        print('Using pressure in kPa from box 3 as concentrations')
        with shelve.open(args['path'] + '/P-data.db') as db:
            C = {}
            for feed in args['feeds']:
                try:
                    C[feed] = db[feed]
                except KeyError:
                    print('Feed {} not in database'.format(feed))

    with shelve.open(args['path'] + '/general-data.db') as db:
        gen_data = {}
        for feed in args['feeds']:
            try:
                gen_data[feed] = db[feed]
            except KeyError:
                print('Feed {} not in database'.format(feed))

    # TODO: reorganize so that can plot isotherms with same x-axis
    for mol in args['mol']:
        loadings = {'mean': [],'95conf':[],'feed':[]}
        concentrations = {'mean': [],'95conf':[],'feed':[]}
        for feed in N.keys():
            for run in N[feed].keys():
                if (args['type'] in run) and (run in C[feed].keys()):
                    numIndep = gen_data[feed][run]['numIndep']
                    if args['units'] == 'molec/uc':
                        qfactor = 1/gen_data[feed][run]['zeolite']['unit cells']
                    elif args['units'] == 'g/g':
                        qfactor = (gen_data[feed][run]['molecular weight'][args['mol']]/N_av)/gen_data[feed]['zeolite']['mass (g)']
                    elif args['units'] == 'mol/kg':
                        qfactor = gen_data[feed][run]['zeolite'][' mol/kg / 1 mlcl adsorbed']
                    loadings['mean'].append( N[feed][run][mol]['box1']['mean']*qfactor )
                    loadings['95conf'].append( N[feed][run][mol]['box1']['stdev']*qfactor
                                                /pow(numIndep, 1/2)*T_values['%i'%numIndep] )
                    loadings['feed'].append( feed )
                    if c_box == 'box2':
                        c_mean, c_stdev = C[feed][run][mol][c_box]['mean'], C[feed][run][mol][c_box]['stdev']
                    elif c_box == 'box3':
                        # gas phase adsorption
                        c_mean, c_stdev = C[feed][run]['box2']['mean'], C[feed][run]['box2']['stdev']
                    concentrations['mean'].append( c_mean )
                    concentrations['95conf'].append( c_stdev
                                                 /pow(numIndep, 1/2)*T_values['%i'%numIndep] )
                    concentrations['feed'].append( feed )
        iso_name = 'isotherm-mol%s-%s'%(mol,args['units'])
        with open(iso_name.replace('/','_'),'a') as f:
            for c, q, dc, dq, name in zip(concentrations['mean'], loadings['mean'], concentrations['95conf'],
                                    loadings['95conf'], loadings['feed']):
                f.write('%e %e %e %e %s\n'%(c, q, dc, dq, name))
