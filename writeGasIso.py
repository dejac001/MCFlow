'''
Write isotherm from previously generated databank files
'''
from runAnalyzer import checkRun
from statistics import calc95conf
from file_formatting.writer import writeAGR
from chem_constants import R, N_av
import math

if __name__ == '__main__':
    from analysis_parsers import Plot
    import json 

    my_parser = Plot()
    # TODO: decide to do kH adsorption if kH provided? -- gas Iso & kH iso in same file?
    my_parser.parser.add_argument('-kH','--henry',help='kH [mean, uncertainty]',
                                  type=float, nargs='+')
    my_parser.isotherm(gas=True)

    args = vars(my_parser.parse_args())
    assert args['units'], 'Units must be defined for isotherm'
    assert len(args['henry']) != 1, 'Must input kH(mean) and kH(stdev)'
    assert args['box'], 'Box needed for pressure in gas phase isotherm'

    N = {}; P = {}; gen_data = {}
    for file, var in zip(['N-data.json','P-data.json', 'general-data.json'],
                         [N, P, gen_data]):
        with open('%s/%s'%(args['path'], file)) as f:
            db = json.load(f)
        for feed in args['feeds']:
            assert feed in db.keys(), 'Feed {} not in database'.format(feed)
            var[feed] = db[feed]

    # mol_data = {}
    # for feed in args['feeds']:
    #     # determine if run has been completed
    #     run = checkRun(args['type'], list(P[feed].keys()), list(N[feed].keys()))
    #     # gen data
    #     numIndep = gen_data[feed][run]['numIndep']
    #     # concentrations
    #     mol, box, conc = getConc(P[feed])
    #     if args['verbosity'] > 0:
    #         print('liquid analysis for feed %s is with mol %s, box %s'%(feed, mol, box))
    #     # initialize variables if needed
    #     mols_adsorbed = sorted([i for i in N[feed][run].keys()
    #                                 if N[feed][run][i]['box1']['mean'] > 1e-06])
    #     for mol in mols_adsorbed:
    #         if mol not in mol_data.keys():
    #             mol_data[mol] = {}
    #             for key in ['loading', 'concentrations']:
    #                 mol_data[mol][key] = {'mean':[],'95conf':[],'feed':[]}
    #         loadings = mol_data[mol]['loading']
    #         concentrations = mol_data[mol]['concentrations']
    #         if args['units'] == 'molec/uc':
    #             qfactor = 1/gen_data[feed][run]['zeolite']['unit cells']
    #         elif args['units'] == 'g/g':
    #             qfactor = (gen_data[feed][run]['molecular weight'][args['mol']]/
    #                             N_av)/gen_data[feed]['zeolite']['mass (g)']
    #         elif args['units'] == 'mol/kg':
    #             qfactor = gen_data[feed][run]['zeolite'][' mol/kg / 1 mlcl adsorbed']
    #         loadings['mean'].append( N[feed][run][mol]['box1']['mean']*qfactor )
    #         loadings['95conf'].append( calc95conf( N[feed][run][mol]['box1']['stdev']*qfactor ) )
    #         loadings['feed'].append( feed )
    #         concentrations['mean'].append( conc['mean'] )
    #         concentrations['95conf'].append( calc95conf(conc['stdev'], numIndep)  )
    #         concentrations['feed'].append( feed )
    #     # write after each feed!
    #     for mol in mol_data.keys():
    #         iso_name = 'isotherm-mols-%s-mol%s-%s'%('_'.join(mols_adsorbed),mol, args['units'])
    #         iso_name = iso_name.replace('/','_')
    #         x_data = mol_data[mol]['concentrations']
    #         y_data = mol_data[mol]['loadings']
    #         writeAGR(x_data['mean'], y_data['mean'],
    #                  x_data['95conf'], y_data['95conf'],
    #                  x_data['feed'], iso_name)

