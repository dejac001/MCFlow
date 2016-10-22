def writeAGR(x, y, dx, dy, names, file, description):
    file = file.replace('/','_')
    #TODO: when writing to file, sort all lines (previously written too) by x-value
    if os.path.isfile(file):
        # don't add description if already written to
        description = ''
    else:
        assert description, 'Description needed before writing to file %s'%file
    with open(file,'a') as f:
        if description:
            f.write('# %s\n'%description)
        for i in sorted(x):
            my_index = x.index(i)
            my_line = ''
            if dx and dy:
                for data in (x, y, dx, dy):
                    my_line += '%e '%data[my_index]
            elif x and y:
                for data in (x, y):
                    my_line += '%e '%data[my_index]
            my_line += '#%s\n'%names[my_index]
            assert my_line.count('\n') == 1, 'Too many line endings %s'%my_line
            f.write(my_line)

def getMolAds(num_molec):
    mols_adsorbed = []
    for mol in map(str,sorted(map(int,num_molec.keys()))):
        mean = num_molec[mol]['box1']['mean']
        ntotal = sum(num_molec[mol][box]['mean'] for box in num_molec[mol].keys())
        if (mean > 1e-06) and (mean < ntotal):
            mols_adsorbed.append(mol)
    return mols_adsorbed

def calculateS_ab(Na, Nb, boxFrom='box2', boxTo='box1'):
    '''
    z = boxTo
    cbox = boxFrom
    S_ab  = [N(z,a)/N(z,tot)/(N(z,b)/N(z,tot))]   /   [N(cbox,a)/N(cbox,tot)/(N(cbox,b)/N(cbox,tot))]
    S_ab  = [N(z,a)/N(z,b)] / [N(cbox,a)/N(cbox,b)] = N(z,a)*N(cbox,b)/(N(z,b)*N(cbox,a))
    '''
    mean = (Na[boxTo]['mean']/Nb[boxTo]['mean']) / (Na[boxFrom]['mean']/Nb[boxTo]['mean'])
    dS_dNa = {}; dS_dNb = {}
    dS_dNa[boxTo]   = Nb[boxFrom]['mean']      /     (Na[boxFrom]['mean']*Nb[boxTo]['mean'])
    dS_dNb[boxFrom] = Na[boxTo]['mean']        /     (Na[boxFrom]['mean']*Nb[boxTo]['mean'])
    dS_dNb[boxTo]   = -Na[boxTo]['mean']*Nb[boxFrom]['mean']/(Na[boxFrom]['mean']*math.pow(Nb[boxTo]['mean'],2))
    dS_dNa[boxFrom] = -Na[boxTo]['mean']*Nb[boxFrom]['mean']/( Nb[boxTo]['mean']*math.pow(Na[boxFrom]['mean'],2))
    stdev = math.sqrt( math.pow(dS_dNa[boxFrom],2)*math.pow(Na[boxFrom]['stdev'],2)
                       + math.pow(dS_dNa[boxTo],2)*math.pow(Na[boxTo]['stdev'],2)
                       + math.pow(dS_dNb[boxFrom],2)*math.pow(Nb[boxFrom]['stdev'],2)
                       + math.pow(dS_dNb[boxTo],2)*math.pow(Nb[boxTo]['stdev'],2)
                       )
    return mean, stdev

def DGvC(dG, N, C, data, feed, xlabel):
    file_description = '%s    dG(kJ/mol)    %s     dG'%(xlabel[0],xlabel[1])
    solutes = sorted([i for i in N.keys()
                                if ((N[i]['box2']['mean'] > 1e-06)
                      and (N[i]['box2']['mean'] < 300))      ])
    nIndep = data['numIndep']
    for mol in solutes:
        file_name = 'dG-mol%s_vs_%s.dat'%(mol, xlabel[0])
        dG_mean, dG_stdev = (dG[mol]['box3--box2']['mean'],
                            dG[mol]['box3--box2']['stdev'])
        writeAGR([C['mean']],[dG_mean],
                 [calc95conf(C['stdev'], nIndep)], [calc95conf(dG_stdev, nIndep)],
                 [feed], file_name, file_description)

def QvX(N, X, data, feed, run, units, xlabel):
    '''
    :param X: either solution concentration (g/mL) or pressure (kPa)
    '''
    file_description = '%s    Q(%s)    %s     dQ'%(xlabel[0], units,  xlabel[1])
    nIndep = data[run]['numIndep']
    mols_adsorbed = getMolAds(N)
    for mol in mols_adsorbed:
        file_name  = 'Qmol%s-%s-w%sin-zeo-vs-%s.dat'%(mol,units,''.join(mols_adsorbed),
                                                        xlabel[0][:xlabel[0].find('(')])
        if units == 'molec/uc':
            qfactor = 1/data[run]['zeolite']['unit cells']
        elif units == 'g/g':
            qfactor = (data[run]['molecular weight'][mol]/
                            N_av)/data[run]['zeolite']['mass (g)']
        elif units == 'mol/kg':
            qfactor = data[run]['zeolite'][' mol/kg / 1 mlcl adsorbed']
        Q_mean, Q_stdev = (N[mol]['box1']['mean']*qfactor, N[mol]['box1']['stdev']*qfactor)
        writeAGR([X['mean']],[Q_mean],
                 [calc95conf(X['stdev'], nIndep)], [calc95conf(Q_stdev, nIndep)],
                 [feed], file_name, file_description)

def SvC(N, C, units, feed, nIndep, xlabel):
    file_description = '%s    S(%s)    dC     %s'%(xlabel[0], units, xlabel[1])
    mols_adsorbed = getMolAds(N)
    mols_adsorbed.reverse() # start from higher number molecules (i.e. solutes)
    for imol, mol1 in enumerate(mols_adsorbed):
        for mol2 in mols_adsorbed[(imol+1):]:
            s_name = '%s/%s'%(mol1,mol2)
            S_mean, S_stdev = calculateS_ab(N[mol1], N[mol2])
            file_name = 'S_%s-vs-%s.dat'%(s_name, xlabel[0])
            writeAGR([C['mean']],[S_mean],
                     [calc95conf(C['stdev'], nIndep)], [calc95conf(S_stdev, nIndep)],
                     [feed], file_name, file_description)

def SvP(N, rho, P,  units, feed, nIndep, T, xlabel):
    '''
    this is a special case of SvC, where we plot P as opposed to C
    and the boxFrom might not be the same for each molecule
    (ex: in the case of IG or kH adsorption)
    :param P: Pressure of x axis (note: not pressure data)
    '''
    file_description = '%s    S(%s)    dC     %s'%(xlabel[0], units, xlabel[1])
    mols_adsorbed = getMolAds(N)
    mols_adsorbed.reverse() # start from higher number molecules (i.e. solutes)
    for imol, mol1 in enumerate(mols_adsorbed):
        for mol2 in mols_adsorbed[(imol+1):]:
            s_name = '%s/%s'%(mol1,mol2)
            my_data = {mol1:{}, mol2:{}}
            for my_mol in [mol1, mol2]:
                # input adsorbed info
                my_data[my_mol]['boxTo'] = N[my_mol]['box1']
                # find box for pressure info
                num_dens = 1. # molec/nm**3
                for box, value in rho[my_mol].items():
                    if ((box != 'box1') and
                            (value['mean'] > 2.25e-12) and
                            (value['mean'] < num_dens)):
                        # find non-zeolite box with minimum number density
                        # that corresponds to a pressure above 1e-10 bar
                        num_dens = value['mean']
                        my_box = box
                try:
                    print('For mol%s, obtaining pressure'
                          ' from num dens in %s'%(my_mol, my_box))
                except UnboundLocalError:
                    print('No gas box found for mol %s'%my_mol)
                    print(rho[my_mol].items())
                # input pressure info
                my_data[my_mol]['boxFrom'] = getPi_ig(rho[my_mol][my_box], T)
            S_mean, S_stdev = calculateS_ab(my_data[mol1],my_data[mol2],
                                            boxFrom='boxFrom',boxTo='boxTo')
            file_name = 'S_%s-vs-%s.dat'%(s_name, xlabel[0])
            writeAGR([P['mean']],[S_mean],
                     [calc95conf(P['stdev'], nIndep)], [calc95conf(S_stdev, nIndep)],
                     [feed], file_name, file_description)

from MCFlow.runAnalyzer import getConc, checkRun, calc95conf, getPi_ig
from MCFlow.chem_constants import N_av
import os, math

if __name__ == '__main__':
    from parser import Plot
    import shelve

    # TODO: find way to have conditional arguments based off of what xaxis and yaxis are
    # (in other words, we have an excess of variables as arguments)
    my_parser = Plot()
    my_parser.axes()
    my_parser.isotherm()

    args = vars(my_parser.parse_args())
    assert args['yaxis'], 'No y axis chosen for plot'
    assert args['xaxis'], 'No x axis chosen for plot'

    # TODO: make general db reader that all files can do
    num_molec = {}; rho = {}; gen_data = {}
    files = ['N-data.db','rho-data.db','general-data.db']
    variables = [num_molec, rho, gen_data]
    if args['xaxis'] == 'C':
        dG = {}; C = {};
        variables.insert(0, C)
        variables.insert(0, dG)
        files.insert(0, 'Conc-data.db')
        files.insert(0, 'dG-data.db')
    # read db
    for file, var in zip(files, variables):
        with shelve.open('%s/%s'%(args['path'], file)) as db:
            for feed in args['feeds']:
                assert feed in db.keys(), 'Feed {} not in database for file {}'.format(feed, file)
                var[feed] = db[feed]

    for feed in args['feeds']:
        # determine if run has been completed
        run = checkRun(args['type'], list(num_molec[feed].keys()), list(rho[feed].keys()))
        # gen data
        numIndep = gen_data[feed][run]['numIndep']
        # concentrations
        if args['xaxis'] == 'C':
            conc, c_mol, c_box = getConc(C[feed])
            x_label = ['C-mol%s(g/mL)'%c_mol,'dC']
        # if args['verbosity'] > 0:
        #     print('liquid analysis for feed %s is with mol %s, box %s'%(feed, mol, c_box))
        # initialize variables if needed
            if (args['yaxis'] == 'Q'):
                QvX(num_molec[feed][run], conc, gen_data[feed], feed,  run, args['units'], x_label)
            elif (args['yaxis'] == 'dG'):
                DGvC(dG[feed][run], num_molec[feed][run], conc, gen_data[feed], feed, x_label)
            elif (args['yaxis'] == 'S'):
                SvC(num_molec[feed][run], C[feed][run], args['units'],
                    args['feed'], gen_data[feed][run]['numIndep'], x_label)
        elif args['xaxis'] == 'Pig':
            # TODO: add alternative way to calculate P w/ P-data.db & using mole fraction in box
            assert args['mol'], 'Mol needed for number density to calculate P assuming I.G.'
            assert args['box'], 'Box needed for number density to calculate P assuming I.G.'
            assert args['Temp'], 'Temperature needed for number density to calculate P assuming I.G.'
            x_label = ['Pig-mol%s(kPa)'%args['mol'],'dP']
            press = getPi_ig(rho[feed][run][args['mol']]['box%s'%args['box']], args['Temp'])
            if args['yaxis'] == 'Q':
                QvX(num_molec[feed][run], press,gen_data[feed], feed,  run, args['units'], x_label)
            elif args['yaxis'] == 'S':
                SvP(num_molec[feed][run], rho[feed][run], press,
                    args['units'], feed, gen_data[feed][run]['numIndep'],
                    args['Temp'], x_label)

