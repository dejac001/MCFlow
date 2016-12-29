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


def calculateS_ab(Na, Nb, boxFrom='box2', boxTo='box1'):
    '''
    z = boxTo
    cbox = boxFrom
    S_ab  = [N(z,a)/N(z,tot)/(N(z,b)/N(z,tot))]   /   [N(cbox,a)/N(cbox,tot)/(N(cbox,b)/N(cbox,tot))]
    S_ab  = [N(z,a)/N(z,b)] / [N(cbox,a)/N(cbox,b)] = N(z,a)*N(cbox,b)/(N(z,b)*N(cbox,a))
    '''
    mean = (Na[boxTo]['mean']/Nb[boxTo]['mean']) / (Na[boxFrom]['mean']/Nb[boxFrom]['mean'])
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

class IdealGasAds:
    def __init__(self, **kwargs):

        # TODO: make general db reader that all files can do
        self.N = {}; self.rho = {}; self.gen_data = {}
        self.files = ['N-data.db','rho-data.db','general-data.db']
        self.variables = [self.N, self.rho, self.gen_data]
        if kwargs:
            assert kwargs['mol'], 'Mol needed for number density to calculate P assuming I.G.'
            assert kwargs['Temp'], 'Temperature needed for number density to calculate P assuming I.G.'
            self.xlabel = ['Pig-mol%s(kPa)'%kwargs['mol'],'dP']
            self.T = kwargs['Temp']
            self.mol = kwargs['mol']
            self.units = kwargs['units']
            self.feeds = kwargs['feeds']
            self.path = kwargs['path']

    def readDBs(self):
        for file, var in zip(self.files, self.variables):
            with shelve.open('%s/%s'%(self.path, file)) as db:
                for feed in self.feeds:
                    assert feed in db.keys(), 'Feed {} not in database for file {}'.format(feed, file)
                    var[feed] = db[feed]

    def getMolAds(self, num_molec):
        mols_adsorbed = []
        for mol in map(str,sorted(map(int,num_molec.keys()))):
            mean, stdev = (num_molec[mol]['box1']['mean'], num_molec[mol]['box1']['stdev'])
            ntotal = sum(num_molec[mol][box]['mean'] for box in num_molec[mol].keys())
            if (mean > 0.) and (mean < ntotal):
                mols_adsorbed.append(mol)
        assert mols_adsorbed, 'No molecules adsorbed. Mean: %3.2e, stdev: %3.2e'%(mean, stdev)
        if not mols_adsorbed: mols_adsorbed = ['1']
        return mols_adsorbed

    def getX(self, number_density, temperature):
        '''
        :param rho: number density in molec/nm**3
        :param T: temperature (K)
        :return:
        '''
        p_mean = number_density['mean']/N_av*R['nm**3*kPa/(mol*K)']*temperature
        p_stdev = number_density['stdev']/N_av*R['nm**3*kPa/(mol*K)']*temperature
        return {'mean':p_mean,'stdev':p_stdev}

    def findVapBox(self, num_dens, mol):
        rho_vap = 1. # molec/nm**3
        for box, value in num_dens[mol].items():
            if ((box != 'box1') and
                    (value['mean'] > 1.25e-12) and
                    (value['mean'] < rho_vap)):
                # find non-zeolite box with minimum number density
                # that corresponds to a pressure above 1e-10 bar
                rho_vap = value['mean']
                my_box = box
        assert my_box, 'Vapor box not found, maybe need to decrease min value'
        print('For mol%s, obtaining pressure'
              ' from num dens in %s'%(mol, my_box))
        return my_box
    
    def QvX(self, feed, run):
        '''
        :var X: either solution concentration (g/mL) or pressure (kPa)
        '''   
        N = self.N[feed][run]
        gen_data = self.gen_data[feed][run]
        file_description = '%s    Q(%s)    %s     dQ'%(self.xlabel[0], self.units,
                                                       self.xlabel[1])
        if 'Pig' in self.xlabel[0]:
            vapor_box = self.findVapBox( self.rho[feed][run], self.mol)
            # todo: need to choose molecule for x axis here
            X = self.getX(self.rho[feed][run][self.mol][vapor_box], self.T)
        elif 'C' in self.xlabel[0]:
            X, cmol, cbox = self.getX(self.C[feed])
        elif 'P-box' in self.xlabel[0]:
            X = self.P[feed][run]['box%s'%self.vapor_box]
        nIndep = gen_data['numIndep']
        mols_adsorbed = self.getMolAds(N)
        for mol in mols_adsorbed:
            file_name  = 'Qmol%s-%s-w%sin-zeo-vs-%s.dat'%(mol,self.units,''.join(mols_adsorbed),
                                                           self. xlabel[0][:self.xlabel[0].find('(')])
            if self.units == 'molec/uc':
                qfactor = 1/gen_data['zeolite']['unit cells']
            elif self.units == 'g/g':
                qfactor = ((gen_data['molecular weight'][mol]/N_av) / 
                            gen_data['zeolite']['mass (g)'] )
            elif self.units == 'mol/kg':
                qfactor = gen_data['zeolite'][' mol/kg / 1 mlcl adsorbed']
            Q_mean, Q_stdev = (N[mol]['box1']['mean']*qfactor, N[mol]['box1']['stdev']*qfactor)
            writeAGR([X['mean']],[Q_mean],
                     [calc95conf(X['stdev'], nIndep)], [calc95conf(Q_stdev, nIndep)],
                     [feed], file_name, file_description)
    
    
    def SvX(self, feed,  run):
        '''
        this is a special case of SvC, where we plot P as opposed to C
        and the boxFrom might not be the same for each molecule
        (ex: in the case of IG or kH adsorption)
        :param P: Pressure of x axis (note: not pressure data)
        '''
        def getPdata(num_molec, num_dens, mola, molb):
            my_data = {mola:{}, molb:{}}
            for my_mol in [mola, molb]:
                # input adsorbed info
                my_data[my_mol]['boxTo'] = num_molec[my_mol]['box1']
                # find box for pressure info
                my_box = IdealGasAds().findVapBox(num_dens, my_mol)
                # input pressure info
                my_data[my_mol]['boxFrom'] = IdealGasAds().getX(num_dens[my_mol][my_box], self.T)
            return my_data
        N = self.N[feed][run]
        rho = self.rho[feed][run]
        if 'Pig' in self.xlabel[0]:
            vapor_box = self.findVapBox( self.rho[feed][run], self.mol)
            X = self.getX(self.rho[feed][run][self.mol][vapor_box], self.T)
        elif 'C' in self.xlabel[0]:
            X, cmol, cbox = self.getX(self.C[feed])
        elif 'P-box' in self.xlabel[0]:
            X = self.P[feed][run]['box%s'%self.vapor_box]
        file_description = '%s    S(%s)    %s     dS'%(self.xlabel[0], self.units, self.xlabel[1])
        mols_adsorbed = self.getMolAds(N)
        mols_adsorbed.reverse() # start from higher number molecules (i.e. solutes)
        nIndep = self.gen_data[feed][run]['numIndep']
        for imol, mol1 in enumerate(mols_adsorbed):
            for mol2 in mols_adsorbed[(imol+1):]:
                s_name = '%s/%s'%(mol1,mol2)
                if 'P' in self.xlabel[0]:
                    box_from, box_to = 'boxFrom','boxTo'
                    s_data = getPdata(N, rho, mol1, mol2)
                elif 'C' in self.xlabel[0]:
                    box_from, box_to = 'box2', 'box1'
                    s_data = N
                elif 'P-box' in self.xlabel[0]:
                    print('P-box not implemented for SvX yet')
                S_mean, S_stdev = calculateS_ab(s_data[mol1],s_data[mol2],
                                                boxFrom=box_from, boxTo=box_to)
                file_name = 'S_%s-vs-%s.dat'%(s_name, self.xlabel[0])
                writeAGR([X['mean']],[S_mean],
                         [calc95conf(X['stdev'], nIndep)], [calc95conf(S_stdev, nIndep)],
                         [feed], file_name, file_description)

class GasAds(IdealGasAds):
    def __init__(self, **kwargs):
        self.N = {}; self.P = {}; self.gen_data = {};
        self.files = ['N-data.db','P-data.db','general-data.db']
        self.variables = [self.N, self.P, self.gen_data]
        if kwargs:
            assert kwargs['box'], 'Box needed to calculate pressure (box)'
            self.xlabel = ['P-box%s(kPa)'%kwargs['box'],'dP']
            self.vapor_box = kwargs['box']
            self.units = kwargs['units']
            self.feeds = kwargs['feeds']
            self.path = kwargs['path']

class LiqAds(IdealGasAds):
    def __init__(self, **kwargs):
        self.dG = {}; self.C = {}; self.N = {}; self.rho = {}; self.gen_data = {}
        self.files = ['dG-data.db','Conc-data.db','N-data.db','rho-data.db','general-data.db']
        self.variables = [self.dG, self.C, self.N, self.rho, self.gen_data]
        if kwargs:
            self.xlabel = ['C(g/mL)','dC']
            self.units = kwargs['units']
            self.feeds = kwargs['feeds']
            self.path = kwargs['path']
    
    def dGvX(self, feed):
        file_description = '%s    dG(kJ/mol)    %s     dG'%(self.xlabel[0],self.xlabel[1])
        N = self.N[feed][run]
        dG = self.dG[feed][run]
        nIndep = self.gen_data[feed][run]['numIndep']
        X, cmol, cbox = self.getX(self.C[feed])
        solutes = sorted([i for i in N.keys()
                                    if ((N[i]['box2']['mean'] > 1e-06)
                          and (N[i]['box2']['mean'] < 300))      ])
        for mol in solutes:
            file_name = 'dG-mol%s_vs_%s.dat'%(mol, self.xlabel[0])
            dG_mean, dG_stdev = (dG[mol]['box3--box2']['mean'],
                                dG[mol]['box3--box2']['stdev'])
            writeAGR([X['mean']],[dG_mean],
                     [calc95conf(X['stdev'], nIndep)], [calc95conf(dG_stdev, nIndep)],
                     [feed], file_name, file_description)

    def getX(self, c_data, mol='', box=''):
        for key, value in sorted(c_data.items()):
            if (len(key) == 1) or (len(key) == 2):
                mol = key
            elif 'box' in key:
                box = key
            if isinstance(value, dict):
                if 'mean' in value.keys():
                    return value, mol, box
                else:
                    return LiqAds().getX(value, mol, box)

from MCFlow.runAnalyzer import checkRun, calc95conf
from MCFlow.chem_constants import N_av, R
from MCFlow.parser import Plot
import os, math, shelve

if __name__ == '__main__':
    # TODO: find way to have conditional arguments based off of what xaxis and yaxis are
    # (in other words, we have an excess of variables as arguments)
    my_parser = Plot()
    my_parser.axes()
    my_parser.isotherm()

    args = vars(my_parser.parse_args())
    assert args['yaxis'], 'No y axis chosen for plot'
    assert args['xaxis'], 'No x axis chosen for plot'


    if args['xaxis'] == 'C':
        my_plotter = LiqAds(**args)
    elif args['xaxis'] == 'Pig':
        # TODO: add alternative way to calculate P w/ P-data.db & using mole fraction in box
        my_plotter = IdealGasAds(**args)
    elif args['xaxis'] == 'Pbox':
        my_plotter = GasAds(**args)
    my_plotter.readDBs()

    for feed in args['feeds']:
        # determine if run has been completed
        run = checkRun(args['type'], my_plotter.variables, feed)
        if args['yaxis'] == 'Q':
            my_plotter.QvX(feed, run)
        elif args['yaxis'] == 'S':
            my_plotter.SvX(feed, run)
        elif args['yaxis'] == 'dG':
            my_plotter.dGvX(feed)
