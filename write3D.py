__author__ = 'dejacor'
def write(data, file, nIndep):
    with open(file,'w') as f:
        for feed in data['x'].keys():
            x, dx = data['x'][feed]['mean'], calc95conf(data['x'][feed]['stdev'], nIndep)
            y, dy = data['y'][feed]['mean'], calc95conf(data['y'][feed]['stdev'], nIndep)
            if feed in data['z'].keys():
                z, dz = data['z'][feed]['mean'], calc95conf(data['z'][feed]['stdev'], nIndep)
                my_line = '%e %e %e %e %e %e \n'%(x,y,z,dx,dy,dz)
            else:
                my_line = '%e %e %e %e\n'%(x,y,dx,dy)
            f.write(my_line)

def getMolAds(num_molec):
    mols_adsorbed = []
    for mol in map(str,sorted(map(int,num_molec.keys()))):
        mean, stdev = (num_molec[mol]['box1']['mean'], num_molec[mol]['box1']['stdev'])
        ntotal = sum(num_molec[mol][box]['mean'] for box in num_molec[mol].keys())
        if (mean > 0.) and (mean < ntotal):
            mols_adsorbed.append(mol)
    assert mols_adsorbed, 'No molecules adsorbed. Mean: %3.2e, stdev: %3.2e'%(mean, stdev)
    if not mols_adsorbed: mols_adsorbed = ['1']
    return mols_adsorbed

def findVapBox(num_dens, mol):
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
    print('For mol%s, vapor box determined to be %s'%(mol, my_box))
    return my_box

class Plotter:
    def __init__(self, **kwargs):
        # TODO: make general db reader that all files can do
        self.N = {}; self.rho = {}; self.gen_data = {}; self.dG = {}; self.K = {}; self.P = {}; self.X = {}
        self.U = {}; self.boxlx = {}; self.C = {}
        self.files = ['N-data.json','rho-data.json','general-data.json','dG-data.json',
                        'K-data.json','P-data.json','X-data.json',
                      'U-data.json','boxlx-data.json','Conc-data.json']
        self.variables = [self.N, self.rho, self.gen_data, self.dG, self.K, self.P, self.X,
                          self.U, self.boxlx, self.C]

        self.path = kwargs['path']
        self.feeds = kwargs['feeds']
        if 'units' in kwargs.keys():
#            assert kwargs['mol'], 'Mol needed for number density to calculate P assuming I.G.'
#            assert kwargs['Temp'], 'Temperature needed for number density to calculate P assuming I.G.'
#            self.xlabel = ['Pig-mol%s_%4.1f'%(kwargs['mol'],kwargs['Temp']),'dP']
#             self.T = kwargs['Temp']
#             self.mol = kwargs['mol']
#             self.indep = kwargs['indep']
            self.units = kwargs['units']
            # self.zeoVolume = kwargs['ZeoVolume']
            # self.box = kwargs['box']
            # self.boxes = kwargs['boxes']
            # self.film = kwargs['film']

    def read_json(self):
        files_to_remove = []
        for file, var in zip(self.files, self.variables):
            with open('%s/%s'%(self.path, file), 'w') as f:
                data = json.load(f)
            for feed in self.feeds:
                if feed not in data.keys():
                    print('Feed {} not in {}'.format(feed, file))
                    files_to_remove.append( self.files.index(file))
                else:
                    var[feed] = data[feed]
        self.files = [self.files[i] for i in range(len(self.files)) if i not in files_to_remove]
        self.variables = [self.variables[i] for i in range(len(self.variables)) if i not in files_to_remove]

    def number_density(self, mol, box):
        num_dens = self.rho[self.feed][self.run]

        return {
            key:num_dens[mol][box][key] for key in ('mean','stdev','raw')
        }

    def deltaG(self, mol, boxFrom, boxTo):
        deltaG = self.dG[self.feed][self.run]
        boxPair = boxFrom+'--'+boxTo
        data = {key:deltaG[mol][boxPair][key] for key in ('mean','stdev','raw')}
        return data

    def Pig(self, mol, box):
        rho = self.number_density(self, mol, box)
        factor = 1/N_av*R['nm**3*kPa/(mol*K)']*self.T
        return {
            key:rho[key]*factor for key in ('mean','stdev','raw')
        }

    def Q(self, mol):
        N = self.N[self.feed][self.run]
        gen_data = self.gen_data[self.feed][self.run]
        zeolite_mass_g = gen_data['zeolite']['mass (g)']
        if self.units == 'molec/uc':
            qfactor = 1/gen_data['zeolite']['unit cells']
        elif self.units == 'g/g':
            qfactor = ((gen_data['molecular weight'][mol]/N_av) /
                        zeolite_mass_g )
        elif self.units == 'mol/kg':
            qfactor = gen_data['zeolite'][' mol/kg / 1 mlcl adsorbed']
        data = {key:N[mol]['box1'][key]*qfactor for key in ('mean','stdev')}
        data['raw'] = [i*qfactor for i in N[mol]['box1']['raw']]
        return data

    def R(self,mol,box):
        N = self.N[self.feed][self.run]
        nIndep = self.gen_data[self.feed][self.run]['numIndep']
        N_tot_raw = [0 for i in range(nIndep)]
        for box, vals in N[mol].items():
            for j in range(len(vals['raw'])):
                N_tot_raw += vals['raw'][j]
        N_tot_raw = [i/nIndep for i in N_tot_raw]
        R_raw = [N[mol][box]['raw'][i]/N_tot_raw[i]*100 for i in range(len(N_tot_raw))]
        return {'mean':np.mean(R_raw), 'stdev':np.std(R_raw), 'raw':R_raw}

    def moleFrac(self,mol,box):
        mol_frac = self.X[self.feed][self.run]
        try:
            data = {key:mol_frac[box][mol][key] for key in ('raw','mean','stdev')}
        except KeyError:
            print('Molecules in feed {} only {}'.format(self.feed, mol_frac[box].keys()))
            data = {key:1e-10 for key in ('raw','mean','stdev')}
            data['raw'] = [1e-15]
        return data

    def Kratio(self, box, mol1, mol2):
        Kfrac = self.K[self.feed][self.run]
        mol_pair = '%s/%s'%(mol1,mol2)
        return {
            key:Kfrac[box][mol_pair][key] for key in ('raw','mean','stdev')
        }

    def Selec(self, boxFrom, boxTo, mol1, mol2):
        K_from = self.Kratio(boxFrom, mol1, mol2)
        K_to = self.Kratio(boxTo, mol1, mol2)
        S_mean = K_to['mean'] / K_from['mean']
        S_stdev = eProp_division(K_to['mean'], K_to['stdev'],
                                      K_from['mean'], K_from['stdev'])
        S_raw = [i/j for i,j in zip(K_to['raw'],K_from['raw'])]
        return {'mean':S_mean,'stdev':S_stdev,'raw':S_raw}

    def Pbox(self, box):
        return self.P[self.feed][self.run][box]

    def dU(self, boxFrom, boxTo):
        U = self.U[self.feed][self.run]
        N = self.N[self.feed][self.run]
        nIndep = self.gen_data[self.feed][self.run]['numIndep']
        N_tot_boxFrom_raw = [0. for i in range(nIndep)]
        N_tot_boxTo_raw = [0. for i in range(nIndep)]
        for mol in N.keys():
            for j in range(len(N[mol][boxFrom]['raw'])):
                N_tot_boxFrom_raw[j] += N[mol][boxFrom]['raw'][j]
                N_tot_boxTo_raw[j] += N[mol][boxTo]['raw']['j']
        dU_raw = []
        for i in range(nIndep):
            du = (U[boxTo]['raw'][i]/N_tot_boxTo_raw[i] - U[boxFrom]['raw'][i]/N_tot_boxFrom_raw[i])*8.314/1000.
            dU_raw.append(du)
        return {
            'mean':np.mean(dU_raw),
            'stdev':np.st(dU_raw),
            'raw':dU_raw
        }

    def Pi(self, mol, box):
        P = self.Pbox(box)
        y = self.moleFrac(mol, box)
        mean = P['mean']*y['mean']
        raw = [Pbox*yi for Pbox,yi in zip(P['raw'],y['raw'])]
        stdev = mean*np.sqrt(
            pow(P['stdev']/P['mean'],2) + pow(y['stdev']/y['mean'],2)
        )
        return {'mean':mean, 'stdev':stdev, 'raw':raw}

    def T(self):
        T = self.gen_data[self.feed][self.run]['temperature']
        nIndep = self.gen_data[self.feed][self.run]['numIndep']
        return {'mean':T,'stdev':0.,
                'raw':[T for i in range(nIndep)]}

    def Conc(self):
        def getC(c_data, mol='', box=''):
            for key, value in sorted(c_data.items()):
                if (len(key) == 1) or (len(key) == 2):
                    mol = key
                elif 'box' in key:
                    box = key
                if isinstance(value, dict):
                    if 'mean' in value.keys():
                        return value
                    else:
                        return getC(value, mol, box)
        return getC(self.C[self.feed])

    def convertChoice(self, choice, mol='',box='',boxFrom='',
                      boxTo='', mol1='',mol2=''):
        if 'box' not in box: box = 'box' + box
        if choice == 'rho':
            return self.number_density(mol, box)
        elif choice == 'dG':
            return self.deltaG(mol, boxFrom, boxTo)
        elif choice == 'Pig':
            return self.Pig(mol, box)
        elif choice == 'Q':
            return self.Q(mol)
        elif choice == 'X':
            return self.moleFrac(mol, box)
        elif choice == 'K':
            return self.Kratio(box, mol1, mol2)
        elif choice == 'R':
            return self.R(mol, box)
        elif choice == 'S':
            return self.Selec(boxFrom, boxTo, mol1, mol2)
        elif choice == 'Pbox':
            return self.Pbox(box)
        elif choice == 'dU':
            return self.dU(boxFrom, boxTo)
        elif choice == 'Pi':
            return self.Pi(mol, box)
        elif choice == 'T':
            return self.T()
        elif choice == 'C':
            return self.Conc()

from MCFlow.runAnalyzer import checkRun, calc95conf
from MCFlow.chem_constants import N_av, R
from MCFlow.calc_tools import eProp_division
#from MCFlow.writeXvY import calculateMoleFrac, calculateIGMolFrac, calculateS_ab
import numpy as np
import os, math, json

from MCFlow.parser import Plot, Main

class Plot3D(Plot):
    def __init__(self):
        Main.__init__(self)
        choices = ['rho','dG','Pig','Q','X','K','S','Pbox','dU','Pi','T','C','','R']
        self.parser.description = 'Plot results in 3D'
        self.parser.add_argument('-mx','--molx',help='Molecule associated with x axis', type=str)
        self.parser.add_argument('-my','--moly',help='Molec w/ y axis',type=str)
        self.parser.add_argument('-mz','--molz',help='Molec w/ z axis',type=str)
        self.parser.add_argument('-bx','--boxx',help='Box associated with x axis', type=str)
        self.parser.add_argument('-by','--boxy',help='Box associated with y axis', type=str)
        self.parser.add_argument('-bz','--boxz',help='Box associated with z axis', type=str)
        self.parser.add_argument('-x','--xaxis',choices=choices)
        self.parser.add_argument('-y','--yaxis',choices=choices)
        self.parser.add_argument('-z','--zaxis',choices=choices)
        self.parser.add_argument('-B','--boxes',help='box numbers to use',
                                 type=str,nargs='+',default=['box3','box2'])
        self.parser.add_argument('-M','--molecules',help='mol numbers to use',
                                 type=str,nargs='+',default=['1','2'])
        self.parser.add_argument('-u','--units',help='choices for units on plot',
                         choices=['molec/uc','g/g','mol/kg','None',
                                  '(mol/mol)/(mol/mol)','(mol/mol)/(kPa/kPa)'])
    # def isotherm(self):
    #     self.parser.add_argument('-kH','--henry',help='Henry constant (g/mL/kPa) '
    #                                                    'fmt: (mean, 95%%conf.) ',
    #                               type=float, nargs= '+')

    #     self.parser.add_argument('-V','--ZeoVolume',help='Volume of zeolite box (\AA**3)',
    #                              type = float, default = 40.044*39.798*40.149 )
    # def kH(self):
    #     self.parser.add_argument('-y', '--yaxis', help='y axis of plot',
    #                              choices=['Q','S','X'])
    #     self.parser.add_argument('-d','--density',help='infinite dilution density',type=float,
    #                              nargs='+',default=[0.981,0.006]
    #                              # Jorgensen and Jenson, J. Comput. Chem.
    #                              )
    #     self.parser.add_argument('-MW','--molecWeight',help='molecular weight of solvent',
    #                              type=float, default=18.02)

if __name__ == '__main__':
    # TODO: find way to have conditional arguments based off of what xaxis and yaxis are
    # (in other words, we have an excess of variables as arguments)
    my_parser = Plot3D()

    args = vars(my_parser.parse_args())
    assert args['yaxis'], 'No y axis chosen for plot'
    assert args['xaxis'], 'No x axis chosen for plot'

    my_plotter = Plotter(**args)
    my_plotter.read_json()


    boxFrom, boxTo = args['boxes']
    mol1, mol2 = args['molecules']
    kwargs = {'boxFrom':boxFrom,'boxTo':boxTo,'mol1':mol1,'mol2':mol2}

    data = {'x':{},'y':{},'z':{}}
    for feed in args['feeds']:
        file_name = ''
        # determine if run has been completed
        my_plotter.run = checkRun(args['type'], my_plotter.variables, feed)
        my_plotter.feed = feed
        for axis in ['x','y','z']:
            kwargs['mol'] = args['mol%s'%axis]
            kwargs['box'] = args['box%s'%axis]
            file_name += 'mol%s-box%s-%s_'%(kwargs['mol'],kwargs['box'],args['%saxis'%axis])
            data[axis][feed] = my_plotter.convertChoice(args['%saxis'%axis], **kwargs)
    write(data,file_name + '.dat', my_plotter.gen_data[my_plotter.feed][my_plotter.run]['numIndep'])
