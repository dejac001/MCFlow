class Main:
    def __init__(self):

        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("-v", "--verbosity", action="count", default=0)
        self.parser.add_argument('-t','--type',
                             help='type of run (equil- or prod-)',
                             type = str, default='equil-',choices = ['equil-','prod-'])
        self.parser.add_argument('-p','--path',
                                 help='main path to feeds (see --feeds)',
                                 type=str,default=os.getcwd())
        self.parser.add_argument('-f','--feeds',
                                 help ='state point (ex: temperature or pressure);'
                                             ' must be in ls ${path}',
                             nargs = '+', type=str)
        self.parser.add_argument('-i','--indep',
                        help='RANGE of independent simulations, if = 0 (for plotting only) plot all data',
                        nargs = '+',type=int, default = range(1,9))

    def other(self):
        self.parser.add_argument('-iv','--interval',
                        help='number of runs interval to analyze',
                        type=int,default = 0)
        self.parser.add_argument('-s','--guessStart',help ='run to guess start at',
                        type = int, default = 1)

    def parse_args(self):
        self.args = self.parser.parse_args()
        assert self.args.feeds, 'No feeds to analyze'
        if self.args.feeds == ['.']:
            self.args.indep = ['.']
        return self.args

class Plot(Main):
    def __init__(self):
        Main.__init__(self)
        self.parser.description = 'Plot results'
        # TODO: reformat mol pass in so pass in by str name (i.e. WATER or 15PDO)
        self.parser.add_argument('-m','--mol',help='Molecule to analyze', type=str)
        self.parser.add_argument('-TK','--Temp',help='Temperature in Kelvin',
                                  type=float)
        self.parser.add_argument('-b','--box',help='box number to analyze', type=str)

    def axes(self):
        self.parser.add_argument('-x','--xaxis', help='x axis of plot',
                                 choices = ['C','Pig', 'Pbox','Pi','Q','rho','x'])
        self.parser.add_argument('-y','--yaxis',help='x axis of plot',
                                 choices = ['Q','dG','S','R',
                                            'X','dH','dU','dHig','Pigxy','Txy'])
        self.parser.add_argument('-B','--boxes',help='box numbers to use',
                                 type=str,nargs='+',default=['box3','box2'])
    def isotherm(self):
        self.parser.add_argument('-u','--units',help='choices for units on plot',
                                 choices=['molec/uc','g/g','mol/kg',
                                          '(mol/mol)/(mol/mol)','(mol/mol)/(kPa/kPa)'])
        self.parser.add_argument('-V','--ZeoVolume',help='Volume of zeolite box (\AA**3)',
                                 type = float, default = 40.044*39.798*40.149 )
    def kH(self):
        self.parser.add_argument('-kH','--henry',help='Henry constant (g/mL/kPa) '
                                                       'fmt: (mean, 95%%conf.) ',
                                  type=float, nargs= '+')

class MultMols(Plot):
    def __init__(self):
        Main.__init__(self)
        Main.other(self)
        self.parser.description = 'Obtain results for multiple molecules'
        self.parser.add_argument('-m','--mol',help='MoleculeS to analyze', type=str, nargs='+')
        self.parser.add_argument('-TK','--Temp',help='Temperature in K', type=float)
        self.parser.add_argument('-b','--box',help='box to analyze', type=str)

class Results(Plot):
    def __init__(self):
        Plot.__init__(self)
        Main.other(self)
        self.parser.description = 'Obtain results for simulations'
    def multi(self):
        self.parser.add_argument('-M','--molecules',help='Molecules to analyze for trajectory',
                                  type = str, nargs = '+')
        self.parser.add_argument('-B','--boxes',help='box numbers to use',
                                 type=str,nargs='+',default=['box2'])
        

class Change(Results):
    def __init__(self):
        Main.__init__(self)
        Main.other(self)
        self.parser.description='Change MC probabilities, and potentially ' \
                                'obtain results from these runs'
        self.parser.add_argument('-T','--time',
                        help='time in seconds for next run',type = int,
                        default=96*60*60-600)
        self.parser.add_argument('-r','--rcut',
                        help='rcut fraction of boxlx to use for vapor box',
                        type=float, default=0.50)
        self.parser.add_argument('-N','--nstep', help='number of MCCs for next run',
                        type=int, default=50000)

class ChangeInput(Main):
    def __init__(self):
        Main.__init__(self)
        self.parser.add_argument('-I','--input',help='name of fort.4 file',type=str,
                        default= 'fort.4')
        self.parser.add_argument('-R','--restart',help='name of restart file',type=str,
                        default= 'fort.77')
    def molecules(self):
        self.parser.description = 'Add molecules of a given type fo a given box'
        self.parser.add_argument('-m','--molID',help='molID of molecule to add in fort.4', type=str)
        self.parser.add_argument('-bR','--boxRemove',help='number of box to remove molecules from',type=str)
        self.parser.add_argument('-bA','--boxAdd',help='number of box to add molecules to')
        self.parser.add_argument('-n','--nAdd',help='number of molecules to add. If ==all, all molecules in box',type=int)

class Structure:
    def __init__(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('-f', '--file',help='input file for analysis',type=str)
        self.parser.add_argument('-b','--bins',help='bin size for analysis',type=float)
        self.parser.add_argument('-abc','--vectors',help='unit cell vectors',type=float,nargs='+')
        self.parser.add_argument('-B','--bead',help='bead to analyze', type = str, default=['COM'])
    def analysis(self):
        self.parser.add_argument('-ref','--reference',help='reference density [for dGmap (molec/nm**3)]'
                                                           ' [for S, Kref]',type=float)
        self.parser.add_argument('-n','--numFrames',help='total number of frames',type=int)
        self.parser.add_argument('-T','--Temp',help='Temperature [ K ]',type=float)
    def parse_args(self):
        return self.parser.parse_args()

import argparse, os
