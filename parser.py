class Main:
    def __init__(self):

        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("-v", "--verbosity", action="count", default=0)
        self.parser.add_argument('-t','--type',
                             help='type of run (equil- or prod-)',
                             type = str, default='equil-')
        self.parser.add_argument('-p','--path',
                                 help='main path to feeds (see --feeds)',
                                 type=str,default=os.getcwd())
        self.parser.add_argument('-f','--feeds',
                                 help ='state point (ex: temperature or pressure);'
                                             ' must be in ls ${path}',
                             nargs = '+', type=str)
        self.parser.add_argument('-i','--indep',
                        help='RANGE of independent simulations',
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
        return self.args

class Plot(Main):
    def __init__(self):
        Main.__init__(self)
        self.parser.description = 'Plot results'
        # TODO: reformat mol pass in so pass in by str name (i.e. WATER or 15PDO)
        self.parser.add_argument('-m','--mol',help='Molecule to analyze', type=str)
        self.parser.add_argument('-TK','--Temp',help='Temperature in Kelvin',
                                  type=float)
        self.parser.add_argument('-b','--box',help='box to analyze', type=str)
    def axes(self):
        self.parser.add_argument('-x','--xaxis', help='x axis of plot',
                                 choices = ['C','Pig'])
        self.parser.add_argument('-y','--yaxis',help='x axis of plot',
                                 choices = ['Q','dG','S'])
    def isotherm(self):
        self.parser.add_argument('-u','--units',help='choices for units on plot',
                                 choices=['molec/uc','g/g','mol/kg',
                                          '(mol/mol)/(mol/mol)','(mol/mol)/(kPa/kPa)'])
    def kH(self):
        self.parser.add_argument('-kH','--henry',help='Henry constant (g/mL/kPa) '
                                                       'fmt: (mean, 95%%conf.) ',
                                  type=float, nargs= '+')




class Results(Plot):
    def __init__(self):
        Plot.__init__(self)
        Main.other(self)
        self.parser.description = 'Obtain results for simulations'

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
        self.parser.add_argument('-b','--box',help='number of box to add molecules to',type=str)
        self.parser.add_argument('-n','--nAdd',help='number of molecules to add',type=int)

import argparse, os
