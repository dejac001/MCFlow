from writeXvY import IdealGasAds
class kH(IdealGasAds):
    def __init__(self, **kwargs):
        self.N={};self.P={};self.rho={};self.gen_data={}
        self.files = ['N-data.db','rho-data.db', 'P-data.db','general-data.db']
        self.variables = [self.N, self.P, self.rho, self.gen_data]
        if kwargs:
            assert kwargs['units'], 'Units must be defined for isotherm'
            assert kwargs['box'], 'Box needed for number density in kH isotherm'
            assert kwargs['Temp'], 'Temperature needed for ideal gas law'
            assert (kwargs['henry'] and
            len(kwargs['henry']) == 2), 'kH(mean) and kH(95% conf) needed'
            assert kwargs['mol'], 'Mol needed for x-axis'
            self.kH_mean, self.kH_95conf = kwargs['henry']
            self.mol = kwargs['mol']
            if 'box' in kwargs['box']:
                self.box = kwargs['box']
            else:
                self.box = 'box%s'%kwargs['box']
        self.xlabel = ['kH(mol%s)' % self.mol, 'dkH']

    def getX(self):
        # pressure info ----------------------
        numIndep = self.gen_data[self.feed][self.run]['numIndep']
        try:
            rho_mean, rho_stdev = (self.rho[self.feed][self.run][self.mol][self.box]['mean'],
                                    self.rho[self.feed][self.run][self.mol][self.box]['stdev'])
            # convert number density to pressure [kPa]
            # (molec/nm**3)*(mol/molec)*(nm**3*kPa/(mol*K))*K = kPa
            p_mean = rho_mean/N_av*R['nm**3*kPa/(mol*K)']*args['Temp']
            p_stdev = rho_stdev/N_av*R['nm**3*kPa/(mol*K)']*args['Temp']
        except KeyError:
            print('No num dens. data for feed {}, using box pressure'.format(feed))
            print(' - This assumes box is unary')
            pressure = P[self.feed][self.run]['box%s'%args['box']]
            p_mean, p_stdev = pressure['mean'], pressure['stdev']
        p_95conf = calc95conf(p_stdev,numIndep)
        C_95conf = p_mean * self.kH_mean * math.pow(
            math.pow(self.kH_95conf / self.kH_mean, 2) +
            math.pow(p_95conf / p_mean, 2), 0.5)
        return {'mean':p_mean*self.kH_mean,'95conf':C_95conf}


'''
Write isotherm from previously generated databank files
'''
from runAnalyzer import checkRun, calc95conf
from writeXvY import writeAGR
from chem_constants import R, N_av
import math

if __name__ == '__main__':
    from parser import Plot

    my_parser = Plot()
    my_parser.isotherm()
    my_parser.kH()

    args = vars(my_parser.parse_args())

    my_parser = Plot()
    my_parser.axes()
    my_parser.isotherm()

    args = vars(my_parser.parse_args())
    assert args['yaxis'], 'No y axis chosen for plot'

    my_plotter = kH(**args)
    my_plotter.readDBs()

    for feed in args['feeds']:
        # determine if run has been completed
        my_plotter.run = checkRun(args['type'], my_plotter.variables, feed)
        my_plotter.feed = feed
        if args['yaxis'] == 'Q':
            my_plotter.QvX()
        elif args['yaxis'] == 'S':
            my_plotter.SvX()


