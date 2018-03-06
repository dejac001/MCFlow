from MCFlow.writeXvY import IdealGasAds
from MCFlow.writeHB import HBvC

class kH(HBvC,IdealGasAds):
    def __init__(self, **kwargs):
        self.N={};self.P={};self.rho={};self.gen_data={};self.K={}; self.X={}; self.dG={}
        self.files = ['N-data.db','P-data.db', 'rho-data.db','general-data.db','K-data.db','X-data.db','dG-data.db']
        self.variables = [self.N, self.P, self.rho, self.gen_data,self.K, self.X,self.dG]
        self.xlabel = ['kH', 'dkH']
        if kwargs['feeds']:
            self.feeds = kwargs['feeds']
        if 'angle' in kwargs.keys() and kwargs['angle']:
            self.angle = kwargs['angle']
        if 'dist' in kwargs.keys() and kwargs['dist']:
            self.dist = kwargs['dist']
        if kwargs:
            self.boxes = kwargs['boxes']
            self.units = kwargs['units']
            assert kwargs['box'], 'Box needed for number density in kH isotherm'
            if 'box' in kwargs['box']:
                self.box = kwargs['box']
            else:
                self.box = 'box%s'%kwargs['box']
            assert kwargs['Temp'], 'Temperature needed for ideal gas law'
            self.T = kwargs['Temp']
            assert (kwargs['henry'] and
                        len(kwargs['henry']) == 2), 'kH(mean) and kH(95% conf) needed'
            self.kH_mean, self.kH_95conf = kwargs['henry']
            assert kwargs['mol'], 'Mol needed for x-axis'
            self.mol = kwargs['mol']
            self.indep = kwargs['indep']
            self.feeds = kwargs['feeds']
            self.path = kwargs['path']
#           assert kwargs['density'], 'Infinite-dilution liquid density needed for kH iso'
            if 'density' in kwargs.keys():
                self.density = kwargs['density']
            if 'film' in kwargs.keys():
                self.film = kwargs['film']
            else:
                self.film = False

    def getX(self):
        # pressure info ----------------------
        numIndep = self.gen_data[self.feed][self.run]['numIndep']
#       try:
            # convert number density to pressure [kPa]
            # (molec/nm**3)*(mol/molec)*(nm**3*kPa/(mol*K))*K = kPa
        rho = self.rho[self.feed][self.run][self.mol]['box3']
        p_mean = rho['mean']/N_av*R['nm**3*kPa/(mol*K)']*self.T
        p_stdev = rho['stdev']/N_av*R['nm**3*kPa/(mol*K)']*self.T
        p_raw = [i/N_av*R['nm**3*kPa/(mol*K)']*self.T for i in rho['raw']]
#       except :
#           pressure = self.P[self.feed][self.run][self.box]
#           print('No num dens. data for feed {}, using box pressure'.format(feed))
#           print(' - This assumes box is unary')
#           p_mean, p_stdev = pressure['mean'], pressure['stdev']
        p_95conf = calc95conf(p_stdev,numIndep)
        C_95conf = p_mean * self.kH_mean * math.pow(
            math.pow(self.kH_95conf / self.kH_mean, 2) +
            math.pow(p_95conf / p_mean, 2), 0.5)
        return {'mean':p_mean*self.kH_mean,'95conf':C_95conf,
                'raw':[i*self.kH_mean for i in p_raw]}

'''
Write isotherm from previously generated databank files
'''
from runAnalyzer import checkRun, calc95conf
from chem_constants import R, N_av
import math

if __name__ == '__main__':
    from parser import Plot

    my_parser = Plot()
    my_parser.isotherm()
    my_parser.kH()
    my_parser.parser.add_argument('-a','--angle',help='angle for db criteria',type=str)
    my_parser.parser.add_argument('-d','--dist',help='distance for db criteria',type=str)

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
        elif args['yaxis'] == 'X':
            my_plotter.XvX()
        elif args['yaxis'] == 'dG':
            my_plotter.dGvX()
        elif args['yaxis'] == 'HB':
            my_plotter.readHB()
            my_plotter.HB_write()
