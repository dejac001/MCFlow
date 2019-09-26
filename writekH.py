from MCFlow.writeXvY import IdealGasAds, writeAGR
from MCFlow.writeHB import HBvC

class kH(HBvC,IdealGasAds):
    def __init__(self, **kwargs):
        self.N={};self.P={};self.rho={};self.gen_data={};self.K={}; self.X={}; self.dG={}; self.dHmixt={}; self.U = {}
        self.files = ['N-data.json','P-data.json', 'rho-data.json','general-data.json','K-data.json','X-data.json','dG-data.json','dH-mixt-data.json', 'U-data.json']
        self.variables = [self.N, self.P, self.rho, self.gen_data,self.K, self.X,self.dG,self.dHmixt, self.U]
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

    def dHigvX(self):
        U = self.U[self.feed][self.run]
        dH_store = self.dHmixt[self.feed][self.run]
        gen_data = self.gen_data[self.feed][self.run]
        file_name = 'dHig_vapor_to_%s.dat'%self.box
        X = self.getX()
        N = self.N[self.feed][self.run]
        P = self.P[self.feed][self.run]
        N1_tot = sum(N[i][self.box]['mean'] for i in N.keys())
        _rho_ = self.rho[self.feed][self.run]
        rho_total = getRhoTotal(_rho_)
        N1 = {'mean':N1_tot,
              'stdev':math.pow( sum(N[i][self.box]['stdev']**2 for i in N.keys()), 0.5)} # mol
        vapor_box = 'box3'
        N2 = {'mean':
    sum(N[i]['box3']['mean'] for i in N.keys())}
        N2['stdev'] = math.pow(
sum(N[i]['box3']['stdev']**2 for i in N.keys())
,  0.5)
        nIndep = gen_data['numIndep']
#       dH = ((U[self.box]['mean']/N1['mean'] - U[vapor_box]['mean']/N2['mean'])*8.314/1000 - 8.314/1000*self.T  -
#               dH_store['box3-->box2']['mean'])
        H_box1 = U[self.box]['mean']/N1['mean']
        H_box2 = P['box2']['mean']*1000./rho_total['box2']*N_av/R['\AA**3*kPa/(mol*K)']
        H_box3 = U[vapor_box]['mean']/N2['mean'] - P['box3']['mean']*1000./rho_total['box3']*N_av/R['\AA**3*kPa/(mol*K)']
        dH = (H_box1 - H_box2  - H_box3)*R['kJ/(mol*K)']
        ddH = math.pow(
            (1/N1['mean'])**2*U['box1']['stdev']**2 +
            (-1*U['box1']['mean']/N1['mean']**2)**2*N1['stdev']**2  +
            (1/N2['mean'])**2*U[vapor_box]['stdev']**2 +
            (-1*U[vapor_box]['mean']/N2['mean']**2)**2*N2['stdev']**2
            + dH_store['box3-->box2']['stdev']**2
            ,0.5)*8.314/1000
        dH_mean, dH_stdev = dH, ddH
        file_description = 'dHig     %s    ddHig     %s'%(self.xlabel[0],
                                                       self.xlabel[1])
        if '95conf' not in X.keys():
            X['95conf'] = calc95conf(X['stdev'], nIndep)
        writeAGR([X['mean']],[dH_mean],
                [X['95conf']], [calc95conf(dH_stdev, nIndep)],
                [self.feed], file_name, file_description)

    def dHvX(self):
        P = self.P[self.feed][self.run]
        U = self.U[self.feed][self.run]
        _rho_ = self.rho[self.feed][self.run]
        rho_total = getRhoTotal(_rho_)
        dH = self.dHmixt[self.feed][self.run]
        N = self.N[self.feed][self.run]
        N1_tot = sum(N[i][self.box]['mean'] for i in N.keys())
        N1 = {'mean':N1_tot,
              'stdev':math.pow( sum(N[i][self.box]['stdev']**2 for i in N.keys()), 0.5)} # mol
        gen_data = self.gen_data[self.feed][self.run]
        file_description = 'Q(%s)     %s    dQ     %s'%(self.units, self.xlabel[0],
                                                       self.xlabel[1])
        if self.boxes:
            boxFrom, boxTo = self.boxes
#       else:
#           boxFrom = self.box
#           boxTo = 'box1'
        file_name = 'dH_%s_mixture.dat'%boxTo
        assert 'box' in boxFrom, 'Wrong notation for box'
        transfer = boxFrom + '-->' + boxTo
        transfer2 = 'box2-->box1'
        X = self.getX()
        nIndep = gen_data['numIndep']
        H_zeo = (U['box1']['mean']/N1_tot + 0.1*1000./rho_total['box1']*N_av/R['\AA**3*kPa/(mol*K)']) * R['kJ/(mol*K)']
        dH_mean, dH_stdev = (dH[transfer]['mean'] + dH[transfer2]['mean'] - H_zeo,
                                   math.pow( dH[transfer]['stdev']**2 + 4*dH[transfer2]['stdev']**2 , 0.5))
        if '95conf' not in X.keys():
            X['95conf'] = calc95conf(X['stdev'], nIndep)
        writeAGR([X['mean']],[dH_mean],
                [X['95conf']], [calc95conf(dH_stdev, nIndep)],
                [self.feed], file_name, file_description)

    def getX(self):
        # pressure info ----------------------
        numIndep = self.gen_data[self.feed][self.run]['numIndep']
#       try:
            # convert number density to pressure [kPa]
            # (molec/nm**3)*(mol/molec)*(nm**3*kPa/(mol*K))*K = kPa
        my_box = self.findVapBox(self.rho[self.feed][self.run], self.mol)
        rho = self.rho[self.feed][self.run][self.mol][my_box]
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
from runAnalyzer import checkRun, getRhoTotal
from statistics import calc95conf
from chem_constants import R, N_av
import math

if __name__ == '__main__':
    from analysis_parsers import Plot

    my_parser = Plot()
    my_parser.isotherm()
    my_parser.kH()
    my_parser.parser.add_argument('-a','--angle',help='angle for db criteria',type=str)
    my_parser.parser.add_argument('-d','--dist',help='distance for db criteria',type=str)

    args = vars(my_parser.parse_args())
    assert args['yaxis'], 'No y axis chosen for plot'

    my_plotter = kH(**args)
    my_plotter.read_json()

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
        elif args['yaxis'] == 'dH':
            my_plotter.dHvX()
        elif args['yaxis'] == 'dHig':
            my_plotter.dHigvX()
        elif args['yaxis'] == 'HB':
            my_plotter.readHB()
            my_plotter.HB_write()
