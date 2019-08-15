
from MCFlow.writeXvY import *
from MCFlow.writekH import kH

class mykH(kH):
    def __init__(self, **kwargs):
        kH.__init__(self, **kwargs)
    def QvX(self):
        '''
        :var X: either solution concentration (g/mL) or pressure (kPa)
        '''
        assert self.units, 'Units must be defined for isotherm'
        N = self.N[self.feed][self.run]
        gen_data = self.gen_data[self.feed][self.run]
        file_description = '%s    Q(%s)    %s     dQ'%(self.xlabel[0], self.units,
                                                       self.xlabel[1])
        X = self.getX()
        nIndep = gen_data['numIndep']
        zeolite_mass_g = gen_data['zeolite']['mass (g)']
        for mol in N.keys():
            for channel, values in N[mol].items():
                if None in N.keys():
                    N['None'] = N.pop(None)
                file_name  = 'Q%s-%s-%s-w%sin-zeo-vs-%s.dat'%(mol,channel,self.units,''.join(sorted(N.keys())),
                                                               self.xlabel[0][:self.xlabel[0].find('(')])
                if self.units == 'molec/uc':
                    qfactor = 1/gen_data['zeolite']['unit cells']
                elif self.units == 'g/g':
                    qfactor = ((gen_data['molecular weight'][mol]/N_av) /
                                zeolite_mass_g )
                elif self.units == 'mol/kg':
                    qfactor = gen_data['zeolite'][' mol/kg / 1 mlcl adsorbed']
                elif self.units == 'None':
                    qfactor = 1.
                if (0 in self.indep) and (len(self.indep) == 1):
                    Q_vals = [i*qfactor for i in values['raw']]
                    Q_stdev = [0. for i in Q_vals]
                    writeAGR( X['raw'], Q_vals, None,
                                None, ['%s/%i'%(self.feed,j) for j in
                                range(1, nIndep+1)], file_name, file_description)
                else:
                    if 'all trans mean' in values.keys():
                        Q_mean, Q_stdev = values['all trans mean']*qfactor, values['all trans stdev']*qfactor
                    else:
                        Q_mean, Q_stdev = values['mean']*qfactor, values['stdev']*qfactor
#                   Q_mean, Q_stdev = values['mean']*qfactor, values['stdev']*qfactor
                    if '95conf' in X.keys():
                        dX = X['95conf']
                    else:
                        dX = calc95conf(X['stdev'], nIndep)
                    writeAGR([X['mean']],[Q_mean],
                             [dX], [calc95conf(Q_stdev, nIndep)],
                             [self.feed], file_name, file_description)

class myLiqAds(mykH, LiqAds):

    def __init__(self, **kwargs):
        LiqAds.__init__(self, **kwargs)

    def getX(self):
        return LiqAds.getX(self)

    

if __name__ == '__main__':
    # TODO: find way to have conditional arguments based off of what xaxis and yaxis are
    # (in other words, we have an excess of variables as arguments)
    my_parser = Plot()
    my_parser.axes()
    my_parser.isotherm()
    my_parser.parser.add_argument('-N-data','--NDATA',help='different file to replace N-data.json with',
                                  type=str)

    args = vars(my_parser.parse_args())
    assert args['yaxis'], 'No y axis chosen for plot'
    assert args['xaxis'], 'No x axis chosen for plot'
    assert args['NDATA'], 'Specialed N-data.jsonneeded to be provided'


    if args['xaxis'] == 'C':
        my_plotter = myLiqAds(**args)
    elif args['xaxis'] == 'Pig':
        # TODO: add alternative way to calculate P w/ P-data.json & using mole fraction in box
        my_plotter = IdealGasAds(**args)
    elif args['xaxis'] == 'rho':
        my_plotter = RhoBoxAds(**args)
    elif args['xaxis'] == 'Pbox':
        my_plotter = GasBoxAds(**args)
    elif args['xaxis'] == 'Pi':
        my_plotter = GasMolAds(**args)
    elif args['xaxis'] == 'Q':
        my_plotter = LoadAds(**args)
    elif args['xaxis'] == 'x':
        my_plotter = MoleFrac(**args)
    elif args['xaxis'] == 'kH':
        my_plotter = mykH(**args)
    old_index = my_plotter.files.index('N-data.json')
    my_plotter.files[old_index] = args['NDATA']
    my_plotter.read_json()

    for feed in args['feeds']:
        # determine if run has been completed
        print(feed, 'before')
        my_plotter.run = checkRun(args['type'], my_plotter.variables, feed)
        print(feed, my_plotter.run)
        my_plotter.feed = feed
        if args['yaxis'] == 'Q':
            my_plotter.QvX()
        elif args['yaxis'] == 'S':
            my_plotter.SvX()
        elif args['yaxis'] == 'dG':
            my_plotter.dGvX()
        elif args['yaxis'] == 'R':
            my_plotter.RvX()
        elif args['yaxis'] == 'X':
            my_plotter.XvX()
        elif args['yaxis'] == 'dH':
            my_plotter.dHvX()
        elif args['yaxis'] == 'dU':
            my_plotter.dUvX()
        elif args['yaxis'] == 'Pigxy':
            my_plotter.Pig_xy()
        elif args['yaxis'] == 'Txy':
            my_plotter.Txy()
