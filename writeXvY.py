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
            if '@' not in description: description = '# ' + description
            f.write('%s\n'%description)
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

def calculateMoleFrac(N, mol, box):
    N_tot = sum(N[i][box]['mean'] for i in N.keys())
    mean = N[mol][box]['mean']/N_tot
    stdev2 = math.pow((1/N_tot - N[mol][box]['mean']/N_tot**2),2)*math.pow(N[mol][box]['mean'],2)
    for molOther in [i for i in N.keys() if i != mol]:
        stdev2 += math.pow(- N[mol][box]['mean']/N_tot**2,2)*math.pow(N[molOther][box]['mean'],2)
    stdev = math.pow(stdev2,0.5)
    return mean, stdev

def calculateIGMolFrac(rho, mol, box, Ptot, T):
    '''

    :param rho: molec/uc
    :param mol:
    :param box:
    :param Ptot: kPa
    :return:
    '''
    print('- calculating mole fraction in %s assuming ideal gas'%(box))
    if 'box' not in box: box = 'box' + box
    number_density = rho[mol][box]
    temperature = T
    p_mean = number_density['mean']/N_av*R['nm**3*kPa/(mol*K)']*temperature
    p_stdev = number_density['stdev']/N_av*R['nm**3*kPa/(mol*K)']*temperature
    y_mean = p_mean/Ptot['mean']
    y_stdev = y_mean*math.sqrt( pow(p_stdev/p_mean,2) +
                                pow(Ptot['stdev']/Ptot['mean'],2))
    return y_mean,y_stdev


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
        self.N = {}; self.rho = {}; self.gen_data = {}; self.dG = {}
        self.files = ['N-data.db','rho-data.db','general-data.db','dG-data.db']
        self.variables = [self.N, self.rho, self.gen_data, self.dG]
        if kwargs:
            assert kwargs['mol'], 'Mol needed for number density to calculate P assuming I.G.'
            assert kwargs['Temp'], 'Temperature needed for number density to calculate P assuming I.G.'
            self.xlabel = ['Pig-mol%s_%4.1f'%(kwargs['mol'],kwargs['Temp']),'dP']
            self.T = kwargs['Temp']
            self.mol = kwargs['mol']
            self.indep = kwargs['indep']
            self.units = kwargs['units']
            self.feeds = kwargs['feeds']
            self.path = kwargs['path']
            self.box = kwargs['box']
            self.boxes = kwargs['boxes']

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

    def dGvX(self):
        if (0 in self.indep) and (len(self.indep) == 1): raise NotImplemented
        assert len(self.boxes) == 2, 'Too many boxes'
        file_description = '%s    dG(kJ/mol)    %s     dG'%(self.xlabel[0],self.xlabel[1])
        N = self.N[self.feed][self.run]
        dG = self.dG[self.feed][self.run]
        nIndep = self.gen_data[self.feed][self.run]['numIndep']
        X = self.getX()
        boxFrom, boxTo = self.boxes
        assert 'box' in boxFrom, 'Wrong notation for box'
        solutes = sorted([i for i in N.keys()
                                    if ((N[i][boxTo]['mean'] > 1e-06)
                          and (N[i][boxTo]['mean'] < 300))      ])
        for mol in solutes:
            file_name = 'dG-mol%s_vs_%s.dat'%(mol, self.xlabel[0])
            dG_mean, dG_stdev = (dG[mol]['--'.join([boxFrom, boxTo])]['mean'],
                                dG[mol]['--'.join([boxFrom, boxTo])]['stdev'])
            writeAGR([X['mean']],[dG_mean],
                     [calc95conf(X['stdev'], nIndep)], [calc95conf(dG_stdev, nIndep)],
                     [self.feed], file_name, file_description)

    def getX(self):
        vapor_box = self.findVapBox( self.rho[self.feed][self.run], self.mol)
        # todo: need to choose molecule for x axis here
        number_density = self.rho[self.feed][self.run][self.mol][vapor_box]
        temperature = self.T
        p_mean = number_density['mean']/N_av*R['nm**3*kPa/(mol*K)']*temperature
        p_stdev = number_density['stdev']/N_av*R['nm**3*kPa/(mol*K)']*temperature
        p_raw = [i/N_av*R['nm**3*kPa/(mol*K)']*temperature for i in number_density['raw']]
        return {'mean':p_mean,'stdev':p_stdev, 'raw':p_raw}

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
        print('For mol%s, vapor box determined to be %s'%(mol, my_box))
        return my_box
    
    def QvX(self):
        '''
        :var X: either solution concentration (g/mL) or pressure (kPa)
        '''   
        N = self.N[self.feed][self.run]
        gen_data = self.gen_data[self.feed][self.run]
        file_description = '%s    Q(%s)    %s     dQ'%(self.xlabel[0], self.units,
                                                       self.xlabel[1])
        X = self.getX()
        nIndep = gen_data['numIndep']
        mols_adsorbed = ['1','8'] #self.getMolAds(N)
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
            if (0 in self.indep) and (len(self.indep) == 1):
                Q_vals = [i*qfactor for i in N[mol]['box1']['raw']]
                Q_stdev = [N[mol]['box1']['stdev']*qfactor for i in Q_vals]
                writeAGR( X['raw'], Q_vals, None,
                            None, ['%s/%i'%(feed,j) for j in 
                            range(1, nIndep+1)], file_name, file_description)
            else:
                Q_mean, Q_stdev = (N[mol]['box1']['mean']*qfactor, N[mol]['box1']['stdev']*qfactor)
                writeAGR([X['mean']],[Q_mean],
                         [calc95conf(X['stdev'], nIndep)], [calc95conf(Q_stdev, nIndep)],
                         [feed], file_name, file_description)

    def RvX(self):
        assert self.mol, 'Molecule needed to calculate recovery'
        assert self.box, 'Box needed to calculate recovery'
        if (0 in self.indep) and (len(self.indep) == 1): raise NotImplemented
        N = self.N[self.feed][self.run]
        # N: feed, run, mol, box, mean/stdev
        X = self.getX()
        nIndep = self.gen_data[self.feed][self.run]['numIndep']
        file_name = 'R_mol%s_box%s-vs-%s.dat'%(self.mol, self.box, self.xlabel[0])
        file_description = '%s    Recovery  %s     dR'%(self.xlabel[0], self.xlabel[1])
        Ni_tot = sum(N[self.mol][i]['mean'] for i in N[self.mol].keys())
        R_mean = N[self.mol]['box%s'%self.box]['mean']/Ni_tot*100
        R_stdev = pow(pow(1/Ni_tot - N[self.mol]['box%s'%self.box]['mean']/pow(Ni_tot,2),2)*pow(N[self.mol]['box%s'%self.box]['stdev'],2),0.5)*100
        writeAGR([X['mean']],[R_mean],
                         [calc95conf(X['stdev'], nIndep)], [calc95conf(R_stdev, nIndep)],
                         [feed], file_name, file_description)

    def XvX(self):
        if (0 in self.indep) and (len(self.indep) == 1): raise NotImplemented
        N = self.N[self.feed][self.run]
        X = self.getX()
        nIndep = self.gen_data[self.feed][self.run]['numIndep']
        file_description = '%s    x (mol/mol)    %s     dx'%(self.xlabel[0],  self.xlabel[1])
        for box in N['1'].keys():
            for mol1 in N.keys():
                if sum(N[mol1][ibox]['mean'] for ibox in N[mol1].keys()) > 0.:
                    file_name = 'X_mol%s_%s-vs-%s.dat'%(mol1, box,self.xlabel[0])
                    x_mean, x_stdev = calculateMoleFrac(N, mol1, box)
                    if ('P-box' in self.xlabel[0]) and (x_mean < 0.75):
                        assert self.T, 'Temperature needed for vapor p'
                        rho = self.rho[self.feed][self.run]
                        x_mean, x_stdev = calculateIGMolFrac(rho, mol1, box, X, self.T)
                    writeAGR([X['mean']],[x_mean],
                         [calc95conf(X['stdev'], nIndep)], [calc95conf(x_stdev, nIndep)],
                         [feed], file_name, file_description)
    
    def SvX(self):
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
        if (0 in self.indep) and (len(self.indep) == 1): raise NotImplemented
        N = self.N[self.feed][self.run]
        rho = self.rho[self.feed][self.run]
        X = self.getX()
        file_description = '%s    S(%s)    %s     dS'%(self.xlabel[0], self.units, self.xlabel[1])
        mols_adsorbed = self.getMolAds(N)
        mols_adsorbed.reverse() # start from higher number molecules (i.e. solutes)
        nIndep = self.gen_data[self.feed][self.run]['numIndep']
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

class RhoBoxAds(IdealGasAds):
    def __init__(self, **kwargs):

        # TODO: make general db reader that all files can do
        self.N = {}; self.rho = {}; self.gen_data = {}
        self.files = ['N-data.db','rho-data.db','general-data.db']
        self.variables = [self.N, self.rho, self.gen_data]
        if kwargs:
            assert kwargs['mol'], 'Mol needed for number density'
            assert kwargs['box'], 'Box needed for number density'
            self.xlabel = ['rho-mol%s_box%s'%(kwargs['mol'],kwargs['box']),'drho']
            self.T = kwargs['Temp']
            self.mol = kwargs['mol']
            self.units = kwargs['units']
            self.feeds = kwargs['feeds']
            self.path = kwargs['path']
            self.indep = kwargs['indep']
            self.box = kwargs['box']
    def getX(self):
        return self.rho[self.feed][self.run][self.mol]['box%s'%self.box]

class GasBoxAds(IdealGasAds):
    def __init__(self, **kwargs):
        self.N = {}; self.P = {}; self.gen_data = {}; self.rho = {}; self.dG = {}
        self.files = ['N-data.db','P-data.db','general-data.db','rho-data.db', 'dG-data.db']
        self.variables = [self.N, self.P, self.gen_data, self.rho, self.dG]
        if kwargs:
            assert kwargs['box'], 'Box needed to calculate pressure (box)'
            self.xlabel = ['P-box%s(kPa)'%kwargs['box'],'dP']
            self.vapor_box = kwargs['box']
            self.box = kwargs['box']
            self.indep = kwargs['indep']
            self.mol = kwargs['mol']
            self.units = kwargs['units']
            self.feeds = kwargs['feeds']
            self.path = kwargs['path']
            self.T = kwargs['Temp']
            self.boxes = kwargs['boxes']
    def getX(self):
        return self.P[self.feed][self.run]['box%s'%self.vapor_box]

class LoadAds(IdealGasAds):
    def __init__(self, **kwargs):
        self.N = {}; self.P = {}; self.gen_data = {}; self.U = {}; self.boxlx = {}; self.dG = {}
        self.files = ['N-data.db','P-data.db','general-data.db','U-data.db','boxlx-data.db','dG-data.db']
        self.variables = [self.N, self.P, self.gen_data, self.U, self.boxlx, self.dG]
        if kwargs:
            # assert ('box' in kwargs['box']
            #         or 'box' in kwargs['boxes'][0]), 'Box needed for enthalpy of adsorption from'
            self.xlabel = ['%s'%kwargs['xaxis'],'d%s'%kwargs['xaxis']]
            self.box = kwargs['box']
            self.feeds = kwargs['feeds']
            self.units = kwargs['units']
            self.path = kwargs['path']
            self.mol = kwargs['mol']
            self.zeoVolume = kwargs['ZeoVolume']
            self.indep = kwargs['indep']
            self.T = kwargs['Temp']
            self.boxes = kwargs['boxes']
    def getX(self):
        if (0 in self.indep) and (len(self.indep) == 1): raise NotImplemented
        N = self.N[self.feed][self.run]
        gen_data = self.gen_data[self.feed][self.run]
        if self.units == 'molec/uc':
            qfactor = 1/gen_data['zeolite']['unit cells']
        elif self.units == 'g/g':
            raise AttributeError
        elif self.units == 'mol/kg':
            qfactor = gen_data['zeolite'][' mol/kg / 1 mlcl adsorbed']
        Q_mean, Q_stdev = (sum(N[i]['box1']['mean'] for i in N.keys())*qfactor,
                           math.pow( sum(N[i]['box1']['stdev']**2 for i in N.keys()), 0.5)*qfactor)
        return {'mean':Q_mean, 'stdev':Q_stdev}

# def dHigvX(self):
#         U = self.U[self.feed][self.run]
#         N = self.N[self.feed][self.run]
#         P = self.P[self.feed][self.run]
#         gen_data = self.gen_data[self.feed][self.run]
#         box_data = self.boxlx[self.feed][self.run][self.box]
#         file_description = 'Q(%s)     %s    dQ     %s'%(self.units, self.xlabel[0],
#                                                        self.xlabel[1])
#         file_name = 'dH_%s_mol%s.dat'%(self.box,self.mol)
#         X = self.getX()
#         # Vz = {'mean':self.zeoVolume*10**(-24), 'stdev':0.} # cm**3
#         # Vv = {'mean':box_data['mean']*10**(-24), 'stdev':box_data['stdev']*10**(-24)} # cm**3
#         Nz_tot = sum(N[i]['box1']['mean'] for i in N.keys())
#         # Nz = {'mean':Nz_tot/N_av,
#         #       'stdev':math.pow( sum(N[i]['box1']['stdev']**2 for i in N.keys()), 0.5)/N_av} # mol
#         Nz = {'mean':Nz_tot,
#               'stdev':math.pow( sum(N[i]['box1']['stdev']**2 for i in N.keys()), 0.5)} # mol
#         Nv_tot = sum(N[i][self.box]['mean'] for i in N.keys())
#         # Nv = {'mean':Nv_tot/N_av,
#         #       'stdev':math.pow( sum(N[i][self.box]['stdev']**2 for i in N.keys()), 0.5)/N_av} # mol
#         Nv = {'mean':Nv_tot,
#               'stdev':math.pow( sum(N[i][self.box]['stdev']**2 for i in N.keys()), 0.5)} # mol
#         nIndep = gen_data['numIndep']
#         R = 8314. #cm**3*kPa / (mol*K)
#         # dH = ( U['box1']['mean'] -
#         #        U[self.box]['mean']) + P[self.box]['mean']*(Vz['mean']/Nz['mean']
#         #                                             - Vv['mean']/Nv['mean'])/R # units of Kelvin
#         # ddH = math.pow(
#         #     U['box1']['stdev']**2 +
#         #     U[self.box]['stdev']**2 +
#         #     (1/R*(Vz['mean']/Nz['mean'] - Vv['mean']/Nv['mean']))**2*P[self.box]['stdev']**2 +
#         #     (P[self.box]['mean']/(R*Nz['mean']))**2*Vz['stdev']**2 +
#         #     (-1*P[self.box]['mean']*Vz['mean']/(Nz['mean']**2*R))**2*Nz['stdev']**2 +
#         #     (-1*P[self.box]['mean']/(Nv['mean']*R))**2*Vv['stdev']**2 +
#         #     (P[self.box]['mean']/(R*Nv['mean']**2))**2*Nv['stdev']**2
#         # , 0.5)
#         dH = (U['box1']['mean']/Nz_tot - U[self.box]['mean']/Nv_tot)*8.314/1000 - 8.314/1000*self.T
#         ddH = math.pow(
#             (1/Nz_tot)**2*U['box1']['stdev']**2 +
#             (-1*U['box1']['mean']/Nz['mean']**2)**2*Nz['stdev']**2 +
#             (1/Nv_tot)**2*U[self.box]['stdev']**2 +
#             (-1*U[self.box]['mean']/Nv['mean']**2)**2*Nv['stdev']**2
#             ,0.5)*8.314/1000
#         # dH_mean, dH_stdev = (8.314/1000*dH, 8.314/1000*ddH)
#         dH_mean, dH_stdev = dH, ddH
#         writeAGR([X['mean']],[dH_mean],
#                 [calc95conf(X['stdev'], nIndep)], [calc95conf(dH_stdev, nIndep)],
#                 [self.feed], file_name, file_description)

    def dUvX(self):
        U = self.U[self.feed][self.run]
        N = self.N[self.feed][self.run]
        gen_data = self.gen_data[self.feed][self.run]
        file_description = 'Q(%s)     %s    dQ     %s'%(self.units, self.xlabel[0],
                                                       self.xlabel[1])
        file_name = 'dU_%s_mixture.dat'%self.box
        X = self.getX()
        N1_tot = sum(N[i]['box1']['mean'] for i in N.keys())
        N1 = {'mean':N1_tot,
              'stdev':math.pow( sum(N[i]['box1']['stdev']**2 for i in N.keys()), 0.5)} # mol
        N2_tot = sum(N[i][self.box]['mean'] for i in N.keys())
        N2 = {'mean':N2_tot,
              'stdev':math.pow( sum(N[i][self.box]['stdev']**2 for i in N.keys()), 0.5)} # mol
        nIndep = gen_data['numIndep']
        dU = (U['box1']['mean']/N1_tot - U[self.box]['mean']/N2_tot)*8.314/1000
        ddU = math.pow(
            (1/N1_tot)**2*U['box1']['stdev']**2 +
            (-1*U['box1']['mean']/N1['mean']**2)**2*N1['stdev']**2 +
            (1/N2_tot)**2*U[self.box]['stdev']**2 +
            (-1*U[self.box]['mean']/N2['mean']**2)**2*N2['stdev']**2
            ,0.5)*8.314/1000
        dU_mean, dU_stdev = dU, ddU
        writeAGR([X['mean']],[dU_mean],
                [calc95conf(X['stdev'], nIndep)], [calc95conf(dU_stdev, nIndep)],
                [self.feed], file_name, file_description)


class GasMolAds(IdealGasAds):
    def __init__(self, **kwargs):
        self.N = {}; self.P = {}; self.gen_data = {}; self.rho = {}; self.dG = {}
        self.files = ['N-data.db','P-data.db','general-data.db','rho-data.db', 'dG-data.db']
        self.variables = [self.N, self.P, self.gen_data, self.rho, self.dG]
        if kwargs:
            assert kwargs['mol'], 'Mol needed to calculate partial pressure'
            self.box = kwargs['box']
            self.xlabel = ['Pi_mol%s(kPa)'%kwargs['mol'],'dP']
            self.mol = kwargs['mol']
            self.units = kwargs['units']
            self.feeds = kwargs['feeds']
            self.indep = kwargs['indep']
            self.path = kwargs['path']
            self.boxes = kwargs['boxes']
    def getX(self):
        # N: feed, run, mol, box, mean/stdev
        # P: feed, run, box, mean/stdev
        vapor_box = self.findVapBox( self.rho[self.feed][self.run], self.mol)
        # todo: need to choose molecule for x axis here
        P, N = self.P[self.feed][self.run], self.N[self.feed][self.run]
        N_tot = sum(N[i][vapor_box]['mean'] for i in N.keys())
        N_tot_raw = np.zeros(len(N['1'][vapor_box]['raw']))
        for mol in N.keys():
            for i, indep in enumerate(N[mol][vapor_box]['raw']):
                N_tot_raw[i] = N_tot_raw[i] + np.array(indep)
        y_raw = [Ni/N for Ni, N in zip(N[self.mol][vapor_box]['raw'], N_tot_raw ) ]
        y_mean = N[self.mol][vapor_box]['mean']/N_tot
        mean = P[vapor_box]['mean']*y_mean
        raw = [Ptot*yi for Ptot, yi in zip(P[vapor_box]['raw'], y_raw)]
        stdev = pow(y_mean,2)*pow(P[vapor_box]['stdev'],2) + pow(P[vapor_box]['mean']/N_tot -
                   P[vapor_box]['mean']*N[self.mol][vapor_box]['mean']/pow(N_tot,2),2)*pow(N[self.mol][vapor_box]['stdev'],2)
        for my_mol in [i for i in N.keys() if i != self.mol]:
            if N[my_mol][vapor_box]['mean'] > 0.:
                stdev += pow(-P[vapor_box]['mean']*N[self.mol][vapor_box]['mean']/pow(N_tot,2),2)*pow(N[my_mol][vapor_box]['stdev'],2)
        stdev = pow(stdev, 0.5)
        return {'mean':mean, 'stdev':stdev, 'raw':raw}


class LiqAds(IdealGasAds):
    def __init__(self, **kwargs):
        self.dG = {}; self.C = {}; self.N = {}; self.rho = {}; self.gen_data = {}
        self.files = ['dG-data.db','Conc-data.db','N-data.db','rho-data.db','general-data.db']
        self.variables = [self.dG, self.C, self.N, self.rho, self.gen_data]
        if kwargs:
            self.xlabel = ['C(g/mL)','dC']
            self.mol = kwargs['mol']
            self.units = kwargs['units']
            self.indep = kwargs['indep']
            self.feeds = kwargs['feeds']
            self.path = kwargs['path']
            self.box = kwargs['box']
            self.boxes = kwargs['boxes']
    


    def getX(self):
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

class MoleFrac(LiqAds):
    def __init__(self, **kwargs):
        self.dG = {}; self.C = {}; self.N = {}; self.rho = {}; self.gen_data = {}
        self.files = ['dG-data.db','Conc-data.db','N-data.db','rho-data.db','general-data.db']
        self.variables = [self.dG, self.C, self.N, self.rho, self.gen_data]
        if kwargs:
            assert kwargs['mol'], 'Mol needed for axes mole fractions'
            assert kwargs['Temp'], 'Temperature needed for IG pressure'
            self.T = kwargs['Temp']
            self.xlabel = ['x (mol/mol)','dx']
            self.mol = kwargs['mol']
            self.units = kwargs['units']
            self.indep = kwargs['indep']
            self.feeds = kwargs['feeds']
            self.path = kwargs['path']
            self.box = kwargs['box']
            self.boxes = kwargs['boxes']

    def getX(self,box):
        rho, N = self.rho[self.feed][self.run], self.N[self.feed][self.run]
        N_tot = sum(N[i][box]['mean'] for i in N.keys())
        x_mean = N[self.mol][box]['mean']/N_tot
        x_stdev = 0.
        for i in N.keys():
            if i == self.mol:
                df_di = (N_tot-N[self.mol][box]['mean'])/pow(N_tot,2)
            else:
                df_di = (-N[self.mol][box]['mean'])/pow(N_tot,2)
            x_stdev += math.pow( df_di,2)*math.pow(N[i][box]['stdev'],2)
        x_stdev = math.sqrt(x_stdev)
        return {'mean':x_mean, 'stdev':x_stdev}


    def Pig_xy(self):
        nIndep = gen_data['numIndep']
        rho, N = self.rho[self.feed][self.run], self.N[self.feed][self.run]
        vapor_box = self.findVapBox( rho, self.mol)
        num_box = len(rho.keys())
        if num_box == 3:
            liquid_box = 'box2'
        else:
            liquid_box = 'box1'
        # get pressure
        p_mean = 0.; p_stdev = 0.
        factor = 1/N_av*R['nm**3*kPa/(mol*K)']*self.T
        for mol in rho.keys():
            mean, stdev = rho[mol][vapor_box]['mean'], rho[mol][vapor_box]['stdev']
            p_mean += mean*factor
            p_stdev += math.pow( stdev, 2)*math.pow(factor, 2)
        p_stdev = math.sqrt(p_stdev)
        file_description = '%s     Pig    %s     dP'%(self.xlabel[0],
                                                       self.xlabel[1])
        file_name = 'Pig_x_mol%s'%self.mol
        X = self.getX(liquid_box)
        writeAGR([X['mean']],[p_mean],
        [calc95conf(X['stdev'], nIndep)], [calc95conf(p_stdev, nIndep)],
        [self.feed], file_name, file_description)
        file_name = 'Pig_y_mol%s'%self.mol
        X = self.getX(vapor_box)
        writeAGR([X['mean']],[p_mean],
        [calc95conf(X['stdev'], nIndep)], [calc95conf(p_stdev, nIndep)],
        [self.feed], file_name, file_description)



from MCFlow.runAnalyzer import checkRun, calc95conf
from MCFlow.chem_constants import N_av, R
from MCFlow.parser import Plot
import numpy as np
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
    elif args['xaxis'] == 'rho':
        my_plotter = RhoBoxAds(**args)
    elif args['xaxis'] == 'Pbox':
        my_plotter = GasBoxAds(**args)
    elif args['xaxis'] == 'Pi':
        my_plotter = GasMolAds(**args)
    elif args['xaxis'] == 'Q':
        my_plotter = LoadAds(**args)
    elif args['xaxis'] == 'X':
        my_plotter = MoleFrac(**args)
    my_plotter.readDBs()

    for feed in args['feeds']:
        # determine if run has been completed
        my_plotter.run = checkRun(args['type'], my_plotter.variables, feed)
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
