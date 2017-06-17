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

def errorPropDiv(num, den, num_error, den_error):
    f = num/den
    stdev = f*np.sqrt( math.pow(num_error/num,2) + math.pow(den_error/den,2) )
    return stdev

class IdealGasAds:
    def __init__(self, **kwargs):

        # TODO: make general db reader that all files can do
        self.N = {}; self.rho = {}; self.gen_data = {}; self.dG = {}; self.K = {}; self.P = {}
        self.files = ['N-data.db','rho-data.db','general-data.db','dG-data.db','K-data.db','P-data.db']
        self.variables = [self.N, self.rho, self.gen_data, self.dG, self.K, self.P]
        if kwargs:
            # assert kwargs['mol'], 'Mol needed for number density to calculate P assuming I.G.'
            # assert kwargs['Temp'], 'Temperature needed for number density to calculate P assuming I.G.'
            # self.xlabel = ['Pig-mol%s_%4.1f'%(kwargs['mol'],kwargs['Temp']),'dP']
            self.T = kwargs['Temp']
            self.mol = kwargs['mol']
            self.indep = kwargs['indep']
            self.units = kwargs['units']
            self.feeds = kwargs['feeds']
            self.path = kwargs['path']
            self.box = kwargs['box']
            self.boxes = kwargs['boxes']
            self.film = kwargs['film']

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
        molecules = sorted([i for i in N.keys()
                                    if ((N[i][boxTo]['mean'] > 0.)
                          and (N[i][boxFrom]['mean'] > 0.))      ])
        for mol in molecules:
            file_name = 'dG-mol%s_vs_%s_%s-->%s.dat'%(mol, self.xlabel[0], boxFrom, boxTo)
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
        assert self.units, 'Units must be defined for isotherm'
        N = self.N[self.feed][self.run]
        gen_data = self.gen_data[self.feed][self.run]
        file_description = '%s    Q(%s)    %s     dQ'%(self.xlabel[0], self.units,
                                                       self.xlabel[1])
        X = self.getX()
        nIndep = gen_data['numIndep']
        mols_adsorbed = self.getMolAds(N)
        if (self.film == True):
            rho = self.rho[self.feed][self.run]
            zeolite_mass_g = gen_data['zeolite']['mass (g)'] + extra_film_mass
        else:
            zeolite_mass_g = gen_data['zeolite']['mass (g)']
        for mol in mols_adsorbed:
            file_name  = 'Qmol%s-%s-w%sin-zeo-vs-%s.dat'%(mol,self.units,''.join(mols_adsorbed),
                                                           self. xlabel[0][:self.xlabel[0].find('(')])
            if (self.film == True):
                # subtract the no. of molecules in vapor space
                vapor_box = self.findVapBox(rho, mol)
                numInVap = rho[mol][vapor_box]['mean']*film_vapor_volume
            if self.units == 'molec/uc':
                qfactor = 1/gen_data['zeolite']['unit cells']
            elif self.units == 'g/g':
                qfactor = ((gen_data['molecular weight'][mol]/N_av) /
                            zeolite_mass_g )
            elif self.units == 'mol/kg':
                if self.film == True:
                    qfactor = (1/N_av)/(zeolite_mass_g/1000.)
                else:
                    qfactor = gen_data['zeolite'][' mol/kg / 1 mlcl adsorbed']
            if (0 in self.indep) and (len(self.indep) == 1):
                if not self.film:
                    Q_vals = [i*qfactor for i in N[mol]['box1']['raw']]
                    Q_stdev = [N[mol]['box1']['stdev']*qfactor for i in Q_vals]
                else:
                    raise NotImplemented
                writeAGR( X['raw'], Q_vals, None,
                            None, ['%s/%i'%(self.feed,j) for j in
                            range(1, nIndep+1)], file_name, file_description)
            else:
                if not self.film:
                    Q_mean, Q_stdev = (N[mol]['box1']['mean']*qfactor, N[mol]['box1']['stdev']*qfactor)
                else:
                    Q_mean, Q_stdev = ((N[mol]['box1']['mean']-numInVap)*qfactor,
                                        (N[mol]['box1']['stdev']-numInVap)*qfactor)
                if '95conf' in X.keys():
                    dX = X['95conf']
                else:
                    dX = X['stdev']
                writeAGR([X['mean']],[Q_mean],
                         [dX], [calc95conf(Q_stdev, nIndep)],
                         [self.feed], file_name, file_description)

    def DensvX(self):
        assert 'box' in self.box, 'Box needed for density'
        if (0 in self.indep) and (len(self.indep) == 1): raise NotImplemented
        X = self.getX()
        rho = self.rho[self.feed][self.run]
        gen_data = self.gen_data[self.feed][self.run]
        nIndep = gen_data['numIndep']
        file_description = '%s    density(g/mL)    %s     dd'%(self.xlabel[0], self.xlabel[1])
        file_name = 'Dens_v_%s_%s.dat'%(self.xlabel[0], self.box)
        density = {'mean':0.,'stdev':0.}
        for mol in rho.keys():
            factor = gen_data['molecular weight'][mol]/N_av*math.pow(10,21)
            if self.box in rho[mol]:
                dens_molec_nm3 = {'mean':rho[mol][self.box]['mean'],
                                    'stdev':rho[mol][self.box]['stdev']}
                density['mean'] += dens_molec_nm3['mean']*factor
                density['stdev'] += math.pow(dens_molec_nm3['stdev']*factor,2)
        density['stdev'] = math.sqrt( density['stdev'] )
        writeAGR([X['mean']],[density['mean']],
                     [calc95conf(X['stdev'], nIndep)], [calc95conf(density['stdev'], nIndep)],
                     [self.feed], file_name, file_description)

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
                         [self.feed], file_name, file_description)

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
                         [self.feed], file_name, file_description)

    def SvX(self):
        '''
        this is a special case of SvC, where we plot P as opposed to C
        and the boxFrom might not be the same for each molecule
        (ex: in the case of IG or kH adsorption)
        :param P: Pressure of x axis (note: not pressure data)
        '''
        if (0 in self.indep) and (len(self.indep) == 1): raise NotImplemented
        rho = self.rho[self.feed][self.run]
        K = self.K[self.feed][self.run]
        if 'P-box' in self.xlabel[0]:
            P = self.P[self.feed][self.run]
        X = self.getX()
        file_description = '%s    S(%s)    %s     dS'%(self.xlabel[0], self.units, self.xlabel[1])
        gen_data = self.gen_data[self.feed][self.run]
        nIndep = gen_data['numIndep']
        for mol_pair in K['box1']:
            K_to = K['box1'][mol_pair]
            mol1, mol2 = map(int,mol_pair.split('/'))
            file_name = 'S_%s-vs-%s.dat' % (mol_pair, self.xlabel[0])
            if (('P' in self.xlabel[0]) or ('C' in self.xlabel[0])):
                if 'P-box' in self.xlabel[0]:
                    vapor_box_a = self.findVapBox(rho, '%i'%mol1)
                    vapor_box_b = self.findVapBox(rho, '%i'%mol2)
                    pa, pa_stdev = P[vapor_box_a]['mean'], P[vapor_box_a]['stdev']
                    pb, pb_stdev = P[vapor_box_b]['mean'], P[vapor_box_b]['stdev']
                    K_from = {
                        'mean':pa/pb,
                              'stdev': errorPropDiv(pa, pb, pa_stdev, pb_stdev)
                    }
                elif 'P' in self.xlabel[0]:
                    vapor_box = self.findVapBox(rho, '%i'%mol1)
                    K_from = K[vapor_box][mol_pair]
                elif 'C' in self.xlabel[0]:
                    K_from = K['box2'][mol_pair]
                S_mean = K_to['mean'] / K_from['mean']
                S_stdev = errorPropDiv(K_to['mean'], K_from['mean'],
                                       K_to['stdev'], K_from['stdev'])
                writeAGR([X['mean']], [S_mean],
                         [calc95conf(X['stdev'],nIndep)], [calc95conf(S_stdev,nIndep)],
                         [self.feed], file_name, file_description)
            elif 'kH' in self.xlabel[0]:
                if mol_pair == '%s/1'%self.mol:
                    K_to['95conf'] = calc95conf(K_to['stdev'],nIndep)
                    MW_W = self.solventMW
                    MW_D = gen_data['molecular weight'][self.mol]
                    S_mean = (
        K_to['mean']*(self.density[0]/X['mean'] - 1)*MW_D/MW_W
                    )
                    w = X['mean']/self.density[0]
                    S_long = K_to['mean']/(
                    (w/MW_D)/((1-w)/MW_W )
                    )
                    dS_dKto = (self.density[0]/X['mean'] - 1)*MW_D/MW_W
                    dS_drho = K_to['mean']/X['mean']*MW_D/MW_W
                    dS_dC = -K_to['mean']*self.density[0]/math.pow(X['mean'],2)*MW_D/MW_W
                    S_95conf = math.sqrt(
                        math.pow(dS_dKto*K_to['95conf'],2) +
                        math.pow(dS_drho*self.density[1],2) +
                        math.pow(dS_dC*X['95conf'],2)
                    )
                    writeAGR([X['mean']], [S_mean],
                             [X['95conf']], [S_95conf],
                             [self.feed], file_name, file_description)
            else:
                print(self.xlabel, 'not considered')

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
            self.film = kwargs['film']
    def getX(self):
        return self.rho[self.feed][self.run][self.mol]['box%s'%self.box]

class GasBoxAds(IdealGasAds):
    def __init__(self, **kwargs):
        IdealGasAds.__init__(self, **kwargs)
        self.vapor_box = self.box
        self.xlabel= ['P-box(kPa)', 'dP-box(kPa)']
        if kwargs:
            assert self.vapor_box, 'Box needed to calculate pressure (box)'
            self.film = kwargs['film']

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
            self.film = kwargs['film']
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
        if self.mol:
            Q_mean, Q_stdev = (N[self.mol]['box1']['mean']*qfactor, N[self.mol]['box1']['stdev']*qfactor)
        else:
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
        self.N = {}; self.P = {}; self.gen_data = {}; self.rho = {}#; self.dG = {}
        self.files = ['N-data.db','P-data.db','general-data.db','rho-data.db']#, 'dG-data.db']
        self.variables = [self.N, self.P, self.gen_data, self.rho]#, self.dG]
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
            self.film = kwargs['film']
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
        self.dG = {}; self.C = {}; self.N = {}; self.rho = {}; self.gen_data = {}; self.K = {}
        self.files = ['dG-data.db','Conc-data.db','N-data.db','rho-data.db','general-data.db','K-data.db']
        self.variables = [self.dG, self.C, self.N, self.rho, self.gen_data,self.K]
        if kwargs:
            self.xlabel = ['C(g/mL)','dC']
            self.mol = kwargs['mol']
            self.units = kwargs['units']
            self.indep = kwargs['indep']
            self.feeds = kwargs['feeds']
            self.path = kwargs['path']
            self.box = kwargs['box']
            self.boxes = kwargs['boxes']
            self.film = kwargs['film']



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
        #self.dG = {};
        self.N = {}; self.rho = {}; self.gen_data = {}
        self.files = [#'dG-data.db',
    'N-data.db','rho-data.db','general-data.db']
        self.variables = [#self.dG,
                self.N, self.rho, self.gen_data]
        if kwargs:
            assert kwargs['mol'], 'Mol needed for axes mole fractions'
            self.T = kwargs['Temp']
            self.xlabel = ['x(mol/mol)','dx']
            self.mol = kwargs['mol']
            self.units = kwargs['units']
            self.indep = kwargs['indep']
            self.feeds = kwargs['feeds']
            self.path = kwargs['path']
            self.box = kwargs['box']
            self.boxes = kwargs['boxes']
            self.film = kwargs['film']

    def XvX(self):
        assert 'box' in self.box, 'Box needed for mole frac plot'
        nIndep = self.gen_data[self.feed][self.run]['numIndep']
        N = self.N[self.feed][self.run]
        x = self.getX(box=self.box)
        file_description = 'x (mol/mol)    x (mol/mol)    dx     dx'
        for mol in N.keys():
            y = self.getX(box=self.box, myMol=mol)
            if y['mean'] > 0.:
                file_name = 'XvX_%s_%s.dat'%(self.mol,mol)
                writeAGR([x['mean']],[y['mean']],
                        [calc95conf(x['stdev'], nIndep)], [calc95conf(y['stdev'], nIndep)],
                        [self.feed], file_name, file_description)

    def getX(self,box=None,myMol=None):
        if myMol==None: myMol = self.mol
        if box == None: box = self.box
        N = self.N[self.feed][self.run]
        N_tot = sum(N[i][box]['mean'] for i in N.keys())
        x_mean = N[myMol][box]['mean']/N_tot
        x_stdev = 0.
        for i in N.keys():
            if i == myMol:
                df_di = (N_tot-N[myMol][box]['mean'])/pow(N_tot,2)
            else:
                df_di = (-N[myMol][box]['mean'])/pow(N_tot,2)
            x_stdev += math.pow( df_di,2)*math.pow(N[i][box]['stdev'],2)
        x_stdev = math.sqrt(x_stdev)
        return {'mean':x_mean, 'stdev':x_stdev}

    def Txy(self):
        nIndep = self.gen_data[self.feed][self.run]['numIndep']
        rho, N = self.rho[self.feed][self.run], self.N[self.feed][self.run]
        try:
            vapor_box = self.findVapBox( rho, self.mol)
        except KeyError:
            vapor_box = 'box2'
        num_box = len(rho['1'].keys())
        if num_box == 3:
            liquid_box = 'box2'
        else:
            liquid_box = 'box1'
        T = self.gen_data[self.feed][self.run]['temperature']
        file_description = '%s     T    %s '%(self.xlabel[0],
                                                       self.xlabel[1])
        file_name = 'Tx_mol%s'%self.mol
        X = self.getX(box=liquid_box)
        writeAGR([X['mean']],[T],
        [calc95conf(X['stdev'], nIndep)], [0.],
        [self.feed], file_name, file_description)
        file_name = 'Ty_mol%s'%self.mol
        X = self.getX(box=vapor_box)
        writeAGR([X['mean']],[T],
        [calc95conf(X['stdev'], nIndep)], [0.],
        [self.feed], file_name, file_description)


    def Pig_xy(self):
        nIndep = self.gen_data[self.feed][self.run]['numIndep']
        self.T = self.gen_data[self.feed][self.run]['temperature']
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
        X = self.getX(box=liquid_box)
        writeAGR([X['mean']],[p_mean],
        [calc95conf(X['stdev'], nIndep)], [calc95conf(p_stdev, nIndep)],
        [self.feed], file_name, file_description)
        file_name = 'Pig_y_mol%s'%self.mol
        X = self.getX(box=vapor_box)
        writeAGR([X['mean']],[p_mean],
        [calc95conf(X['stdev'], nIndep)], [calc95conf(p_stdev, nIndep)],
        [self.feed], file_name, file_description)


from MCFlow.runAnalyzer import checkRun, calc95conf
from MCFlow.chem_constants import N_av, R
from MCFlow.parser import Plot
import numpy as np
import os, math, shelve
extra_film_mass = 24/N_av*1.0079 + 36*2/N_av*1.0079
film_vapor_volume = 4.*4.*4.

if __name__ == '__main__':
    # TODO: find way to have conditional arguments based off of what xaxis and yaxis are
    # (in other words, we have an excess of variables as arguments)
    my_parser = Plot()
    my_parser.axes()
    my_parser.isotherm()
    my_parser.parser.add_argument('-NS','--film',help='whether or not the zeolite box is a nanosheet',type=bool,default=False)

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
    elif args['xaxis'] == 'x':
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
        elif args['yaxis'] == 'Txy':
            my_plotter.Txy()
        elif args['yaxis'] == 'dens':
            my_plotter.DensvX()
