N_av = 6.022*10**23
R = {}
R['kJ/(mol*K)'] = 8.314/1000 # kJ/(mol*K)
R['cm**3*kPa/(mol*K)'] = 8314.
R['nm**3*kPa/(mol*K)'] = R['cm**3*kPa/(mol*K)']*10**21
R['\AA**3*kPa/(mol*K)'] = R['nm**3*kPa/(mol*K)']*1000
R['\AA**3*MPa/(mol*K)'] = R['\AA**3*kPa/(mol*K)']/1000
