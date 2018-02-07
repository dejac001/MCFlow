
N_av = 6.022*10**23
R = {}
R['kJ/(mol*K)'] = 8.314/1000 # kJ/(mol*K)
R['cm**3*kPa/(mol*K)'] = 8314.
R['nm**3*kPa/(mol*K)'] = R['cm**3*kPa/(mol*K)']*10**21
R['\AA**3*kPa/(mol*K)'] = R['nm**3*kPa/(mol*K)']*1000
R['\AA**3*MPa/(mol*K)'] = R['\AA**3*kPa/(mol*K)']/1000
g_uc_MFI = (2304/N_av*15.9994 + 1152/N_av*28.055)/(2*2*3)


from MCFlow.chem_constants import R, N_av

import sys
assert len(sys.argv) == 4, 'usage: python %s N p(kPa) T(K)'%sys.argv[0]

N, p, T = map(float, sys.argv[-3:])

moles = N/N_av

V = R['\AA**3*kPa/(mol*K)']*T*moles/p

import math
print('boxlength should be {}'.format(math.pow(V, 1/3)))
