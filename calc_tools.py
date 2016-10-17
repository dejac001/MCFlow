def g_mL(nmlcl, boxlx, MW=0.):
    Concs = [   (n/N_av*MW) / ((math.pow(b,3))*10**(-24)) for (n,b) in zip(nmlcl, boxlx)]
    return Concs

import math
from chem_constants import N_av