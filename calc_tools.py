def g_mL(nmlcl, boxlx, MW=0.):
    Concs = [   (n/N_av*MW) / ((math.pow(b,3))*10**(-24)) for (n,b) in zip(nmlcl, boxlx)]
    return Concs

def calcRelDiff(val1, val2):
    '''
    define rel. diff. as abs difference over mean
    '''
    return abs( (val1 - val2) / ( (val1+val2)/2 ) )

import math
from .chem_constants import N_av
