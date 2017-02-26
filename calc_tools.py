def g_mL(nmlcl, boxlx, MW=0.):
    Concs = [   (n/N_av*MW) / ((math.pow(b,3))*10**(-24)) for (n,b) in zip(nmlcl, boxlx)]
    return Concs

def calcRelDiff(val1, val2):
    '''
    define rel. diff. as abs difference over mean
    '''
    return abs( (val1 - val2) / ( (val1+val2)/2 ) )

def fold(xyz, boxlengths):
    def pbc(coord, boxlength):
        if coord > boxlength:
            folded = coord - boxlength
        elif coord < 0:
            folded = coord + boxlength
        else:
            folded = coord
        return folded
    new_coords = []
    for i, boxl in zip(xyz, boxlengths):
        new_coords.append( pbc(i,boxl) )
    return new_coords

def get_vector(xyz1, xyz2, abc):
    vector = [xyz1[i] - xyz2[i] for i in range(len(xyz1))]
    for i in range(len(abc)):
        # implement minimum image
        if vector[i] > abc[i]/2.:
            # positive
            vector[i] = abc[i] - vector[i]
        elif vector[i] < -1.*abc[i]/2.:
            vector[i] = abc[i] + vector[i]
    return vector

def calculate_distance(xyz1, xyz2, abc):
    my_vector = get_vector(xyz1, xyz2, abc)
    return np.linalg.norm(my_vector, 2, 0)

def calculate_distance2(xyz1, xyz2, abc):
    my_vector = get_vector(xyz1, xyz2, abc)
    return sum(r*r for r in my_vector)

def convertg_mL_to_x(c, MW):
    '''
    c: concentration in g/mL
    MW: molecular weight of solute in g/mol
    assume solvent is water (MW=18.02g/mol), density is 1 g/mL
    '''
    rho_solution = 1.0 # g/mL -- assumption
    v = 100 # basis
    mol_solute = c*v/MW
    mol_water = rho_solution*v/18.02
    return mol_solute/(mol_solute + mol_water)

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, math.sqrt(variance))

def calculate_angle(c1,c2,c3,abc):
    '''

    :param c1: xyz coords of bead 1
    :param c2: xyz coords of middle bead
    :param c3: xyz coords of final bead
    :param abc:
    :return: angle in radians!!
    '''
    a = get_vector(c2,c1,abc)
    b = get_vector(c2,c3,abc)
    theta = np.arccos(
        np.dot(a,b) / ( np.linalg.norm(a,2,0)*np.linalg.norm(b,2,0) )
    )
    return theta


import numpy as np
import math
from MCFlow.chem_constants import N_av
import random
abc = [40., 39., 40.]
x1, y1, z1 = [random.random()*i for i in abc]
x2, y2, z2 = [random.random()*i for i in abc]
vector = [x2-x1, y2-y1, z2-z1]

if __name__ == '__main__':
    import timeit
    print(timeit.timeit('ca(vector)',setup="from __main__ import ca, vector"))
    print(timeit.timeit('cb(vector)',setup="from __main__ import cb, vector"))
    print(timeit.timeit('np.linalg.norm(vector, 2, 0)',setup='''from __main__ import vector, np'''))
