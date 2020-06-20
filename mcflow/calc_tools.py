def g_mL(nmlcl, boxlx, MW=0.):
    Concs = [(n / N_av * MW) / ((math.pow(b, 3)) * 10 ** (-24)) for (n, b) in zip(nmlcl, boxlx)]
    return Concs


def calcRelDiff(val1, val2):
    '''
    define rel. diff. as abs difference over mean
    '''
    return abs((val1 - val2) / ((val1 + val2) / 2))


def pbc(coord, boxlength):
    if coord >= boxlength:
        folded = coord - boxlength
    elif coord < 0:
        folded = coord + boxlength
    else:
        folded = coord
    return folded


def fold(xyz, boxlengths):
    new_coords = []
    for i, boxl in zip(xyz, boxlengths):
        new_coords.append(pbc(i, boxl))
    return new_coords


def fold_coords(xyz, boxlengths):
    new_coords = []
    for coords in xyz:
        new_coords.append(fold(coords, boxlengths))
    return np.matrix(new_coords)


def get_vector(xyz1, xyz2, abc):
    vector = [xyz1[i] - xyz2[i] for i in range(len(xyz1))]
    for i in range(len(abc)):
        # implement minimum image
        if vector[i] > abc[i] / 2.:
            # positive
            vector[i] = abc[i] - vector[i]
        elif vector[i] < -1. * abc[i] / 2.:
            vector[i] = -1 * (abc[i] + vector[i])
    return vector


def calculate_distance(xyz1, xyz2, abc):
    my_vector = get_vector(xyz1, xyz2, abc)
    return np.linalg.norm(my_vector, 2, 0)


def calculate_distance2(xyz1, xyz2, abc):
    my_vector = get_vector(xyz1, xyz2, abc)
    return sum(r * r for r in my_vector)


def convertg_mL_to_x(c, MW):
    '''
    c: concentration in g/mL
    MW: molecular weight of solute in g/mol
    assume solvent is water (MW=18.02g/mol), density is 1 g/mL
    '''
    rho_solution = 1.0  # g/mL -- assumption
    v = 100  # basis
    mol_solute = c * v / MW
    mol_water = rho_solution * v / 18.02
    return mol_solute / (mol_solute + mol_water)


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values - average) ** 2, weights=weights)  # Fast and numerically precise
    return (average, math.sqrt(variance))


def calculate_angle(c1, c2, c3, abc):
    '''

    :param c1: xyz coords of bead 1
    :param c2: xyz coords of middle bead
    :param c3: xyz coords of final bead
    :param abc:
    :return: angle in radians!!; limits are from 0 to pi
    '''
    a = get_vector(c2, c1, abc)
    b = get_vector(c2, c3, abc)
    theta = np.arccos(
        np.dot(a, b) / (np.linalg.norm(a, 2, 0) * np.linalg.norm(b, 2, 0))
    )
    return theta


def determine_gauche(phi):
    if (phi <= 120. * np.pi / 180.) or (phi >= 240. * np.pi / 180):
        return True  # gauche
    else:
        return False  # Trans


def calculate_torsion_poor(vectors):
    # phi convention trappe website; 180 is trans
    b1, b2, b3 = vectors
    b1_c_b2 = np.cross(b1, b2)
    b2_c_b3 = np.cross(b2, b3)
    denominator = np.linalg.norm(b1_c_b2, ord=2) * np.linalg.norm(b2_c_b3, ord=2)
    numerator = np.dot(b1_c_b2, b2_c_b3)
    return np.arccos(numerator / denominator)


def calculate_torsion(vectors):
    # phi convention trappe website; 180 is trans
    b1, b2, b3 = vectors
    b1_c_b2 = np.cross(b1, b2)
    b2_c_b3 = np.cross(b2, b3)
    atan = np.arctan2(np.dot(np.cross(b1_c_b2, b2_c_b3), b2 / np.linalg.norm(b2, ord=2)),
                      np.dot(b1_c_b2, b2_c_b3))
    # if abs(atan) < 1e-8:
    #     print(atan)
    #     print(b1,b2,b3)
    #     print(np.cross(b1,b2),np.linalg.norm(b1-b2,1))
    #     print(np.linalg.norm(b1,ord=1)+np.linalg.norm(b2,ord=1))
    #     quit()
    if atan < 0.:
        return atan + 2 * np.pi
    else:
        return atan


def ePropC2K(mean, error):
    '''
    do error propagation when converting celcius to kelvin
    '''
    df_dT = -1 / math.pow(mean + 273.15, 2)
    error_new = abs(df_dT * error)
    return error_new


def eProp_division(A, sA, B, sB):
    '''
    error propagation when dividing A by B
    '''
    #   abs_f = abs(A/B)
    #   factor = math.sqrt(
    #       math.pow(sA/A,2) + math.pow(sB/B,2)
    #   )
    #   return abs_f*factor
    B2 = B * B
    return math.sqrt(
        1 / B2 * sA * sA +
        math.pow(-A / B2, 2) * sB * sB
    )


import numpy as np
import math
from chem_constants import N_av
import random

abc = [40., 39., 40.]
x1, y1, z1 = [random.random() * i for i in abc]
x2, y2, z2 = [random.random() * i for i in abc]
vector = [x2 - x1, y2 - y1, z2 - z1]

if __name__ == '__main__':
    import timeit

    print(timeit.timeit('ca(vector)', setup="from __main__ import ca, vector"))
    print(timeit.timeit('cb(vector)', setup="from __main__ import cb, vector"))
    print(timeit.timeit('np.linalg.norm(vector, 2, 0)', setup='''from __main__ import vector, np'''))
