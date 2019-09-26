import math


def calc95conf(stdev, numIndep):
    T_values = {'8': 2.365, '4': 3.182, '16': 2.131, '32': 2.04, '2': 4.303,
                '6': 2.571, '7': 2.447, '1': 1.0, '5': 2.776}
    assert '%i' % numIndep in T_values.keys(), 'No T-value stored for %i indep' % numIndep
    return stdev / math.pow(numIndep, 0.5) * T_values['%i' % numIndep]
