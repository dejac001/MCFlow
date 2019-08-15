from MCFlow.chem_constants import N_av

def themGhosts(vapor_volume,T,P):
    '''
    assume vapor volume will take size of ghost molecules if ideal gas
    :param vapor_volume: in \AA**3
    :value 83.14 is cm**3*bar/mol/K                                                                                      
    :param P in bar
    :return: number of ghost molecules
    '''
    nGhost = int(   vapor_volume*10**(-24)*P/(83.14*T)*N_av    )
    return nGhost

if __name__ == '__main__':

    V = 32/(3.334962e-06 )*1000
    import math
    boxlx = math.pow(V, 1/3)
    print('boxlength', boxlx)
    print(themGhosts(V, 323., 3.0))
