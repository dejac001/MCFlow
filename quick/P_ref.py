if __name__ == '__main__':
    import argparse
    from MCFlow.chem_constants import R, N_av
    parser = argparse.ArgumentParser(description='Calculate rho (molec/nm**3) for IG given p, T')
    parser.add_argument('-T','--temperature',help='temperature in Kelvin',type=float)
    parser.add_argument('-r','--rho',help='number density in molec / nm**3',type=float)
    args = vars(parser.parse_args())

    # PV = nRT
    # P = rho*RT
#   P = args['pressure']*100 # convert from bar to kPa
    P = args['rho']/N_av*R['nm**3*kPa/(mol*K)']*args['temperature']
    print('for rho  = %e molec/nm**3, T = %e K, P = %e kPa'% (args['rho'],
            args['temperature'], P))
