if __name__ == '__main__':
    import argparse
    from MCFlow.chem_constants import R, N_av
    parser = argparse.ArgumentParser(description='Calculate rho (molec/nm**3) for IG given p, T')
    parser.add_argument('-T','--temperature',help='temperature in Kelvin',type=float)
    parser.add_argument('-r','--rho',help='density in molec/ angst.**3',type=float)
    args = vars(parser.parse_args())

    # PV = nRT
    # P/RT = n / V
    rho = args['rho']/N_av # mol / nm**3
    P = rho *R['nm**3*kPa/(mol*K)'] *args['temperature']
    print('for rho = %e molec/anst.**3, T = %e K, p = %e kPa'%(args['rho'], args['temperature'], P))
