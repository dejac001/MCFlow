if __name__ == '__main__':
    import argparse
    from MCFlow.chem_constants import R, N_av
    parser = argparse.ArgumentParser(description='Calculate rho (molec/nm**3) for IG given p, T')
    parser.add_argument('-T','--temperature',help='temperature in Kelvin',type=float)
    parser.add_argument('-p','--pressure',help='pressure in bar',type=float)
    args = vars(parser.parse_args())

    # PV = nRT
    # P/RT = n / V
    P = args['pressure']*100 # convert from bar to kPa
    rho = P    / ( R['nm**3*kPa/(mol*K)'] * args['temperature'] ) # mol / nm**3
    rho = rho * N_av
    print('for p = %e bar, T = %e K, rho = %e molec/ nm**3'% (args['pressure'],
            args['temperature'], rho))
