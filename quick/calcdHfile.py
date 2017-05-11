import numpy as np
from MCFlow.runAnalyzer import calc95conf
import sys
if __name__ == '__main__':
    import os
    T = sys.argv[-1]
    if 'PycharmProjects' in os.getcwd():
        os.chdir('/Users/dejacor/Documents/BETO-sim'
                 'ulations/desorption/vapor-desorption/363K/')
    numIndep = 8
    with open('dH-on-the-fly.dat','w') as f:
        for L in [i for i in os.listdir() 
                        if os.path.isfile('%s/1/dH_raw.dat'%i)]:
            dH = {'raw data':[]}
            for seed in range(1,numIndep+1):
                dH12 = []
                for line in open('%s/%i/dH_raw.dat'%(L,seed)):
                    dH12.append( float(line.split()[-3]))
                f.write('%e #%s/%s/%i\n'%(np.mean(dH12), T,L,seed))
#               dH['raw data'].append( np.mean(dH12) )
#           mean = np.mean( dH['raw data'] )
#           stdev = np.std( dH['raw data'] )
#           f.write('%e %e #%s\n'%(mean, calc95conf(stdev, numIndep), L))
