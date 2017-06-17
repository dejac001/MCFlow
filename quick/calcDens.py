from MCFlow.runAnalyzer import calc95conf
from MCFlow.chem_constants import R, N_av
import shelve
import sys, os, math

assert (( len(sys.argv) == 2 )), 'usage: python3 %s feed '%sys.argv[0]

feed = sys.argv[1]

with shelve.open('rho-data.db') as db:
    for run, data in db[feed].items():
        if run != 'time':
            str1 = 'run:%s'%run
            for smol, vals in data.items():
                with shelve.open('general-data.db') as gdb:
                    MW = gdb[feed][run]['molecular weight'][smol]
                conv = 1/N_av*math.pow(10,21)*MW
                str2 = '%s mol:%s'%(str1,smol)
                for box, values in vals.items():
                    nIndep = len(values['raw'])
                    error = calc95conf(values['stdev'],nIndep)
                    mean, error = values['mean']*conv, error*conv
                    str3 = '%s %s ==> %e +/- %e'%(
                            str2,box,mean,error
                                )
                    print(str3)
