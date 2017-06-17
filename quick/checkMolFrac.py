import sys, os

assert (len(sys.argv) == 2), 'usage: python3 %s feed'%sys.argv[0]

feed_test = sys.argv[1]
import shelve, pprint

with shelve.open('X-data.db') as db:
    for run in db[feed_test].keys():
        if 'time' not in run:
            for box, vals in db[feed_test][run].items():
                x = 0.
                for mol in vals.keys():
                    x += vals[mol]['mean']
                assert abs(x-1.) < 1e-5, 'Mole fractions for feed %s and %s are not consistent'%(feed_test,box)
