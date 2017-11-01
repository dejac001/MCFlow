# '669': {'mean': 0.028571428571428574,
#                              'raw': [0.0, 0.04, 0.0, 0.02, 0.04, 0.1, 0.0],
#                              'stdev': 0.033563828927059232},
#                      'sum': 764.7600000000009},
#             'total frames': 299},

def rho_v_r(data, uc_vectors, center, name):
    symbols = ['o','x','s','<','>','^','D','v']
    uc_center = [i/2 for i in uc_vectors]
    key = [i for i in data.keys() if 'prod-' in i][0]
    for mol, values in data[key].items():
        begin = True
        if 'mol' not in mol: continue
        LMN = list(map(int,str(max(int(k) for k in values.keys() if k != 'sum'))))
        supercell = [i + i*j for (i,j) in zip(uc_vectors,LMN)]
        figure = plt.figure()
        ax = figure.add_subplot(111)
        rho_means = []
        rho_stdev = []
        r = []
        for lmn, vals in values.items():
            if lmn == 'sum': continue
            l, m, n = map(int,list(lmn))
            my_center = [uc_center[0]+l*uc_vectors[0],
                         uc_center[1]+m*uc_vectors[1],
                         uc_center[2]+n*uc_vectors[2]]
            center_center_distance = calculate_distance(center, my_center, supercell)
            rho_means.append(vals['mean'])
            rho_stdev.append(vals['stdev'])
            r.append(center_center_distance)
            if begin == True:
                begin = False
                rho = [[] for i in vals['raw']]
                r_raw = [[] for i in vals['raw']]
            for i, val in enumerate(vals['raw']):
                rho[i].append(val)
                r_raw[i].append(center_center_distance)
#       for x, y, sym in zip(r_raw, rho, symbols):
#           ax.plot(x,y, sym)
        ax.errorbar(r, rho_means, yerr=rho_stdev, fmt='ro')
        plt.subplots_adjust(right=0.97)
        set_x_ticks(ax, [0,20, 40, 60, 80, 100, 120])
        ax.set_xlabel('$r_{\mathrm{uc-uc}}\mathrm{\;[\;\AA\;]}$')
        ax.set_ylabel('$Q\mathrm{\;[\;molec\;/\;uc\;]}$')
        figure.savefig(name)


from MCFlow.calc_tools import calculate_distance
from plotting.util import set_x_ticks
import matplotlib.pyplot as plt
import numpy as np


if __name__ == '__main__':
    import argparse, shelve
    parser = argparse.ArgumentParser()
    parser.add_argument('-fi','--file',help='db file to open',type=str)
    parser.add_argument('-feed','--feed',help='feed',type=str)
    parser.add_argument('-n','--name',type=str, default='fig.pdf')
    parser.add_argument('-c','--center',help='center to calculate distance from',type=float,nargs='+')
    parser.add_argument('-abc','--vectors', help='unit cell vectors. For folding coordinates into a unit cell',
                                      type = float, nargs = '+', default = [20.022,19.899,13.383])
    args = vars(parser.parse_args())

    with shelve.open(args['file']) as db:
        data = db[args['feed']]

    rho_v_r(data, args['vectors'], args['center'], args['name'])
