# '669': {'mean': 0.028571428571428574,
#                              'raw': [0.0, 0.04, 0.0, 0.02, 0.04, 0.1, 0.0],
#                              'stdev': 0.033563828927059232},
#                      'sum': 764.7600000000009},
#             'total frames': 299},

def rho_v_r(data, feed):
    key = [i for i in data.keys() if 'prod-' in i][0]
    figure = plt.figure()
    values = data[key]['mol8']
    LMN = list(map(int,str(max(int(k) for k in values.keys() if k != 'sum'))))
    supercell = [i + i*j for (i,j) in zip(uc_vectors,LMN)]
    H = np.zeros((7,7,10))
    S = np.zeros((7,7,10))
    for lmn, vals in values.items():
        if lmn == 'sum': continue
        l, m, n = map(int,list(lmn))
        my_center = [l+2,
                     m+2,
                     n+2]
        if my_center[0] > 6: my_center[0] =  my_center[0] - 7
        if my_center[1] > 6: my_center[1] = my_center[1] - 7
        if my_center[2] > 9: my_center[2] = my_center[2] - 10
        i,j,k = my_center
        H[i,j,k] = vals['mean']
        S[i,j,k] = vals['stdev']
    print('max S:',np.max(S))
    print('min S:',np.min(S))
    print('max H:',np.max(H))
    print('min H:',np.min(H))
    nx, ny, nz = np.shape(H)
    Z, Y = np.meshgrid(list(range(nz+1)), list(range(ny+1)))
     # xy histogram
    nrows, ncolumns= nx, 2
    hmin, hmax = 0., np.max(H)
    smin, smax = 0., np.max(S)
    xlabel, ylabel = '$\mathbf{c}$', '$\mathbf{b}$'
    for x in range(nx):
        ax_mean = figure.add_subplot(nrows,ncolumns,x*ncolumns+1)
        ax_stdev = figure.add_subplot(nrows,ncolumns,x*ncolumns+2)
        h = H[x,:,:]
        s = S[x,:,:]
        image_H  = ax_mean.pcolormesh(Z,Y,h,vmin=hmin,vmax=hmax,cmap='plasma')
        image_S  = ax_stdev.pcolormesh(Z,Y,s,vmin=smin,vmax=smax,cmap='plasma')
        for ax in [ax_mean, ax_stdev]:
            ax.tick_params(axis='both',which='both',direction='out')
            set_x_ticks(ax, [0,2,4,6,8,10])
            set_y_ticks(ax, [0,1,2,3,4,5,6,7])
            ax.tick_params(axis='y',which='minor',left='off',right='off')
        label_subplot_axes(ax_mean, x*ncolumns+1, ncolumns, nrows, xlabel, ylabel, tight=True)
        label_subplot_axes(ax_stdev, x*ncolumns+2, ncolumns, nrows, xlabel, ylabel, tight=True)
        ax_stdev.tick_params(axis='both',which='both',left='off',right='off')
        minimum_ticks(ax_stdev, x*ncolumns+2, ncolumns, nrows)
        minimum_ticks(ax_mean, x*ncolumns+1, ncolumns, nrows)
    my_bottom = 0.055
    my_top=0.99
    plt.subplots_adjust(left=0.11,right=0.7,hspace=0.,wspace=0.,
                        bottom=my_bottom, top=my_top)
    figure.set_size_inches(3.25,7.0)


    if '450' in feed:
        ticks=[0,0.4,0.8,1.2,1.6,2.0]
        s_ticks=np.arange(0.0,2.0,0.1)
        cbar_ax = figure.add_axes([0.82, my_bottom, 0.06, my_top-my_bottom])
    else:
        cbar_ax = figure.add_axes([0.8, my_bottom, 0.08, my_top-my_bottom])
        ticks=[0,1,2,3,4,5,6,7,8]
        s_ticks=np.arange(0.0,2.0,0.2)
    cbar = figure.colorbar(image_H, cax=cbar_ax, orientation='vertical')
    cax2 = cbar_ax.twinx()
    cbar.set_ticks(ticks)
    cbar.set_label('$Q\mathrm{\;[\;molec\;/\;uc\;]}$')
    cbar.ax.yaxis.set_label_position("left")
#   cbar_ax.tick_params(which='both',direction='out')
#   cbar_ax.tick_params(which='minor',left='off')
#   cax2.set_yticks(s_ticks)
    set_y_ticks(cax2, s_ticks)
#   cax2.set_yticks([0.1+i for i in s_ticks], minor=True)
    cax2.set_ylim([smin, smax])
    cax2.set_ylabel('$\sigma_Q\mathrm{\;[\;molec\;/\;uc\;]}$')
    figure.savefig('FIG-test.pdf')


from MCFlow.calc_tools import get_vector
from plotting.util import set_x_ticks, set_y_ticks, label_subplot_axes, minimum_ticks
import matplotlib.pyplot as plt
import numpy as np

uc_vectors = [20.022,19.899,13.383]

if __name__ == '__main__':
    import argparse, shelve
    parser = argparse.ArgumentParser()
    parser.add_argument('-feed','--feed',help='feed',type=str,default='7x7x10')
    parser.add_argument('-f','--file',help='feed',type=str,default='uc-dens-data.db')
    args = vars(parser.parse_args())

    with shelve.open(args['file']) as db:
        data = db[args['feed']]

    rho_v_r(data, args['feed'])
