def plotColorHist(data):
    fig, ax = plt.subplots()

    my_min = 10**8; my_max = -10**8
    for i in data:
        imin,imax = min(i),max(i)
        if imin < my_min:
            my_min = imin
        if imax > my_max:
            my_max = imax
    print('min and max are {} and {}'.format(my_min,my_max))
    ibins = np.linspace(-40,5)
    colors = ['cyan','magenta']
    legend = ['ButOH','Water']
    lines = ['solid','dashed']
    for i, values in enumerate(data):
        print('number of colored values is %i'%len(values))
        n, bins  = np.histogram(values, bins=ibins, density=False)#, log=True)
        n, new_bins = hist_norm_height(n, bins, max(n))
        ax.step(new_bins, n, color=colors[i],label=legend[i],linestyle=lines[i])
    ax.set_xlabel('$\Delta G\mathrm{\;[\;kJ\;/\;mol\;]}$')
    ax.set_ylabel('$h(\Delta G)$')
    ax.legend(loc=0,fontsize='x-small')
    ax.tick_params(colors='black',width=0.5)
    ax.tick_params(direction='in',which='both',left=True,right=True)
    ax.tick_params(axis='x',direction='out',bottom=True, top=True,which='both')
    for dirn in ['top','bottom','left','right']:
        ax.spines[dirn].set_color('black')
        ax.spines[dirn].set_linewidth(0.5)
    y_ticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    y_ticks_minor = [0.1,0.3,0.5,0.7,0.9]
    ax.set_yticks(y_ticks)
    ax.set_yticks(y_ticks_minor,minor=True)
    x_ticks = [-40,-30,-20,-10,0]
    x_ticks_minor = [-35,-25,-15,-5,5]
    ax.set_xticks(x_ticks)
    ax.set_xticks(x_ticks_minor,minor=True)
    ax.set_xlim([ibins[0],ibins[-1]])
    ax.set_ylim([0., 1.])
    # ax.set_ylim([0.9,200.])
    plt.subplots_adjust(left=0.08,right=0.995, top=0.95,bottom=0.21)
    fig = plt.gcf()
    fig.set_size_inches(4.7,2.)
    fig.savefig('dG_hist.png',dpi=600)
    # plt.show()

import numpy as np
import matplotlib.pyplot as plt
from MCFlow.selectivity_histogram import hist_norm_height #,SHist
from MCFlow.file_formatting.writer import vtkRectilinearMesh
from MCFlow.parser import Structure
from MCFlow.file_formatting import reader
from MCFlow.probhistogram import getCoords, getFileName, dG3D, Hist3D

if __name__ == '__main__':
    my_parser = Structure()
    my_parser.parser.add_argument('-n', '--numFrames', help='total number of frames', type=int)
    my_parser.parser.add_argument('-T', '--Temp', help='Temperature [ K ]', type=float)
    my_parser.parser.add_argument('-f1','--files',help='files',type=str,nargs='+')
    my_parser.parser.add_argument('-off','--offAxis',help='axis to average along',
                                  default=[],nargs='+')
    my_parser.parser.add_argument('-ref','--reference',help='reference density ',type=float,nargs='+')
    args = vars(my_parser.parse_args())
    if args['vectors'] == None:
        args['vectors'] = [80.088, 19.899, 13.383]
        print('using default vectors of {}'.format(args['vectors']))
    axis_conf = {'x':0,'y':1,'z':2}

    colors = []
    for f, ref in zip(args['files'],args['reference']):
        coordsA = reader.xyz(f)
        xyz_data = getCoords(coordsA, args['bead'])
        for axis in args['offAxis']:
            xyz_data = np.array(xyz_data, dtype=np.int8)
            xyz_data[:, axis_conf[axis]] = 0.

        if not args['reference']:
            hist = Hist3D(args['vectors'], args['bins'])
            hist.makeHist(xyz_data)
        else:
            hist = dG3D(args['vectors'], args['bins'])
            hist.dGmap(args['numFrames'], xyz_data, ref, args['Temp'])

        new_file = getFileName(f) + '_bead%s.vtk'%('-'.join(args['bead']))
        vtkRectilinearMesh(new_file, hist.edges, hist.histogram)
        hist.removeDummyValues()
        colors.append(hist.color_values)
    plotColorHist(colors)
