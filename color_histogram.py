def plotColorHist(data):
    fig, ax = plt.subplots()
    ibins = np.logspace(0,6,50)
    colors = ['cyan','magenta']
    legend = ['PC','NS']
    lines = ['solid','dashed']
    for i, values in enumerate(data):
        print('number of colored values is %i'%len(values))
        n, bins  = np.histogram(values, bins=ibins, density=False)#, log=True)
        n, new_bins = hist_norm_height(n, bins, max(n))
        ax.step(new_bins, n, color=colors[i],label=legend[i],linestyle=lines[i],linewidth=0.4)
    ax.set_xlabel('ButOH Selectivity',fontsize=8)
    ax.set_ylabel('$h(S)$',fontsize=8)
    ax.legend(loc=0,fontsize='x-small')
    ax.semilogx()
    ax.tick_params(colors='black',width=0.5)
    ax.tick_params(direction='in',which='both',left=True,right=True, labelsize=8)
    ax.tick_params(axis='x',direction='out',bottom=True, top=True,which='both')
    for dirn in ['top','bottom','left','right']:
        ax.spines[dirn].set_color('black')
        ax.spines[dirn].set_linewidth(0.5)
    y_ticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    y_ticks_minor = [0.1,0.3,0.5,0.7,0.9]
    ax.set_yticks(y_ticks)
    ax.set_yticks(y_ticks_minor,minor=True)
    ax.set_xlim([1., 10**6])
    ax.set_ylim([0., 1.])
    # ax.set_ylim([0.9,200.])
    plt.subplots_adjust(left=0.235,right=0.9425, top=0.95,bottom=0.27)
    fig = plt.gcf()
    fig.set_size_inches(1.8,1.54)
    fig.savefig('S_hist.png',dpi=600)
    # plt.show()

import numpy as np
import matplotlib.pyplot as plt
from MCFlow.selectivity_histogram import SHist, hist_norm_height
from MCFlow.file_formatting.writer import vtkRectilinearMesh
from MCFlow.parser import Structure
from MCFlow.file_formatting import reader
from MCFlow.probhistogram import getCoords, getFileName

if __name__ == '__main__':
    my_parser = Structure()
    my_parser.parser.add_argument('-f1','--file1',help='file 2 for selectivity analysis',
                                  type=str,nargs='+')
    my_parser.parser.add_argument('-f2','--file2',help='file 2 for selectivity analysis',
                                  type=str,nargs='+')
    my_parser.parser.add_argument('-off','--offAxis',help='axis to average along',
                                  default=[],nargs='+')
    my_parser.parser.add_argument('-ref','--reference',help='reference density  [for S, Kref]',type=float)
    args = vars(my_parser.parse_args())
    vectors = [[20.022, 19.899, 13.383],[80.088, 19.899, 13.383]]
    args = vars(my_parser.parse_args())
    axis_conf = {'x':0,'y':1,'z':2}

    colors = []
    for f1, f2, abc in zip(args['file1'],args['file2'],vectors):
        coordsA = reader.xyz(f1)
        xyz_dataA = getCoords(coordsA, args['bead'])
        coordsB = reader.xyz(f2)
        xyz_dataB = getCoords(coordsB, args['bead'])
        for axis in args['offAxis']:
            xyz_dataA = np.array(xyz_dataA, dtype=np.int8)
            xyz_dataB = np.array(xyz_dataB, dtype=np.int8)
            xyz_dataA[:, axis_conf[axis]] = 0.
            xyz_dataB[:, axis_conf[axis]] = 0.

        hist = SHist(abc, args['bins'])
        hist.Smap(xyz_dataA, xyz_dataB, args['reference'])
        new_file = getFileName(f1) + getFileName(f2) + '_bead%s.vtk'%('-'.join(args['bead']))
        vtkRectilinearMesh(new_file, hist.edges, hist.histogram)
        hist.removeDummyValues()
        colors.append(hist.color_values)
    plotColorHist(colors)
