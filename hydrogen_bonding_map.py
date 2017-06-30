from MCFlow.hydrogen_bonding import HydrogenBond, HB

class HB_map(HydrogenBond):
    def __init__(self, file_name, *args):
        HydrogenBond.__init__(self, file_name, *args)

    def makeMap(self, box, numFrames):
        assert numFrames < self.nframes, 'Not enough frames'
        if numFrames <= 0: numFrames = self.nframes
        self.getBeads( box)
        my_box = 'box%s'%box
        self.histogram = {}
        for iframe, HB_data in enumerate(self.HB_info):
            abc = self.boxlengths[iframe][my_box]
            try:
                if (iframe+1)%(numFrames//4) == 0:
                    print('%5.1f %% of frames analyzed for HB'%(100*(iframe+1)/numFrames))
            except ZeroDivisionError:
                print('%5.1f %% of frames analyzed for HB'%(100*(iframe+1)/numFrames))
            if iframe + 1 > numFrames:
                break
            for molType1 in HB_data.keys():
                for molType2 in HB_data.keys():
                    mols_int = sorted(map(int,[molType1.strip('mol'), molType2.strip('mol')]))
                    mols = ['mol%i'%i for i in mols_int]
                    pair = '-'.join(mols)
                    if pair not in self.histogram.keys():
                        self.histogram[pair] = {'distance':[],'angle':[]}
                    for i1, mol1 in enumerate(HB_data[molType1]):
                        for i2, mol2 in enumerate(HB_data[molType2]):
                            # mol1 and mol2 to are {'H':[[x1,y1,z1]...],'O':[[x1,y1,z1]...]}
                            for O1 in mol1['O']:
                                for H1 in [i for i in mol1['H'] if calculate_distance2(O1,i,abc) < 1.0]:
                                    for O2 in mol2['O']:
                                        rOO = calculate_distance(O1, O2, abc)
                                        rOH = calculate_distance(H1, O2, abc)
                                        if (rOH <= r_max) and (rOO > 0.01):
                                            aOHO = calculate_angle(O1, H1, O2, abc)/np.pi*180.
                                            self.histogram[pair]['distance'].append( rOH )
                                            self.histogram[pair]['angle'].append( aOHO )

    def storeHist(self, feed, path, type, box, indep):
        import time
        nextRun = runAnalyzer.findNextRun('%s/1/'%path, type)
        if len(indep) == 1:
            path = '%s/%i/'%(path,indep[0])
        with shelve.open(path + '/HB-map.db',writeback=True) as db:
            if feed not in list(db.keys()):
                db[feed] = {}
            run = '%s%i'%(type, nextRun-1)
            if run not in db[feed].keys():
                db[feed][run] = {}
            db[feed][run]['box%s'%box] = self.histogram
            db[feed]['time'] = time.time()

def plotHist(histogram, box, path, rMax, aMin):
    def plot(x_data,y_data, x_edges, y_edges):
        fig = plt.figure()
        ax = fig.add_subplot(111)

    xmin_plot = 1.0
    xmax_plot = rMax # maximum radius
    ymin_plot = aMin # minimum angle
    ymax_plot = 180.
    xedges = np.linspace(xmin_plot, xmax_plot,100)
    yedges = np.linspace(ymin_plot, ymax_plot,100)
    all_x = []; all_y = []
    for pair in histogram.keys():
        x = []; y = []
        for (i,j) in zip(histogram[pair]['distance'], histogram[pair]['angle']):
            if (i <= x_max_plot) and (j >= ymin_plot):
                x.append(i)
                y.append(j)
        all_x += x
        all_y += y
#       if len(x) > 0:
#           fig = plt.figure()
#           ax = fig.add_subplot(111)
#           H, xedges, yedges = np.histogram2d(x,y,bins=(xedges,yedges),normed=False)
#           H = H/np.max(H)
#           image =  ax.contour(H.T, extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]],
#                               linewidths=3, levels=[0.2, 0.4, 0.8, 0.9, 1.0])
#           hist, xedges, yedges, image = ax.hist2d(x, y, normed=False,bins=(xbins,ybins))
            # format Axis
#           plt.subplots_adjust(left=0.15,right=0.96,top = 0.98, bottom=0.02)
#           cbar =fig.colorbar(image, orientation='horizontal',pad=0.13)#,fraction = 0.046)
#           cbar.set_label('$h(r,\\angle)$',fontdict=font)
#           cbar.ax.tick_params(labelsize=14)
#           fig.set_size_inches(5.0, 5.0)
#           fig.savefig('%s/HB-map-box%s_%s_rOH%2.1f.png'%(path,box,pair,r_max), dpi=600)
    pair = 'all'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    H, xedges, yedges = np.histogram2d(all_x,all_y,bins=(xedges,yedges),normed=False)
    H = H/np.max(H)
    image =  ax.contour(H.T, extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]],
                                linewidths=3, levels=[0.2, 0.4, 0.8, 0.9, 1.0])
    plt.subplots_adjust(left=0.15,right=0.96,top = 0.98, bottom=0.02)
    cbar =fig.colorbar(image, orientation='horizontal',pad=0.13)#,fraction = 0.046)
    cbar.set_label('$h(r,\\angle)$',fontdict=font)
    cbar.ax.tick_params(labelsize=14)
    fig.set_size_inches(5.0, 5.0)
    fig.savefig('%s/HB-map-box%s_%s_rOH%2.1f.png'%(path,box,pair,r_max), dpi=600)
#   plt.show()


class HB_format_map(HB):

    def __init__(self):
        HB.__init__(self)

    def getArgs(self):
        self.parser.parser.add_argument('-n','--numFrames',help='number of frames to analyze',type=int,default=0)
        my_args = vars(self.parser.parse_args())
        self.args = my_args
        self.checks()

    def checks(self):
        self.analysis_class = HB_map
        assert self.args['box'], 'Box needed for hydrogen bonding'

    def myCalcs(self, H):
        H.makeMap(self.args['box'], self.args['numFrames'])
        directory = self.args['path'] + '/' + self.feed
        H.storeHist(self.feed, directory, self.args['type'], self.args['box'], self.args['indep'])
        if os.name != 'posix':
            H.plotHist(self.args['box'], directory)

from MCFlow.calc_tools import calculate_distance, calculate_angle, calculate_distance2
import numpy as np
import shelve
from MCFlow import runAnalyzer

font = {'size':18}
r_max = 5.0
import os
if os.name != 'posix':
    import matplotlib.pyplot as plt
    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True)

if __name__ == '__main__':
    H = HB_format_map()
    H.main()
