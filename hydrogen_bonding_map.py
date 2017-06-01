from MCFlow.hydrogen_bonding import HydrogenBond, HB

class HB_map(HydrogenBond):
    def __init__(self, file_name, *args):
        HydrogenBond.__init__(self, file_name, *args)

#    def getBeads(self, box):
#        '''
#        '''
#        my_box = 'box%s'%box
#        self.HB_info = []
#        OTypes = {'62':'alkanol','114':'water','178':'zeo','181':'silanol'}
#        HTypes ={'61':'alkanol','115':'water','182':'silanol'}
#        for iframe, FRAME_DATA in enumerate(self.frame_data):
#            # store H and O for all other mols
#            HB_mols = {}
#            for molType in FRAME_DATA[my_box].keys():
#                if molType not in HB_mols.keys():
#                    HB_mols[molType] = []
#                for imol, each_molecule in enumerate(FRAME_DATA[my_box][molType]):
#                    beads = {'H':[],'O':[]}
#                    for bead in each_molecule.keys():
#                        if bead in OTypes.keys():
#                            for coord in each_molecule[bead]:
#                                beads['O'].append( list(map(float,coord)) )
#                        elif bead in HTypes.keys():
#                            for coord in each_molecule[bead]:
#                                beads['H'].append( list(map(float,coord)) )
#                    HB_mols[molType].append( beads )
#            self.HB_info.append( HB_mols )

    def makeMap(self, box):
        self.getBeads( box)
        my_box = 'box%s'%box
        self.histogram = {}
        for iframe, HB_data in enumerate(self.HB_info):
            abc = self.boxlengths[iframe][my_box]
            try:
                if (iframe+1)%(self.nframes//4) == 0:
                    print('%5.1f %% of frames analyzed for HB'%(100*(iframe+1)/self.nframes))
            except ZeroDivisionError:
                print('%5.1f %% of frames analyzed for HB'%(100*(iframe+1)/self.nframes))
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

    def formatAx(my_ax):
        xlabel, ylabel= ['$r_{\mathrm{HO}}$', '$\\angle~\mathrm{OHO}$']
        my_ax.set_xlabel(xlabel,fontdict=font)
        my_ax.set_ylabel(ylabel,fontdict=font)
        # set y ticks
        y_ticks = [0., 20., 40., 60., 80., 100., 120., 140., 160., 180.]
        my_ax.set_yticks(y_ticks)
        y_minor = [(y_ticks[i] + y_ticks[i+1])/2. for i in range(len(y_ticks)-1)]
        my_ax.set_yticks(y_minor,minor=True)
        my_ax.tick_params(axis='y',direction='out',which='both',labelsize=14,left=True,right=True)
        # set x ticks
        x_ticks = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        my_ax.set_xticks(x_ticks)
        my_ax.tick_params(axis='x',direction='out',which='both',labelsize=14,bottom=True,top=True)
        x_minor = [(x_ticks[i] + x_ticks[i+1])/2. for i in range(len(x_ticks)-1)]
        my_ax.set_xticks(x_minor, minor=True)

    def plotHist(self, box, path):
        xmin_plot = 1.0
        xmax_plot = r_max
        ymin_plot = 0.
        ymax_plot = 180.
        xedges = np.linspace(xmin_plot, xmax_plot,100)
        yedges = np.linspace(ymin_plot, ymax_plot,100)
        all_x = []; all_y = []
        for pair in self.histogram.keys():
            x, y = (self.histogram[pair]['distance'], self.histogram[pair]['angle'])
            all_x += x
            all_y += y
            assert len(x) == len(y), 'x and y lengths not equal'
            if len(x) > 0:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                H, xedges, yedges = np.histogram2d(x,y,bins=(xedges,yedges),normed=False)
                H = H/np.max(H)
                image =  ax.contour(H.T, extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]],
                                    linewidths=3, levels=[0.2, 0.4, 0.8, 0.9, 1.0])
#               hist, xedges, yedges, image = ax.hist2d(x, y, normed=False,bins=(xbins,ybins))
                # format Axis
                plt.subplots_adjust(left=0.15,right=0.96,top = 0.98, bottom=0.02)
                cbar =fig.colorbar(image, orientation='horizontal',pad=0.13)#,fraction = 0.046)
                cbar.set_label('$h(r,\\angle)$',fontdict=font)
                cbar.ax.tick_params(labelsize=14)
                fig.set_size_inches(5.0, 5.0)
                fig.savefig('%s/HB-map-box%s_%s.png'%(path,box,pair), dpi=600)
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
        fig.savefig('%s/HB-map-box%s_%s.png'%(path,box,pair), dpi=600)
#       plt.show()
                

class HB_format_map(HB):

    def __init__(self):
        HB.__init__(self)

    def checks(self):
        self.analysis_class = HB_map
        assert self.args['box'], 'Box needed for hydrogen bonding'

    def myCalcs(self, H):
        H.makeMap(self.args['box'])
        directory = self.args['path'] + '/' + self.feed 
        H.plotHist(self.args['box'], directory)

from MCFlow.calc_tools import calculate_distance, calculate_angle, calculate_distance2
from MCFlow.plotter import formatAxis
import numpy as np
import matplotlib.pyplot as plt

font = {'size':18}
r_max = 4.0
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

if __name__ == '__main__':
    H = HB_format_map()
    H.main()
