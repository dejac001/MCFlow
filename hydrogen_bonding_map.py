from MCFlow.hydrogen_bonding import HydrogenBond, HB
import json

class HB_map(HydrogenBond):
    def __init__(self, file_name, *args):
        HydrogenBond.__init__(self, file_name, *args)
        self.nchain_count = {}

    def makeMap(self, box, numFrames):
        assert numFrames < self.nframes, 'Not enough frames'
        if numFrames <= 0: numFrames = self.nframes
        if 'box' not in box:
            my_box = 'box%s'%box
        else:
            my_box = box
        self.getBeads( my_box)
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
                if molType1 not in self.nchain_count.keys():
                    self.nchain_count[molType1] = 0
                for molType2 in HB_data.keys():
                    mols_int = sorted(map(int,[molType1.strip('mol'), molType2.strip('mol')]))
                    mols = ['mol%i'%i for i in mols_int]
                    pair = '-'.join(mols)
                    if pair not in self.histogram.keys():
                        self.histogram[pair] = {'distance':[],'angle':[]}
                    for i1, mol1 in enumerate(HB_data[molType1]):
                        self.nchain_count[molType1] += 1
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
        db = {}
        run = '%s%i'%(type, nextRun-1)
        db[feed] = {run:  {}}
        db[feed][run]['box%s'%box] = self.histogram
        if 'nchain count' not in db[feed][run].keys():
            db[feed][run]['nchain count'] = {}
        db[feed][run]['nchain count']['box%s'%box] = self.nchain_count
        db[feed]['time'] = time.time()
        with open(path + '/HB-map-%s.json' % feed, 'w') as f:
            json.dump(db, f)

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
        assert self.args['box'].isdigit(), 'Box needed for hydrogen bonding'

    def myCalcs(self, H):
        H.makeMap(self.args['box'], self.args['numFrames'])
        directory = self.args['path'] + '/' + self.feed
        H.storeHist(self.feed, directory, self.args['type'], self.args['box'], self.args['indep'])

from MCFlow.calc_tools import calculate_distance, calculate_angle, calculate_distance2
import numpy as np
import shelve
from MCFlow import runAnalyzer

r_max = 5.0

if __name__ == '__main__':
    H = HB_format_map()
    H.main()
