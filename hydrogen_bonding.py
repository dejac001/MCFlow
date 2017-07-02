
from MCFlow.structure_analysis import Struc
aOHO_min = 120./180*3.1415926535897931 # degrees
rOH_max = 6.25 # 2.5*2.5
rOO_max = 10.89 # 3.3*3.3

def findHydroxylHydrogen(Oxyz, Hcoords, abc):
    for Hxyz in Hcoords:
        if calculate_distance2(Oxyz, Hxyz, abc) < 1.0:
            return Hxyz

def hy_bond_from_DB(hmap_data, box, htype, name, feed):
    run = 'prod-'
    my_db = {run:{}}
    for pair, values in hmap_data[box].items():
        mol1, mol2 = [i.strip('mol') for i in pair.split('-')]
        pair1 = mol1 + '->' + mol2
        pair2 = mol2 + '->' + mol1
        nHB = 0
        for rOH, aOHO in zip(values['distance'], values['angle']):
            if htype == 'strict':
                if (rOH < rOH_max) and (aOHO > aOHO_min):
                    nHB += 1
        n1 = hmap_data[box][mol1]
        n2 = hmap_data[box][mol2]
        my_db[run][pair1] = {box:{'mean':nHB/n1,'stdev':0.}}
        my_db[run][pair2] = {box:{'mean':nHB/n2,'stdev':0.}}
    with shelve.open(name +'-'+ htype + '-'+box+'-data.db',writeback=True) as db:
        db[feed] = my_db

def findHB(beadsFrom, beadsTo, abc, criteria):
    nHB = 0
    for O1 in beadsFrom['O']:
        for O2 in [j for j in beadsTo['O'] if calculate_distance2(O1,j,abc) > 0.1]:
            # iterate through all H on molecule 1 that are bonded to O1
            for H1 in [i for i in beadsFrom['H'] if calculate_distance2(O1,i,abc) < 1.0]:
                # look for O1--H1...O2 hbonds
                if criteria == 'loose':
                    if ((calculate_distance2(O1, O2, abc) < rOO_max) and
                        (calculate_distance2(O2, H1, abc) < rOH_max)):
                        nHB += 1
                elif criteria == 'strict':
                    if ((calculate_distance2(H1,O2,abc) < rOH_max) and
                        (calculate_angle(O1,H1,O2,abc) > aOHO_min)):
                            nHB +=1
            # iterate through all H on molecule 2 that are bonded to O2
            for H2 in [i for i in beadsTo['H'] if calculate_distance2(O2,i,abc) < 1.0]:
                # look for O2--H2...O1 hbonds
                if criteria == 'loose':
                    if ((calculate_distance2(O1, O2, abc) < rOO_max) and
                        (calculate_distance2(O1, H2, abc) < rOH_max)):
                        nHB += 1
                elif criteria == 'strict':
                    if ((calculate_distance2(H2,O1,abc) < rOH_max) and
                        (calculate_angle(O2,H2,O1,abc) > aOHO_min)):
                            nHB +=1
    return nHB

def readDBs(path, my_feed, my_type, boxes):
    my_hist_data = {}
    with shelve.open('%s/%s/HB-map.db'%(path,my_feed)) as db:
        if my_feed not in db.keys(): return my_hist_data
        for key, value in db[my_feed].items():
            if my_type in key:
                my_hist_data = value
                if 'nchain count' not in value.keys():
                    my_hist_data['nchain count'] = {}
                for box in boxes:
                    if box not in value.keys():
                        # we need to search other dbs, its not in this one
                        print(box, 'not in db for my_feed',my_feed)
                        my_hist_data['nchain count'][box] = {}
                        sims = [i for i in os.listdir(my_feed)
                                    if (os.path.isdir( my_feed +'/' + i) and i.isdigit())]
                        my_hist_data[box] = {}
                        for sim in sims:
                            with shelve.open('%s/%s/%s/HB-map.db'%(path,my_feed,sim)) as db2:
                                for k2, v2 in db2[my_feed].items():
                                    if my_type in k2:
                                        if box in v2.keys():
                                            for pair, val in v2[box].items():
                                                if pair not in my_hist_data[box].keys():
                                                    my_hist_data[box][pair] = {'distance':[],'angle':[]}
                                                my_hist_data[box][pair]['distance'] += val['distance']
                                                my_hist_data[box][pair]['angle'] += val['angle']
                                        if ('nchain count' in v2.keys()):
                                            for molty, m_count in v2['nchain count'][box].items():
                                                if molty not in my_hist_data['nchain count'][box].keys():
                                                    my_hist_data['nchain count'][box][molty] = m_count
                                                else:
                                                    my_hist_data['nchain count'][box][molty] += m_count
                    elif box not in my_hist_data['nchain count'].keys():
                        my_hist_data['nchain count'][box] = {}
    return my_hist_data

from MCFlow.file_formatting.reader import Movie

class HydrogenBond(Movie):
    def __init__(self, file_name, *args):
        Movie.__init__(self, file_name, *args)


    def getBeads(self, my_box):
        '''
        '''
        self.HB_info = []
        OTypes = {'62':'alkanol','114':'water','178':'zeo','181':'silanol'}
        HTypes ={'61':'alkanol','115':'water','182':'silanol'}
        for iframe, FRAME_DATA in enumerate(self.frame_data):
            # store H and O for all other mols
            HB_mols = {}
            for molType in FRAME_DATA[my_box].keys():
                if molType not in HB_mols.keys():
                    HB_mols[molType] = []
                for imol, each_molecule in enumerate(FRAME_DATA[my_box][molType]):
                    beads = {'H':[],'O':[]}
                    for bead in each_molecule.keys():
                        if bead in OTypes.keys():
                            for coord in each_molecule[bead]:
                                beads['O'].append( list(map(float,coord)) )
                        elif bead in HTypes.keys():
                            for coord in each_molecule[bead]:
                                beads['H'].append( list(map(float,coord)) )
                    HB_mols[molType].append( beads )
            self.HB_info.append( HB_mols )

    def calcHB(self, my_box):
        self.getBeads(my_box)
        self.HB = []
        for iframe, HB_data in enumerate(self.HB_info):
            try:
                if (iframe+1)%(self.nframes//4) == 0:
                    print('%5.1f %% of frames analyzed for HB'%(100*(iframe+1)/self.nframes))
            except ZeroDivisionError:
                print('%5.1f %% of frames analyzed for HB'%(100*(iframe+1)/self.nframes))
            self.HB.append( {my_box:{}  } )
            for mol1 in sorted(HB_data.keys()):
                for mol2 in sorted(HB_data.keys()):
                    pair = mol1.strip('mol') + '->' + mol2.strip('mol')
                    self.HB[iframe][my_box][pair] = []
                    for HB_from_beads in HB_data[mol1]:
                        my_from_HB = 0
                        for HB_to_beads in HB_data[mol2]:
                            # determine hydrogen bonds with molecules of interest
                            my_nHB = findHB(HB_from_beads, HB_to_beads, self.boxlengths[iframe][my_box], self.criteria)
                            my_from_HB += my_nHB
                        self.HB[iframe][my_box][pair].append( my_from_HB )

    def filterHB(self):
        '''
        we want to only keep positions of O and H involved in a hydrogen bond
        '''
        for iframe, FRAME_DATA in enumerate(self.frame_data):
            for box, values in self.HB[iframe].items():
                for pair in values.keys():
                    indices_to_keep = []
                    for ipair, nHB in enumerate(values[pair]):
                        if nHB > 0:
                            indices_to_keep.append( ipair )
                    FRAME_DATA[box][pair] = [value for i, value in enumerate(FRAME_DATA[box][pair])
                                                   if i in indices_to_keep]
    def countHB(self, nIndep, feed, data):
        new = self.__class__('HB-data')
        for attr, value in self.__dict__.items():
            new.__dict__[attr] = value
        new.countMols(nIndep, feed, data)
        return new

class HB(Struc):
    def __init__(self):
        self.parser = MultMols()
        #TODO: make able to do multiple boxes at same time
        self.parser.parser.add_argument('-ID','--name',help='Name of db for molecule number counting',
                               type=str,default = '')
        self.parser.parser.add_argument('-H','--htype',help='hydrogen bonding criteria type',
                                      type=str, choices = ['loose','strict'],default='strict')
        self.getArgs()

    def getArgs(self):
        my_args = vars(self.parser.parse_args())
        self.args = my_args
        if 'box' not in self.args['box']:
            self.args['box'] = 'box%s'%self.args['box']
        self.checks()

    def checks(self):
        assert self.args['name'], 'Output ID name needed to output local structure info'
        self.analysis_class = HydrogenBond

    def myCalcs(self, D ):
        D.criteria = self.args['htype']
        hist_data = readDBs(self.args['path'], self.feed, self.args['type'], [self.args['box']])
        if (('nchain count' in hist_data.keys()) and
            (self.args['box'] in hist_data['nchain count'].keys()) and
            (len(hist_data['nchain count'][self.args['box']].keys()) > 0)):
            print('get HB rfrom db')
            # make HB from database
            HB = hy_bond_fromDB(hist_data, self.args['box'], self.args['htype'],
                            self.args['name'], self.feed)
            # output DB somehow
        else:
            D.calcHB(self.args['box'])
            HB = D.countHB(self.args['indep'],self.feed, D.HB)
            outputDB(self.args['path'],[self.feed],self.args['type'],
                {self.args['name']+'-' + self.args['htype'] + '-'+self.args['box']:HB})

from MCFlow.runAnalyzer import what2Analyze
from MCFlow.file_formatting.writer import xyz
from MCFlow.calc_tools import calculate_distance2, calculate_angle
from MCFlow.parser import MultMols
from MCFlow.getData import outputDB
import shelve, os

if __name__ == '__main__':
    M = HB()
    M.main()
