from structure_analysis import Struc, output_json
rOO_max = 10.89 # 3.3*3.3


def findHydroxylHydrogen(Oxyz, Hcoords, abc):
    for Hxyz in Hcoords:
        if calculate_distance2(Oxyz, Hxyz, abc) < 1.0:
            return Hxyz


def hy_bond_from_DB(hmap_data, box, htype, name, feed):
    run = 'prod-'
    my_data = {run:{}}
    for pair, values in hmap_data[box].items():
        mol1, mol2 = [i.strip('mol') for i in pair.split('-')]
        pair1 = mol1 + '->' + mol2
        pair2 = mol2 + '->' + mol1
        nHB = 0
        for rOH, aOHO in zip(values['distance'], values['angle']):
            if htype == 'strict':
                if (rOH < rOH_max) and (aOHO > angle):
                    nHB += 1
        n1 = hmap_data[box][mol1]
        n2 = hmap_data[box][mol2]
        my_data[run][pair1] = {box:{'mean':nHB/n1,'stdev':0.}}
        my_data[run][pair2] = {box:{'mean':nHB/n2,'stdev':0.}}
    with open(name +'-'+ htype + '-'+box+'='+ feed + '-data.json', 'w') as f:
        json.dump(my_data, f)

def findHB(beadsFrom, beadsTo, abc, criteria, angle, dist_sq, molNum1, molNum2):
    def addData(molDonor, molAcceptor, O_Donor, O_Acceptor):
        HB_pairs['mol donor'].append(molDonor)
        HB_pairs['oxygen donor'].append(O_Donor)
        HB_pairs['mol acceptor'].append(molAcceptor)
        HB_pairs['oxygen acceptor'].append(O_Acceptor)
    nHB = 0
    HB_pairs = {'mol donor':[],
                         'mol acceptor':[],
                         'oxygen donor':[],
                         'oxygen acceptor':[]}
    for nO1, O1 in enumerate(beadsFrom['O']):
        for nO2, O2 in enumerate(beadsTo['O']):
            if calculate_distance2(O1,O2,abc) < 0.1: continue
            # iterate through all H on molecule 1 that are bonded to O1
            for H1 in [i for i in beadsFrom['H'] if calculate_distance2(O1,i,abc) < 1.0]:
                # look for O1--H1...O2 hbonds
                if criteria == 'loose':
                    if ((calculate_distance2(O1, O2, abc) < rOO_max) and
                        (calculate_distance2(O2, H1, abc) < dist_sq)):
                        nHB += 1
                        addData(molNum2, molNum1, nO2+1, nO1+1)
                elif criteria == 'strict':
                    if ((calculate_distance2(H1,O2,abc) < dist_sq) and
                        (calculate_angle(O1,H1,O2,abc) > angle)):
                        nHB +=1
                        addData(molNum2, molNum1, nO2+1, nO1+1)
            # iterate through all H on molecule 2 that are bonded to O2
            for H2 in [i for i in beadsTo['H'] if calculate_distance2(O2,i,abc) < 1.0]:
                # look for O2--H2...O1 hbonds
                if criteria == 'loose':
                    if ((calculate_distance2(O1, O2, abc) < rOO_max) and
                        (calculate_distance2(O1, H2, abc) < dist_sq)):
                        nHB += 1
                        addData(molNum1, molNum2, nO1+1, nO2+1)
                elif criteria == 'strict':
                    if ((calculate_distance2(H2,O1,abc) < dist_sq) and
                        (calculate_angle(O2,H2,O1,abc) > angle)):
                        nHB +=1
                        addData(molNum1, molNum2, nO1+1, nO2+1)
    return nHB, HB_pairs


def read_json(path, my_feed, my_type, boxes):
    """ read databases

    """
    my_hist_data = {}
    with open('%s/%s/HB-map.json'%(path,my_feed)) as f:
        db = json.load(f)
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
                        with open('%s/%s/%s/HB-map.json'%(path,my_feed,sim)) as f:
                            db2 = json.load(f)
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


from mcflow.file_formatting.reader import Movie


class HydrogenBond(Movie):
    def __init__(self, file_name, *args):
        Movie.__init__(self, file_name, *args)

    def getBeads(self, my_box):
        '''
        '''
        self.HB_info = []
         #OTypes = {'62':'alkanol','114':'water','178':'zeo','181':'silanol'}
        # HTypes ={'61':'alkanol','115':'water','182':'silanol'}
        OTypes = {'62':'alkanol','114':'water','178':'silanol','111':'furan','302':'acidO=C','303':'acidO-C'} # '502':zeo
        HTypes ={'61':'alkanol','115':'water','42':'silanol','304':'acidH'}
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

    def calcHB(self, my_box, verbosity):
        self.getBeads(my_box)
        self.HB = []
        self.HB_chains = []
        for iframe, HB_data in enumerate(self.HB_info):
            if verbosity > 0:
                try:
                    if (iframe+1)%(self.nframes//4) == 0:
                        print('%5.1f %% of frames analyzed for HB'%(100*(iframe+1)/self.nframes))
                except ZeroDivisionError:
                    print('%5.1f %% of frames analyzed for HB'%(100*(iframe+1)/self.nframes))
            self.HB.append( {my_box:{}  } )
            self.HB_chains.append( {my_box:{}} )
            self.HB_pairs = {}
            for mol1 in sorted(HB_data.keys()):
                for mol2 in sorted(HB_data.keys()):
                    pair = mol1.strip('mol') + '->' + mol2.strip('mol')
                    # if pair in self.HB[iframe][my_box].keys(): continue
                    self.HB[iframe][my_box][pair] = []
                    self.HB_pairs[pair] = {
                        'mol donor':[],
                         'mol acceptor':[],
                         'oxygen donor':[],
                         'oxygen acceptor':[]
                                     }
                    for n1, HB_from_beads in enumerate(HB_data[mol1]):
                        my_from_HB = 0
                        for n2, HB_to_beads in enumerate(HB_data[mol2]):
                            # determine hydrogen bonds with molecules of interest
                            my_nHB, HB_pairs = findHB(HB_from_beads, HB_to_beads, 
                                                    self.boxlengths[iframe][my_box], self.criteria,
                                                    self.angle_radians, self.dist_squared,
                                            n1+1, n2+1)
                            my_from_HB += my_nHB
                            for key, value in HB_pairs.items():
                                self.HB_pairs[pair][key] += value
                        self.HB[iframe][my_box][pair].append( my_from_HB )
            for pair, values in self.HB_pairs.items():
                if sum(self.HB[iframe][my_box][pair]) == 0.:
                    self.HB[iframe][my_box][pair+' oxygen path length'] = []
                    self.HB[iframe][my_box][pair+' cluster size'] = []
                    continue
                mol1, mol2 = pair.split('->')
                if mol1 != mol2 or mol1 == '1': continue
                nMol = len(HB_data['mol'+mol1])
                # make distance matrix for pair
                G = nx.Graph()
                for md, ma, od, oa in zip(values['mol donor'],values['mol acceptor'],
                                          values['oxygen donor'], values['oxygen acceptor']):
                    i = md-1 + nMol*(od-1)
                    j = ma-1 + nMol*(oa-1)
                    G.add_edge(i,j)
                    G.add_edge(j,i)
                for I in range(nMol):
                    G.add_edge(I, I+nMol)
                    G.add_edge(I+nMol, I)
                graphs = list(nx.connected_component_subgraphs(G))
#               path_lengths = []
#               for m in graphs:
#                   end_nodes = [i for i in m if m.degree(i)==1]
#                   for i, source in enumerate(end_nodes):
#                       for j in range(i+1,len(end_nodes)):
#                           my_path = nx.shortest_path(m,source=source,
#                                                      target=end_nodes[j])
#                           path_lengths.append(len(my_path))
#               self.HB[iframe][my_box][pair+' oxygen path length'] = path_lengths
                self.HB[iframe][my_box][pair+' cluster size'] = [len(m.nodes()) for
                                                                 m in graphs]

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
        self.parser.parser.add_argument('-R','--readDB',help='whether or not to read databases',
                                      type=bool, default=False)
        self.parser.parser.add_argument('-a','--minAngle',help='minimum OHO angle (degrees) for strict criteria',
                                            type=float,default=120)
        self.parser.parser.add_argument('-d','--minDist',help='minimum OH dist (angstrom) for strict criteria',
                                            type=float,default=2.5)
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
        D.angle_radians, D.dist_squared = self.args['minAngle']/180.*np.pi, self.args['minDist']*self.args['minDist']
        if self.args['readDB']:
            hist_data = read_json(self.args['path'], self.feed, self.args['type'], [self.args['box']])
        if self.args['readDB'] and (('nchain count' in hist_data.keys()) and
                (self.args['box'] in hist_data['nchain count'].keys()) and
                (len(hist_data['nchain count'][self.args['box']].keys()) > 0)):
            raise NotImplementedError
        else:
            D.calcHB(self.args['box'], self.args['verbosity'])
            HB = D.countHB(self.args['indep'],self.feed, D.HB)
            output_json(self.args['path'],self.args['type'],
                {self.feed: HB}, 'Yes')


nx_options = {
 'node_color': 'black',
'node_size': 8,
'width': 1}

import networkx as nx
from mcflow.calc_tools import calculate_distance2, calculate_angle
from analysis_parsers import MultMols
import numpy as np
import json, os

if __name__ == '__main__':
    M = HB()
    M.main()
