from MCFlow.file_formatting.reader import Movie
from MCFlow.structure_analysis import Struc


class FindRegion(Struc):
    def __init__(self,vectors,channelMap={'straight':{'D'}, 'zig-zag':{'F'}, 'intersection':{'H'}}):
        self.boxLengths = vectors
        self.gridName = '/Users/dejacor/Documents/butanol-water/323K/MFI/energy_grid_114_O.out'
        self.mergedGridName = '/Users/dejacor/Documents/butanol-water/323K/MFI/energy-grids-merged/new_channels-4.0.xyz'
        self.channelMap = channelMap
        self.defaultLabel = 'other'
        self.initGrids()

    def initGrids(self):
        self.missedCoords = False
        self.missedLocations = []
        self.regionNames = []
        for channelName in self.channelMap.keys():
            self.regionNames += list(self.channelMap[channelName])
        self.regionNames += [self.defaultLabel]
        with open(self.gridName) as g:
            line = next(g)
            self.numBins = tuple(map(int, line.split()))
            self.labels = np.chararray(shape=self.numBins, unicode=True, itemsize=5)
            for i in range(self.numBins[0]):
                for j in range(self.numBins[1]):
                    for k in range(self.numBins[2]):
                        self.labels[i,j,k] = self.defaultLabel
        with open(self.mergedGridName) as f:
            next(f)
            next(f)
            for line in f:
                binName = line.split()[0]
                xyz = list(map(float, line.split()[1:]))
                i, j, k = self.findBin(xyz)
                if self.labels[i, j, k] == self.defaultLabel:
                    self.labels[i, j, k] = binName
                else:
                    print('Attempt made to label bin twice, quitting')
                    print('Check that unit cell vectors are correct')
                    print(binName, xyz, i,j,k)
                    quit()

    def findBin(self, xyz):
        (x, y, z) = xyz
        (xmax, ymax, zmax) = self.boxLengths
        (nxbin, nybin, nzbin) = self.numBins
        xbin = int(round( x/xmax * nxbin))
        ybin = int(round( y/ymax * nybin))
        zbin = int(round( z/zmax * nzbin))
        # implement periodic boundaries on bins
        if xbin == nxbin: xbin = 0
        if ybin == nybin: ybin = 0
        if zbin == nzbin: zbin = 0
        return xbin, ybin, zbin

    def getRegion(self, xyz):
        (i, j, k) = self.findBin(xyz)
        return self.labels[i,j,k]

    def getChannel(self, xyz):
        region = self.getRegion(xyz)
        for channel, vals in self.channelMap.items():
            if region in vals:
                return channel
        assert True, 'channel not found'


class Channel(FindRegion):

    def __init__(self, regionMap):
        my_parser = Results()
        my_parser.parser.add_argument('-ID','--name',help='Name of db for molecule number counting',
                               type=str)
        my_parser.parser.add_argument('-B','--bead',help='bead for structure analysis',default=['COM'],
                                     type = str, nargs = '+')
        my_parser.parser.add_argument('-abc','--vectors', help='unit cell vectors. For folding coordinates into a unit cell',
                                      type = float, nargs = '+', default = [20.022,19.899,13.383])
        my_parser.parser.add_argument('-g','--gridName', help='energy grid file name',
                        default = '/Users/dejacor/Documents/butanol-water'
                                  '/323K/MFI/energy_grid_114_O.out', type=str)
        my_parser.parser.add_argument('-M', '--mergedGridName', help='energy grid merged file name',
                        default = '/Users/dejacor/Documents/butanol-water/323K/MFI/'
                                    'energy-grids-merged/new_channels.xyz',type=str)
        my_args = vars(my_parser.parse_args())
        self.args = my_args
        self.defaultLabel = 'other'
        self.boxLengths = self.args['vectors']
        self.mergedGridName = self.args['mergedGridName']
        self.numIndep = self.args['indep']
        self.gridName = self.args['gridName']
        self.channelMap = regionMap
        self.checks()
        self.initGrids()
        self.count = {name:[0 for i in self.numIndep]
                            for name in self.regionNames}
        self.data = {region:[] for region in self.regionNames}
        self.averages = {region:{} for region in self.regionNames}
        self.N = dict()

    def checks(self):
        assert self.args['name'], 'Output ID name needed to output local structure info'
        assert 'box' in self.args['box'], 'Need to specify box for structure analysis'
        self.analysis_class = Movie

    def findBinEnergy(self, bins):
        (i, j, k) = bins
        with open(self.gridName) as g:
            for line in g:
                (xbin, ybin, zbin) = list(map(int, line.split()[:3]))
                if xbin == i and ybin == j and zbin == k:
                    print('bin {} with energy {} not labeled as basin'
                          ' in merged grids'.format(bins, line.split()[-1]))
                    return

    def __str__(self):
        rep = 'REGION      LOADING (mlcl)\n'
        for realRegion in self.channelMap.keys():
            rep += '%12s%12.4f +/-%12.4f\n'%(realRegion, self.averages[realRegion]['mean'],
                                            self.averages[realRegion]['stderr'])
            for subRegion in self.channelMap[realRegion]:
                rep += '    %12s%12.4f +/-%12.4f\n'%(subRegion, self.averages[subRegion]['mean'],
                                            self.averages[subRegion]['stderr'])
        if self.defaultLabel in self.averages.keys():
            rep += '%12s%12.4f +/-%12.4f\n'%(self.defaultLabel, self.averages[self.defaultLabel]['mean'],
                                            self.averages[self.defaultLabel]['stderr'])
        return rep

    def __repr__(self):
        Rep = 'REGION      LOADING (mlcl)\n'
        for realRegion in self.channelMap.keys():
            Rep += '   %-9s & %s\n'%(realRegion, tabularNumbers(self.averages[realRegion]['mean'],
                                                      self.averages[realRegion]['stderr']))
            for subRegion in self.channelMap[realRegion]:
                Rep += '       %-9s & %s\n'%(subRegion, tabularNumbers(self.averages[subRegion]['mean'],
                                                      self.averages[subRegion]['stderr']))
        if self.defaultLabel in self.averages.keys():
            Rep += '   %-9s & %s\n'%(self.defaultLabel, tabularNumbers(self.averages[self.defaultLabel]['mean'],
                                        self.averages[self.defaultLabel]['stderr']))
        return Rep

    def regionToChannel(self):
        for mol in self.N.keys():
            for channel in self.channelMap.keys():
                self.N[mol][channel] = np.zeros(len(self.numIndep))
                for region in self.channelMap[channel]:
                    values = self.N[mol].pop(region)
                    self.N[mol][channel] = self.N[mol][channel] + np.array(values)

    def avgIndep(self):
        self.averages[self.feed] = {}
        for mol in self.N.keys():
            self.averages[self.feed][mol] = {}
            N_total = 0.
            for channel, values in self.N[mol].items():
                mean = np.mean(values)
                self.averages[self.feed][mol][channel] = {}
                self.averages[self.feed][mol][channel]['mean'] = mean
                self.averages[self.feed][mol][channel]['stdev'] = np.std(values)
                self.averages[self.feed][mol][channel]['raw'] = values
                N_total = N_total + mean
            if self.averages[self.feed][mol][self.defaultLabel]['mean'] == 0.:
                self.averages[self.feed][mol].pop(self.defaultLabel)
            if N_total < 1e-8:
                self.averages[self.feed].pop(mol)
            else:
                print(mol, self.averages[self.feed][mol])

    def normalize(self, frame_seeds):
        indepRange = self.args['indep']
        frame_by_seed = [frame_seeds.count(i) for i in indepRange]
        for mol in self.N.keys():
            for region in self.N[mol].keys():
                for i, val in enumerate(frame_by_seed):
                    self.N[mol][region][i] = self.N[mol][region][i]/val

    def countBins(self, frame_data, frame_seeds):
        box, indepRange = self.args['box'], self.args['indep']
        total_frames = -1
        for FRAME_DATA in frame_data:
            if total_frames == -1:
                print(box)
                N = {mol:{region:[0 for i in indepRange] 
                                    for region in self.regionNames}
                                        for mol in FRAME_DATA[box].keys()}
            total_frames += 1
            seed_index = frame_seeds[total_frames] - 1
            for mol in FRAME_DATA[box].keys():
                for each_molecule in FRAME_DATA[box][mol]:
                    beads = list(set(each_molecule.keys()) &
                                set(self.args['bead']))
                    assert len(beads) == 1, 'Ambiguous beads'
                    bead = beads[0]
                    for each_coord in each_molecule[bead]:
                        (i, j, k) = self.findBin(each_coord)
                        region = self.labels[i, j, k]
                        if region == self.defaultLabel:
                            self.missedLocations.append(each_coord)
                        N[mol][region][seed_index] += 1
        self.N = N
        self.normalize(frame_seeds)
        self.regionToChannel()
        self.avgIndep()

    def output_results(self):
        from runAnalyzer import findNextRun
        import time
        import json
        import os
        next_run = findNextRun('%s/%s/1' % (self.args['path'], self.feed), self.args['type'])
        data_to_save = {
            '%s%i' % (self.args['type'], next_run - 1): self.averages[self.feed],
            'time': time.time()
        }
        file_name = self.args['path'] + '/%s-%s-data.json' % (self.args['name'], self.feed)
        if os.path.isfile(file_name):
            os.rename(file_name, self.args['path'] + '/old-%s-%s-data.json' % (self.args['name'], self.feed))

        with open(file_name, 'w') as f:
            json.dump(data_to_save, f)

    def myCalcs(self, D ):
        if self.args['vectors']: D.foldMovieToUC(self.args['vectors'])
        self.countBins(D.frame_data, D.frame_seed)
        if self.missedLocations:
            data = {'coords':self.missedLocations,'atoms':['N' for i in self.missedLocations]}
            xyz(self.args['path'] + '/'+ self.feed +'/'+ 'gridMissedLocs.xyz', data)
        self.output_results()

    def main(self):
        for feed in self.args['feeds']:
            self.feed = feed
            if self.args['verbosity'] > 0: print('-'*12 + 'Dir is %s'%self.feed + '-'*12)
            analysis = self.read_movies()
            self.myCalcs(analysis)

from MCFlow.parser import Results
import numpy as np
from MCFlow.file_formatting.writer import xyz

if __name__ == '__main__':
#   mapToExp = {'straight':{'H','Li'}, 'zig-zag':{'Be','B','He','D'}}
    mapToExp = {'straight':{'D'}, 'zig-zag':{'F'}, 'intersection':{'H'}}
    M = Channel(mapToExp)
    M.main()
