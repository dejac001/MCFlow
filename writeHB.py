from MCFlow.writeXvY import LoadAds, MoleFrac, LiqAds
class HB:
    def __init__(self):
        pass

    def readHB(self):
        assert 'box' in self.box, 'Box needed for HB plot'
        self.HB = {}
        for feed in self.feeds:
            my_dir = '%s/%s'%(self.path,feed)
            files_list = []
            for i in os.listdir(my_dir):
                if i.startswith('HB-%sdeg-%sAngst-'%(self.angle, self.dist)) and (self.box in i):
                    if i[-5:] == '.json':
                        files_list.append(i)
            files_list = set(files_list)
            for file in files_list:
                with open('%s/%s'%(my_dir, file)) as f:
                    data = json.load(f)
                if (feed in data.keys()) and (len(data[feed].keys()) > 1):
                    self.HB[feed] = data[feed]
        not_in_json = []
        for feed in self.feeds:
            if feed not in self.HB.keys():
                not_in_json.append(feed)
        assert len(not_in_json) == 0, 'Feeds {} not in databases'.format(not_in_json)
        for i in not_in_json:
            print('not plotting for feed %s, %sdeg-%sAngst'%(i,self.angle,self.dist))
            self.feeds.remove( i )

    def HB_write(self):
        assert self.mol, 'Mol needed  for HB plot x axis'
        nIndep = self.gen_data[self.feed][self.run]['numIndep']
        my_run = [i for i in self.HB[self.feed].keys() if i[:-1] in self.run][0]
        HB = self.HB[self.feed][my_run]
        file_description = 'Q       N_HB        dQ          dN_HB'
        if self.box == 'box3':
            self.box = 'box2'
            X = self.getX()
            self.box = 'box3'
        else:
            X = self.getX()
        for pair in HB.keys():
            if type(HB[pair]) == type(int(1)):
                continue
            elif HB[pair][self.box]['mean'] > 0.:
                file_name = 'HB_%s_v_%s_%s.dat'%(pair.replace('->','-').replace(' ',''), self.mol, self.xlabel[0])
                if '95conf' in X.keys():
                    err = X['95conf']
                    print(X)
                else:
                    err = calc95conf(X['stdev'],nIndep)
                writeAGR([X['mean']], [HB[pair][self.box]['mean']],
                         [err], [calc95conf(HB[pair][self.box]['stdev'], nIndep)],
                         [self.feed], file_name, file_description)

class HBvQ(LoadAds, HB):
    def __init__(self, **kwargs):
        LoadAds.__init__(self, **kwargs)
        assert self.units, 'Units needed for plot'
        if kwargs['feeds']:
            self.feeds = kwargs['feeds']
        if kwargs['angle']:
            self.angle = kwargs['angle']
        if kwargs['dist']:
            self.dist = kwargs['dist']

class HBvC(LiqAds, HB):
    def __init__(self, **kwargs):
        LiqAds.__init__(self, **kwargs)
        if kwargs['feeds']:
            self.feeds = kwargs['feeds']
        if kwargs['angle']:
            self.angle = kwargs['angle']
        if kwargs['dist']:
            self.dist = kwargs['dist']

class HBvX(MoleFrac, HB):
    def __init__(self, **kwargs):
        MoleFrac.__init__(self, **kwargs)
        if kwargs['feeds']:
            self.feeds = kwargs['feeds']
        if kwargs['angle']:
            self.angle = kwargs['angle']
        if kwargs['dist']:
            self.dist = kwargs['dist']
        

def parse_input():
    my_parser = Plot()
    my_parser.axes()
    my_parser.isotherm()
    my_parser.parser.add_argument('-a','--angle',help='angle for json criteria',type=str)
    my_parser.parser.add_argument('-d','--dist',help='distance for json criteria',type=str)

    kwargs = vars(my_parser.parse_args())
    assert kwargs['xaxis'], 'No x axis chosen for plot'
    return kwargs
import shelve, os
from MCFlow.parser import Plot
from MCFlow.runAnalyzer import calc95conf, checkRun
from MCFlow.writeXvY import writeAGR

if __name__ == '__main__':
    # TODO: find way to have conditional arguments based off of what xaxis and yaxis are
    # (in other words, we have an excess of variables as arguments)
    args = parse_input()

    if args['xaxis'] == 'Q':
        my_plotter = HBvQ(**args)
    elif args['xaxis'] == 'x':
        my_plotter = HBvX(**args)
    elif args['xaxis'] == 'C':
        my_plotter = HBvC(**args)
    my_plotter.readDBs()
    my_plotter.readHB()

    for feed in my_plotter.feeds:
        # determine if run has been completed
        my_plotter.run = checkRun(args['type'], my_plotter.variables, feed)
        my_plotter.feed = feed
        my_plotter.HB_write()
