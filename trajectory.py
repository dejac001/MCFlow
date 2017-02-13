class Trajectory():
    def __init__(self, boxList, molList, feeds, indepRange):
        self.N_traj = { '%s'%f: {B: { M :{ k: []
                                           for k in indepRange}
                                      for M in molList }
                                 for B in boxList}
                        for f in feeds}
        self.traj_steps =  { '%s'%f: { k: [] for k in indepRange}
                             for f in feeds}
    def addSim(self, N, feed, sim, tag, old_begin, molList, boxList, numSkip):
        offset = int(next(open('%s%i/config.%s%i'%(tag, old_begin, tag, old_begin))))
        if 'prod-' in tag:
            offset = 0
        index = 0
        begin_index = 0
        for index in range(len(N[molList[0]][boxList[0]])):
            if index % numSkip == 0:
                # append average values
                for mol in molList:
                    for box in boxList:
                        self.N_traj['%s'%feed][box][mol][sim].append(
                            np.mean(N[mol][box][begin_index: (index + 1)])
                        )
                self.traj_steps['%s'%feed][sim].append(index+offset)
                begin_index = index

import numpy as np
from plotter import plotTrajectory
if __name__ == '__main__':
    from file_formatting import reader
    import runAnalyzer
    import matplotlib.pyplot as plt
    from parser import Results

    my_parser = Results()
    my_parser.parser.description = 'Plot trajectory from inputting molecules and boxes'
    my_parser.parser.add_argument('-sk','--numSkip',help ='number of MCCs to skip', type = int, default = 100)

    args = vars(my_parser.parse_args())
    assert args['molecules'], 'Mols needed for trajectory'
    assert args['boxes'], 'Boxes needed for trajectory'
    boxes = args['boxes']
    mols = args['molecules']
    assert 'box' in boxes[0], 'Box needed in string'

    for feed in args['feeds']:
        print('starting trajectory for feed {}.....'.format(feed))
        for sim in args['indep']:
            iPath = '%s/%s/%i'%(args['path'], feed,  sim)
            (old_begin, nfiles) = runAnalyzer.what2Analyze(iPath, args['type'], args['guessStart'], args['interval'])
            print('old_begin = {}, nfiles = {}'.format(old_begin, nfiles))
            N, P, boxLengths, ncycle_old, molWeights, U = reader.read_fort12(iPath, old_begin, nfiles, tag=args['type'])
            if feed == args['feeds'][0] and sim == args['indep'][0]:
                # initialize trajectory
                myTraj = Trajectory(boxes, mols, args['feeds'], args['indep'])
            myTraj.addSim(N, feed, sim, args['type'], old_begin, mols, boxes, args['numSkip'])
    plotTrajectory(1, myTraj.N_traj, myTraj.traj_steps)
    plt.show()
