class Trajectory():
    def __init__(self, boxList, feeds, indepRange):
        self.N_traj = { '%s'%f: { B :{ k: []
                                       for k in indepRange}
                                  for B in boxList }
                         for f in feeds}
        self.traj_steps =  { '%s'%f: { k: [] for k in indepRange}
                             for f in feeds}
    def addSim(self, N, feed, sim, tag, old_begin, molNum, boxList, numSkip):
        offset = int(next(open('%s%i/config.%s%i'%(tag, old_begin, tag, old_begin))))
        index = 0
        begin_index = 0
        for index in range(len(N[molNum][boxList[0]])):
            if index % numSkip == 0:
                for box in boxList:
                    self.N_traj['%s'%feed][box][sim].append( np.mean(N[molNum][box][begin_index: (index + 1)]) )
                self.traj_steps['%s'%feed][sim].append(index+offset)
                begin_index = index

from plotter import *

if __name__ == '__main__':
    from file_formatting import reader
    import runAnalyzer
    import matplotlib.pyplot as plt
    from parser import Results

    my_parser = Results()
    my_parser.parser.description = 'Plot trajectory from inputting molecules and boxes'
    my_parser.parser.add_argument('-sk','--numSkip',help ='number of MCCs to skip', type = int, default = 100)

    args = vars(my_parser.parse_args())

    boxes = ['box%s'%args['box']]

    for feed in args['feeds']:
        print('starting trajectory for feed {}.....'.format(feed))
        for sim in args['indep']:
            iPath = '%s/%s/%i'%(args['path'], feed,  sim)
            (old_begin, nfiles) = runAnalyzer.what2Analyze(iPath, args['type'], args['guessStart'], args['interval'])
            print('old_begin = {}, nfiles = {}'.format(old_begin, nfiles))
            N, P, boxLengths, ncycle_old, molWeights, U = reader.read_fort12(iPath, old_begin, nfiles, tag=args['type'])
            if feed == args['feeds'][0] and sim == args['indep'][0]:
                myTraj = Trajectory(boxes, args['feeds'], args['indep'])
            myTraj.addSim(N, feed, sim, args['type'], old_begin, args['mol'], boxes, args['numSkip'])
    suptitle='%s, %s-(%i-%i); mol#%s'% (args['path'][(args['path'].rfind('/')+1):],
                                        args['type'] ,old_begin, old_begin+nfiles-1,args['mol'])
    plotTrajectory(1, suptitle, myTraj.N_traj, myTraj.traj_steps)
    plt.show()
