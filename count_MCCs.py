if __name__ == '__main__':
    from analysis_parsers import Main
    from runAnalyzer import what2Analyze
    import numpy as np
    my_parser = Main()
    my_parser.other()
    args = vars(my_parser.parse_args())

    all_MCCS = []
    for feed in args['feeds']:
        feed_MCCS = []
        if args['verbosity'] > 0: print('-'*12 + 'Dir is %s'%feed + '-'*12)
        for sim in args['indep']:
            sim_MCCS = 0
            my_dir = '%s/%s/%i'%(args['path'],feed,sim)
            (old_begin, nfiles) = what2Analyze(my_dir, args['type'],
                                               args['guessStart'], args['interval'])
            for run in range(1, old_begin+nfiles):
                try:
                    with open('%s/%s%i/fort12.%s%i'%(my_dir,args['type'],run,args['type'],run)) as f:
                        file = f.readlines()
                        nbox = int(file[0].split()[2])
                        my_MCCS = int((len(file) - 1)/nbox)
                        sim_MCCS += my_MCCS
                except FileNotFoundError:
                    continue
            feed_MCCS.append(sim_MCCS)
        print(' ' *5 + 'MCCs was %5e +/- %5e for feed:%s'%( np.mean(feed_MCCS), np.std(feed_MCCS), feed ) )
        all_MCCS.append( np.mean(feed_MCCS) )
    print('MCCs for all feeds was %5e +/- %5e'%(np.mean(all_MCCS), np.std(all_MCCS)))
    print('MCCs for all feeds ranged from %5e to %5e'%(np.min(all_MCCS), np.max(all_MCCS)))
