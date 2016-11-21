def outputDB(path, feeds, type, data):
    import time
    for data_type in data.keys():
        with shelve.open(path + '/%s-data.db'%data_type,writeback=True) as db:
            for feed in feeds:
                nextRun = runAnalyzer.findNextRun('%s/%s/1/'%(path, feed), type)
                if (feed in db.keys()) and (type in db[feed].keys()):
                    # some data has already been input into db
                    if (('%s%i'%(type,nextRun-1) not in db[feed].keys())
                        or (db[feed]['time'] < os.stat('%s/%s/1/%s%i'%(path, feed,type,
                                                                       nextRun-1)))):
                        # replace old data
                        for run in db[feed].keys():
                            if type in run:
                                db[feed].pop(run)
                        db[feed]['%s%i'%(type, nextRun-1)] = data[data_type].averages[feed]
                        db[feed]['time'] = time.time()
                else:
                    db[feed] = {}
                    db[feed]['%s%i'%(type, nextRun-1)] = data[data_type].averages[feed]
                    db[feed]['time'] = time.time()

def outputGenDB(path, feeds, type, general_data):
    import time
    with shelve.open(path + '/general-data.db', writeback=True) as db:
        for feed in feeds:
            nextRun = runAnalyzer.findNextRun('%s/%s/1/'%(path, feed), type)
            if (feed in db.keys()) and (type in db[feed].keys()):
                # some data has already been input into db
                if (('%s%i'%(type,nextRun-1) not in db[feed].keys())
                        or (db[feed]['time'] < os.stat('%s/%s/1/%s%i'%(path, feed,type,
                                                                       nextRun-1)))):
                        # replace old run
                        for run in db[feed].keys():
                            if type in run:
                                db[feed].pop(run)
                        db[feed]['%s%i'%(type, nextRun-1)] = {}
                        for key in general_data[feed].keys():
                            db[feed]['%s%i'%(type, nextRun-1)][key] = general_data[feed][key]
                        db[feed]['time'] = time.time()
            else:
                db[feed] = {}
                db[feed]['%s%i'%(type, nextRun-1)] = {}
                for key in general_data[feed].keys():
                    db[feed]['%s%i'%(type, nextRun-1)][key] = general_data[feed][key]
                db[feed]['time'] = time.time()



import os, shelve
from MCFlow import runAnalyzer

if __name__ == '__main__':
    from parser import Results
    from runAnalyzer import getFileData, NoFilesAnalyzed

    my_parser = Results()
    my_parser.parser.add_argument('-F','--force',help='whether or not to force update of '
                                                      'feed info', type=bool,
                                  default=False)
    my_parser.parser.add_argument('-l','--liq',help='whether or not liqid',type=bool,
                                  default=False)
    args = vars(my_parser.parse_args())

    old_feeds = []
    # determine which feeds we actually need to analyze
    if not args['force']:
        for testDB in [i for i in os.listdir(args['path']) if '.db' in i]: # only need to do once
            with shelve.open(testDB) as db:
                for feed in [i for i in db.keys() if i in args['feeds']]:
                    nextRun = runAnalyzer.findNextRun('%s/%s/1/'%(args['path'], feed), args['type'])
                    if (('%s%i'%(args['type'], nextRun-1) in db[feed].keys()) and
                            (db[feed]['time'] >= os.stat('%s/%s/1/%s%i'%
                                                             (args['path'], feed, args['type'], nextRun-1)).st_mtime)):
                        # we have already analyzed this run
                        args['feeds'].remove(feed)
            break # no need to go through loop second time
        assert args['feeds'], 'No feeds to analyze'

    args_to_send = args
    for key in ('rcut', 'nstep', 'time', 'force', 'box'):
        if key in args_to_send.keys():
            args_to_send.pop(key)

    feeds = args['feeds']
    for feed in feeds:
        try:
            args['feeds'] = [feed]
            data, gen_data = getFileData(**args_to_send)
            outputDB(args['path'], args['feeds'], args['type'], data)
            outputGenDB(args['path'], args['feeds'], args['type'], gen_data)
        except NoFilesAnalyzed:
            print('No files to analyzed for feed %s'%feed)
