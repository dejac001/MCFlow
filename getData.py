def output_json(path, feeds, type, data):
    """Save data in json. If previous data exists, copy it to an old file"""
    for data_type in data.keys():
        # convert data to json format
        data_to_save = {}
        for feed in feeds:
            nextRun = runAnalyzer.findNextRun('%s/%s/1/'%(path, feed), type)
            assert nextRun is not None, 'No run found'
            data_to_save[feed] = {
                    '%s%i'%(type, nextRun-1) : data[data_type].averages[feed],
                    'time' : time.time()
            }

        file_name = path + '/%s-data.json'% data_type
        if os.path.isfile(file_name):
            os.rename(file_name, path + '/old-%s-data.json' % data_type)
        
        with open(file_name, 'w') as f:
            json.dump(data_to_save, f)


def outputGen_json(path, feeds, run_type, general_data):
    """Save data in json. If previous data exists, copy it to an old file"""
    data_type = 'general'
    data_to_save = {}
    for feed in feeds:
        nextRun = runAnalyzer.findNextRun('%s/%s/1/'%(path, feed), run_type)
        assert nextRun != None, 'No run found'
        run_key ='%s%i'%(run_type, nextRun-1) 
        data_to_save[feed] = {
                run_key : {},
                'time' : time.time()
        }
        for key, val in general_data[feed].items():
            data_to_save[feed][run_key][key] = val


    file_name = path + '/%s-data.json'% data_type
    if os.path.isfile(file_name):
        os.rename(file_name, path + '/old-%s-data.json' % data_type)

    with open(file_name, 'w') as f:
        json.dump(data_to_save, f)


import json
import os
import time

from MCFlow import runAnalyzer

if __name__ == '__main__':
    from parser import Results
    from runAnalyzer import getFileData, NoFilesAnalyzed

    my_parser = Results()
    my_parser.parser.add_argument('-l','--liq',help='whether or not liqid',type=bool,
                                  default=False)
    args = vars(my_parser.parse_args())

    args_to_send = args
    for key in ('rcut', 'nstep', 'time'): #
        if key in args_to_send.keys():
            args_to_send.pop(key)

    feeds = args['feeds']
    for feed in feeds:
        try:
            args['feeds'] = [feed]
            data, gen_data = getFileData(**args_to_send)
            output_json(args['path'], args['feeds'], args['type'], data)
            outputGen_json(args['path'], args['feeds'], args['type'], gen_data)
        except NoFilesAnalyzed:
            print('No files to analyzed for feed %s'%feed)
