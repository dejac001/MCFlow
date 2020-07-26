import json
import os
import time

from mcflow import runAnalyzer


def output_json(path, run_type, data, save_old='Yes'):
    """Save data in json. If previous data exists, copy it to an old file

    :param path: path to *feed* directories
    :type path: str
    :param run_type: type of run tag
    :type run_type: str
    :param data: data to be saved
    :type data: dict
    :param save_old: whether or not to save old data, defaults to Yes
    :type save_old: str, optional
    """
    first_feed = list(data.keys())[0]
    for data_type in data[first_feed].keys():
        # convert data to json format
        data_to_save = {}
        for feed, vals in data.items():
            nextRun = runAnalyzer.findNextRun(os.path.join(path, feed, '1'), run_type)
            assert nextRun is not None, 'No run found'
            data_to_save[feed] = {
                '%s%i' % (run_type, nextRun - 1): vals[data_type].averages[feed],
                'time': time.time()
            }

        file_name = os.path.join(path, '%s-data.json' % data_type)
        if save_old == 'Yes' and os.path.isfile(file_name):
            old_file = os.path.join(path, 'old-%s-data.json' % data_type)
            if os.path.isfile(old_file):
                os.remove(old_file)
            os.rename(file_name, old_file)

        with open(file_name, 'w') as f:
            json.dump(data_to_save, f)


def outputGen_json(path, run_type, data, save_old='Yes'):
    """Save data in json. If previous data exists, copy it to an old file

    :param path: path to *feed* directories
    :type path: str
    :param run_type: type of run tag
    :type run_type: str
    :param data: data to be saved
    :type data: dict
    :param save_old: whether or not to save old data, defaults to Yes
    :type save_old: str, optional
    """
    data_type = 'general'
    data_to_save = {}
    for feed, val in data.items():
        nextRun = runAnalyzer.findNextRun('%s/%s/1/' % (path, feed), run_type)
        assert nextRun != None, 'No run found'
        run_key = '%s%i' % (run_type, nextRun - 1)
        data_to_save[feed] = {
            run_key: {},
            'time': time.time()
        }
        for key, val2 in val.items():
            data_to_save[feed][run_key][key] = val2

    file_name = path + '/%s-data.json'% data_type
    if save_old == 'Yes' and os.path.isfile(file_name):
        old_file = os.path.join(path, 'old-%s-data.json' % data_type)
        if os.path.isfile(old_file):
            os.remove(old_file)
        os.rename(file_name, path + '/old-%s-data.json' % data_type)

    with open(file_name, 'w') as f:
        json.dump(data_to_save, f)