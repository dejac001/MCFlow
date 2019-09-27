import json
import numpy as np
import math
R = 8.314/1000.


def clean_dG(file_name):
    with open(file_name, 'r') as f:
        data = json.load(f)

    clean_data = {}
    for feed, val in data.items():
        try:
            clean_data[feed] = {
                                'mean': val['prod-2'][moltype]['box2--box1']['mean'],
                                'stdev': val['prod-2'][moltype]['box2--box1']['stdev']
            }
        except KeyError:
            continue

    with open('dG.json', 'w') as f:
        json.dump(clean_data, f, indent=4)

    with open('dG.csv', 'w') as f:
        f.write('"Temperature [K]","Vapor Phase Boxlength [Angstrom]","dG mean [kJ/mol]","dG standard deviation [kJ/mol]"\n')
        for key, val in clean_data.items():
            f.write('%s,%s,%e,%e\n'%(key.split('/')[0][:-1], key.split('/')[1], val['mean'], val['stdev']))

def clean_dU(file_name):
    with open(file_name, 'r') as f:
        data = json.load(f)

    clean_data = {}
    for feed, val in data.items():
        if np.isnan(val['prod-2']['box2-->box1']['mean']):
            continue
        try:
            clean_data[feed] = {
                                'mean': val['prod-2']['box2-->box1']['mean'],
                                'stdev': val['prod-2']['box2-->box1']['stdev']
            }
        except KeyError:
            continue

    with open('dU.json', 'w') as f:
        json.dump(clean_data, f, indent=4)

    with open('dU.csv', 'w') as f:
        f.write('"Temperature [K]","Vapor Phase Boxlength [Angstrom]","dU mean [kJ/mol]","dU standard deviation [kJ/mol]"\n')
        for key, val in clean_data.items():
            f.write('%s,%s,%e,%e\n'%(key.split('/')[0][:-1], key.split('/')[1], val['mean'], val['stdev']))


moltype='8'


if __name__ == '__main__':
    clean_dG('dG-data.json')
    clean_dU('dU-data.json')
