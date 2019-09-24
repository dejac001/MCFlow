import json
import math


def clean_dG(file_name):
    with open(file_name, 'r') as f:
        data = json.load(f)

    clean_data = {}
    for feed, val in data.items():
        try:
            clean_data[feed] = {
                                'mean': val['prod-2']['1']['box2--box1']['mean'],
                                'stdev': val['prod-2']['1']['box2--box1']['stdev']
            }
        except KeyError:
            continue

    with open('dG.json', 'w') as f:
        json.dump(clean_data, f, indent=4)

    with open('dG.csv', 'w') as f:
        f.write('"Temperature [K]","Vapor Phase Boxlength [Angstrom]","dG mean [kJ/mol]","dG standard deviation [kJ/mol]"\n')
        for key, val in clean_data.items():
            f.write('%s,%s,%e,%e\n'%(key.split('/')[0][:-1], key.split('/')[1], val['mean'], val['stdev']))

def clean_U(file_name):
    with open(file_name, 'r') as f:
        data = json.load(f)

    clean_data = {}
    for feed, val in data.items():
        try:
            box1_mean = val['prod-2']['box1']['mean']
            box2_mean = val['prod-2']['box2']['stdev']
            box1_std  = val['prod-2']['box1']['mean']
            box2_std  = val['prod-2']['box2']['stdev']
            clean_data[feed] = {
                                'mean': (box1_mean - box2_mean)*8.314,
                                'stdev': math.sqrt(box1_std*box1_std + box2_std*box2_std)*8.314
            }
        except KeyError:
            continue

    with open('dU.json', 'w') as f:
        json.dump(clean_data, f, indent=4)

    with open('dU.csv', 'w') as f:
        f.write('"Temperature [K]","Vapor Phase Boxlength [Angstrom]","dU mean [kJ/mol]","dU standard deviation [kJ/mol]"\n')
        for key, val in clean_data.items():
            f.write('%s,%s,%e,%e\n'%(key.split('/')[0][:-1], key.split('/')[1], val['mean'], val['stdev']))

if __name__ == '__main__':
    clean_dG('dG-data.json')
    clean_U('U-data.json')
