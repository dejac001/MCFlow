import json
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

def clean_U(U_file_name, N_file_name):
    with open(U_file_name, 'r') as f:
        data = json.load(f)
    with open(N_file_name, 'r') as g:
        N_data = json.load(g)

    clean_data = {}
    for feed, val in data.items():
        try:
            N_box1_mean = N_data[feed]['prod-2'][moltype]['box1']['mean']
            N_box2_mean = N_data[feed]['prod-2'][moltype]['box2']['mean']
            N_box1_std = N_data[feed]['prod-2'][moltype]['box1']['stdev']
            N_box2_std = N_data[feed]['prod-2'][moltype]['box2']['stdev']
            box1_mean = val['prod-2']['box1']['mean']
            box2_mean = val['prod-2']['box2']['stdev']
            box1_std  = val['prod-2']['box1']['mean']
            box2_std  = val['prod-2']['box2']['stdev']
            if N_box1_mean > 0. and N_box2_mean > 0.:
                clean_data[feed] = {
                                    'mean': (box1_mean/N_box1_mean - box2_mean/N_box2_mean)*R,
                                    'stdev': math.sqrt(box1_std*box1_std/N_box1_mean/N_box1_mean + box2_std*box2_std/N_box2_mean/N_box2_mean
    +N_box1_std*N_box1_std*(box1_mean/N_box1_mean/N_box1_mean)**2 +
    N_box2_std*N_box2_std*(box2_mean/N_box2_mean/N_box2_mean)**2 
    )*R
            }
        except KeyError:
            continue

    with open('dU.json', 'w') as f:
        json.dump(clean_data, f, indent=4)

    with open('dU.csv', 'w') as f:
        f.write('"Temperature [K]","Vapor Phase Boxlength [Angstrom]","dU mean [kJ/mol]","dU standard deviation [kJ/mol]"\n')
        for key, val in clean_data.items():
            f.write('%s,%s,%e,%e\n'%(key.split('/')[0][:-1], key.split('/')[1], val['mean'], val['stdev']))


moltype='1'


if __name__ == '__main__':
    clean_dG('dG-data.json')
    clean_U('U-data.json', 'N-data.json')
