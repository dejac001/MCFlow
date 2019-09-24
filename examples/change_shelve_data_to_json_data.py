import json
import shelve
import math
import numpy as np

def convert_data(data):
    if isinstance(data, dict):
        for key in data.keys():
            convert_data(data[key])
    elif isinstance(data, np.ndarray):
        data = data.tolist()


def convert(file_name):
    data = {}
    with shelve.open(file_name) as d:
        for key, value in d.items():
            data[key] = value

    convert_data(data)

    with open(file_name[:-3] + '.json', 'w') as f:
        json.dump(data, f)


if __name__ == '__main__':
    convert('dG-data.db')
    convert('U-data.db')
