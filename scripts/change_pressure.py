"""
Change Pressure
===============
Change the pressure of a set of input files to a new set of input files in a different directory

"""
import os
import sys

# my dir
my_dir = os.path.dirname(os.path.abspath(__file__))

# add parent directory to pythonpath
sys.path.append(os.path.join(my_dir, '..'))
import shutil
from mcflow.file_formatting.reader import read_fort4
from mcflow.file_formatting.writer import write_fort4


def set_pressure(data, val):
    for box in data['SIMULATION_BOX'].keys():
        data['SIMULATION_BOX'][box]['pressure'] = '%3.2e' % val
    return data


def change_pressure_dirs(old_path, new_path, P_new_MPa):
    """

    :param old_path: name of old working directory including all input files
    :type old_path: str
    :param new_path: name of new working directory
    :type new_path: str
    :param P_new_MPa: new pressure in MPa
    :type P_new_MPa: float
    """
    assert os.path.isdir(old_path), 'Old path {} not found'.format(old_path)
    assert old_path != new_path, 'Both paths are the same!'

    old = {
        'restart': os.path.join(old_path, 'fort.77'),
        'input': os.path.join(old_path, 'fort.4'),
        'struc': os.path.join(old_path, 'input_struc.xyz'),
        'zeo struc': os.path.join(old_path, 'zeolite.cif'),
        'zeo tab pot': os.path.join(old_path, 'zeolite.ztb'),
        'parameter': os.path.join(old_path, 'topmon.inp'),
    }
    for key in ('input', 'restart', 'parameter'):
        assert os.path.isfile(old[key]), 'Old %s file not found!' % key

    if not os.path.isdir(new_path):
        os.makedirs(new_path)

    new = {
        'restart': os.path.join(new_path, 'fort.77'),
        'input': os.path.join(new_path, 'fort.4'),
        'struc': os.path.join(new_path, 'input_struc.xyz'),
        'zeo struc': os.path.join(new_path, 'zeolite.cif'),
        'zeo tab pot': os.path.join(new_path, 'zeolite.ztb'),
        'parameter': os.path.join(new_path, 'topmon.inp'),
    }

    # copy input files
    # restart
    shutil.copy(old['restart'], new['restart'])

    for f_name in ('struc', 'zeo struc', 'zeo tab pot', 'parameter'):
        if os.path.islink(old[f_name]):
            source = os.readlink(old[f_name])
            try:
                os.symlink(source, new[f_name])
            except FileExistsError:
                os.remove(new[f_name])
                os.symlink(source, new[f_name])
        elif os.path.isfile(old[f_name]):
            try:
                os.symlink(old[f_name], new[f_name])
            except FileExistsError:
                os.remove(new[f_name])
                os.symlink(old[f_name], new[f_name])

    # change pressure
    data = read_fort4(old['input'])
    data = set_pressure(data, P_new_MPa)
    write_fort4(data, new['input'])


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-od', '--olddir', help='old directory of files')
    parser.add_argument('-nd', '--newdir', help='new directory of files')
    parser.add_argument('-P', '--pressure', help='New pressure [MPa]', type=float)
    args = vars(parser.parse_args())
    change_pressure_dirs(args['olddir'], args['newdir'], args['pressure'])
