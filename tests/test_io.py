import os
import numpy as np
path = os.path.abspath(os.path.dirname(__file__))
from scripts.change_pressure import set_pressure


def test_fort4():
    from mcflow.file_formatting.reader import read_fort4
    from mcflow.file_formatting.writer import write_fort4
    data = read_fort4(os.path.join(path, 'test-data', 'fort.4'))
    assert int(data['&mc_shared']['seed']) == 1, 'Incorrect seed'
    assert int(data['&mc_shared']['nbox']) == 2, 'Incorrect box'
    assert int(data['&mc_shared']['nchain']) == 600, 'Incorrect nchain'
    assert np.isclose(float(data['SIMULATION_BOX']['box1']['pressure']), 0.001), 'Incorrect pressure'
    data = set_pressure(data, 0.01)
    write_fort4(data, os.path.join('test-data', 'fort.4.new'))
    new_file_name = os.path.join(path, 'test-data', 'fort.4.new')
    new_data = read_fort4(new_file_name)
    assert np.isclose(float(new_data['SIMULATION_BOX']['box1']['pressure']), 0.01), 'Incorrect pressure'
    assert int(new_data['&mc_shared']['seed']) == 1, 'Incorrect seed'
    assert int(new_data['&mc_shared']['nbox']) == 2, 'Incorrect box'
    assert int(new_data['&mc_shared']['nchain']) == 600, 'Incorrect nchain'
    os.remove(new_file_name)



if __name__ == '__main__':
    test_fort4()