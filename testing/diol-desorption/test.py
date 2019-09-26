import getData


def default_args():
    args = getData.my_parser()
    args['type'] = 'prod-'
    args['feeds'] = ['E-MFI-343K-256-1']
    args['indep'] = [1]
    args['verbosity'] = 1
    args['save_old_data'] = 'No'
    return args


def test_production_analysis():
    getData.main(default_args())


def test_production_analysis_with_energies():
    args = default_args()
    args['energies'] = 'Yes'
    getData.main(args)


if __name__ == '__main__':
    test_production_analysis_with_energies()
