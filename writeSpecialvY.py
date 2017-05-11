
from MCFlow.writeXvY import *

if __name__ == '__main__':
    # TODO: find way to have conditional arguments based off of what xaxis and yaxis are
    # (in other words, we have an excess of variables as arguments)
    my_parser = Plot()
    my_parser.axes()
    my_parser.isotherm()
    my_parser.parser.add_argument('-N-data','--NDATA',help='different file to replace N-data.db with',
                                  type=str)

    args = vars(my_parser.parse_args())
    assert args['yaxis'], 'No y axis chosen for plot'
    assert args['xaxis'], 'No x axis chosen for plot'
    assert args['NDATA'], 'Specialed N-data.db needed to be provided'


    if args['xaxis'] == 'C':
        my_plotter = LiqAds(**args)
    elif args['xaxis'] == 'Pig':
        # TODO: add alternative way to calculate P w/ P-data.db & using mole fraction in box
        my_plotter = IdealGasAds(**args)
    elif args['xaxis'] == 'rho':
        my_plotter = RhoBoxAds(**args)
    elif args['xaxis'] == 'Pbox':
        my_plotter = GasBoxAds(**args)
    elif args['xaxis'] == 'Pi':
        my_plotter = GasMolAds(**args)
    elif args['xaxis'] == 'Q':
        my_plotter = LoadAds(**args)
    elif args['xaxis'] == 'x':
        my_plotter = MoleFrac(**args)
    old_index = my_plotter.files.index('N-data.db')
    my_plotter.files[old_index] = args['NDATA']
    my_plotter.readDBs()

    for feed in args['feeds']:
        # determine if run has been completed
        my_plotter.run = checkRun(args['type'], my_plotter.variables, feed)
        print(feed, my_plotter.run)
        my_plotter.feed = feed
        if args['yaxis'] == 'Q':
            my_plotter.QvX()
        elif args['yaxis'] == 'S':
            my_plotter.SvX()
        elif args['yaxis'] == 'dG':
            my_plotter.dGvX()
        elif args['yaxis'] == 'R':
            my_plotter.RvX()
        elif args['yaxis'] == 'X':
            my_plotter.XvX()
        elif args['yaxis'] == 'dH':
            my_plotter.dHvX()
        elif args['yaxis'] == 'dU':
            my_plotter.dUvX()
        elif args['yaxis'] == 'Pigxy':
            my_plotter.Pig_xy()
        elif args['yaxis'] == 'Txy':
            my_plotter.Txy()
