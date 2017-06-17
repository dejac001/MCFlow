def make_grace(my_ax):
    my_ax.tick_params(axis='y',direction='in',which='both',labelsize=10,left=True,right=True)
    for dirn in ['top','bottom','left','right']:
        my_ax.spines[dirn].set_color('black')
        my_ax.spines[dirn].set_linewidth(1)

    my_ax.tick_params(colors='black',size=5.,width=1.)

def plot_transfer_acceptance(data, feed, run):
    # get list of move types and box transfers
    move_names = []
    boxes = []
    legend = []
    for key, value in data.items():
        # get list of boxes involved
        for boxes, results in value.items():
            if results['accepted']['mean'] > 0.0:
                legend.append( key + '-' + boxes)
    # setup dimensions of plot
    ind = np.arange(len(legend)) # the x locations for the groups
    width = 0.15   # the width of the bars

    fig, ax = plt.subplots()
    vals = []
    yerr = []
    for move in sorted(legend):
        molec, box = move.split('-')
        move_data = data[molec][box]
        percent_accepted = move_data['accepted']['mean'] / move_data['attempted']['mean'] * 100.
        stdev = percent_accepted*np.sqrt(
            math.pow(move_data['accepted']['stdev'] / move_data['accepted']['mean'],2)
            +
            math.pow(move_data['attempted']['stdev']/ move_data['attempted']['mean'],2)
        )
        vals.append(percent_accepted)
        yerr.append(stdev)
    ax.bar(ind, vals, width,yerr=yerr,color=np.random.rand(3,1))
    make_grace(ax)
    ax.set_ylabel('Percent Accepted')
    ax.set_xticks(ind + width)
    ax.set_xticklabels(sorted(legend),rotation=75.)
    ax.set_yscale('log')
    plt.subplots_adjust(left=0.11,right=0.97,top = 0.97, bottom=0.5)
    fig = plt.gcf()
    fig.set_size_inches(6.5, 6.5)
    fig.savefig('%s-%s.png'%(feed,run), dpi=300)
    plt.show()
import math, random
import numpy as np
import matplotlib.pyplot as plt
if __name__ == '__main__':
    import argparse, os, shelve
    from MCFlow.runAnalyzer import checkRun
    parser = argparse.ArgumentParser(description='plot swap and swatch accptances for given feed')
    parser.add_argument('-f','--feed',help='feed to analyze swaps and swatches for')
    parser.add_argument('-t','--type',help='type of run to analyze',default='equil-')
    args = vars(parser.parse_args())

    assert os.path.isfile('SWAP-data.db'), 'No SWAP data found'
    with shelve.open('SWAP-data.db') as db:
        assert args['feed'] in db.keys(), 'Feed not in database'
        data = {args['feed']: db[args['feed']] }
    my_run = checkRun(args['type'],[data],args['feed'])
    plot_transfer_acceptance(data[args['feed']][my_run],args['feed'],my_run)
