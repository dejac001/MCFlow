import os
import shutil
import sys

# my dir
from mcflow.calculate_probs import NoSwapsAccepted, analyzeTransfers

my_dir = os.path.dirname(os.path.abspath(__file__))

# add parent directory to pythonpath
sys.path.append(os.path.join(my_dir, '..'))
from mcflow.save_data import output_json, outputGen_json
from mcflow.runAnalyzer import getFileData, findNextRun
from mcflow.file_formatting import writer, reader


def main(feeds, verbosity=0, type='equil-', path='', indep=None, interval=0, guessStart=1, time=345000,
         rcut=0.5, nstep=50000, liq=False, mol=None, energies=None, box=None, debug=True):
    if indep is None:
        indep = range(1, 9)
    data = {}
    gen_data = {}

    for feed in feeds:
        data[feed], gen_data[feed] = getFileData(feeds=[feed], indep=indep, path=path, type=type, guessStart=guessStart,
                                                 interval=interval, verbosity=verbosity, liq=liq, mol=mol,
                                                 energies=energies, box=box, debug=debug)
        nbox = len(data[feed]['rho'].averages[feed].keys())
        # TODO: format return from analyzeTransfers to fit well with fort4 data dict
        try:
            (newSwaps, newSwatches, pctAct,
             nActCycle, normSwaps, normSwatches,
             pmvol, pswatch_norm, pswap_norm) = analyzeTransfers(data[feed]['SWAP'].averages[feed],
                                                                 gen_data[feed][feed]['ncycle'],
                                                                 # nActPerCycle=3656/200, tavol=0.3,
                                                                 gen_data[feed][feed]['compositions'])
        except NoSwapsAccepted:
            print('no swaps accepted for {}'.format(feed))
            continue
        print(newSwaps)
        # read and write
        for sim in indep:
            my_path = os.path.join(path, feed, '%i' % sim)
            try:
                input_data = reader.read_fort4(os.path.join(my_path, 'fort.4'))
            except FileNotFoundError:
                input_data = reader.read_fort4(os.path.join(my_path, 'old-fort4'))
            # add in swaps
            if input_data['&mc_swap']['pmswap'][-2:] == 'd0':
                input_data['&mc_swap']['pmswap'] = input_data['&mc_swap']['pmswap'][:-2]
            pswap_old = float(input_data['&mc_swap']['pmswap'])
            input_data['&mc_swap']['pmswap'] = '%e' % pswap_norm
            input_data['&mc_shared']['seed'] = '%i' % sim
            nmolty = int(input_data['&mc_shared']['nmolty'])
            assert len(newSwaps) == nmolty, 'Incorrect swaps: too many'
            total_pswap = 0.
            for molNum in range(1, nmolty + 1):
                for mol in newSwaps.keys():
                    if mol.split()[0] == '%i' % molNum:
                        key = mol
                        val = newSwaps[key]
                mol = 'mol%s' % (key.split()[0])
                # total value for molecule type
                total = val['total']
                if (total < 1e-12):
                    my_tot = '-1.0'
                else:
                    total_pswap = total_pswap + total
                    my_tot = '%4e' % total_pswap
                input_data['&mc_swap']['pmswmt'][mol] = my_tot
                # add in swap directions
                dirs = [i for i in val.keys() if i != 'total']
                input_data['MC_SWAP'][mol]['nswapb'] = len(dirs)
                input_data['MC_SWAP'][mol]['pmswapb'] = []
                input_data['MC_SWAP'][mol]['box1 box2'] = []
                my_sum = 0.
                for drn in dirs:
                    my_sum = my_sum + val[drn]  # numpy object, can't add to itself
                    input_data['MC_SWAP'][mol]['pmswapb'].append(my_sum)
                    input_data['MC_SWAP'][mol]['box1 box2'].append(list(map(int, drn)))
            # take out extraneous swatches
            unique_swatches = []
            for iswatch, info in input_data['MC_SWATCH'].items():
                in_line_1 = []
                in_line_2 = []
                splist = []
                for i, line in enumerate(info.split('\n')):
                    if not line.startswith('!'):
                        if len(in_line_1) == 0:
                            in_line_1 = list(line.split())
                        elif len(in_line_2) == 0:
                            in_line_2 = list(line.split())
                        else:
                            splist.append(list(map(int, line.split())))
                            if len(splist) == int(in_line_1[2]): break
                swatch_label = ('! moltyp1<->moltyp2 nsampos 2xncut\n' + ' '.join(in_line_1) + '\n'
                                + '! gswatc 2x(ifrom, iprev)\n' + ' '.join(in_line_2) + '\n' +
                                '! splist\n')
                for sameBead in splist:
                    swatch_label += '%i %i\n' % (sameBead[0], sameBead[1])
                if swatch_label not in unique_swatches:
                    unique_swatches.append(swatch_label)
            # add in swatches
            input_data['MC_SWATCH'] = {}
            input_data['&mc_swatch'] = {'pmsatc': {}, 'nswaty': '0', 'pmswat': '0'}
            if pswatch_norm - pmvol < 1e-8:
                input_data['&mc_swatch']['pmswat'] = '-1.0'
            else:
                input_data['&mc_swatch']['pmswat'] = '%e' % pswatch_norm
            nSwatch = 0
            pmswatch_total = 0.
            for key, value in newSwatches.items():
                mol_names = key.split(' and ')
                mol1_num, mol2_num = [i.split()[0] for i in mol_names]
                for swatch in unique_swatches:
                    mol1, mol2 = swatch.split('\n')[1].split()[:2]
                    if mol1 == mol1_num and mol2 == mol2_num:
                        nSwatch += 1
                        s_num = 'swatch%i' % nSwatch
                        info = swatch + '! nswtcb pmswtcb\n'
                        # add in direstions
                        dirs = [i for i in value.keys() if i != 'total']
                        info += '%i' % len(dirs)
                        my_sum = 0.
                        for drn in sorted(dirs):
                            my_sum += value[drn]
                            info += ' %e' % (my_sum)
                        info += '\n! box numbers\n'
                        for drn in sorted(dirs):
                            info += drn[0] + ' ' + drn[1] + '\n'
                        input_data['MC_SWATCH'][s_num] = info
                        pmswatch_total += value['total']
                        input_data['&mc_swatch']['pmsatc'][s_num] = '%e' % pmswatch_total
            input_data['&mc_swatch']['nswaty'] = '%i' % nSwatch
            # add in other stuff
            #   general info
            input_data['&mc_shared']['time_limit'] = '%i' % time
            input_data['&mc_shared']['nstep'] = '%i' % nstep
            input_data['&mc_shared']['iratio'] = '500'
            input_data['&mc_shared']['rmin'] = '1.0'
            iblock = int(nstep / 10)
            if iblock > 1000: iblock = 1000
            if '&analysis' not in input_data.keys():
                input_data['&analysis'] = {}
            input_data['&analysis']['iblock'] = '%i' % iblock
            input_data['&mc_volume']['pmvol'] = '%e' % pmvol
            #   other probabilities
            if input_data['&mc_cbmc']['pmcb'][-2:] == 'd0':
                input_data['&mc_cbmc']['pmcb'] = input_data['&mc_cbmc']['pmcb'][:-2]
            pmcb_old = float(input_data['&mc_cbmc']['pmcb']) - pswap_old
            pswap_new = pswap_norm
            if pmcb_old <= 0.:
                print('pmcb was zero previously, keeping it that way')
                pmcb_old = 0.
                pmtra_old = (1 - pswap_old) / 2.
                pmcb = -1.0
                pmtra = (1 - pswap_new) * pmtra_old / (1 - pswap_old) + pswap_new
            else:
                pmcb = (1 - pswap_new) * pmcb_old / (1 - pswap_old) + pswap_new
                pmtra_old = (1 - pswap_old - pmcb_old) / 2
                pmtra = (1 - pswap_new) * pmtra_old / (1 - pswap_old) + pmcb
            assert pmtra < 1., 'Pmtra %e too large!' % pmtra
            input_data['&mc_cbmc']['pmcb'] = '%e' % pmcb
            input_data['&mc_simple']['pmtra'] = '%e' % pmtra
            #   box lengths
            for box in range(1, int(input_data['&mc_shared']['nbox']) + 1):
                if 'boxlx' not in data[feed].keys():
                    continue
                boxLengths = data[feed]['boxlx'].averages[feed]['box%i' % box]['mean']
                input_data['SIMULATION_BOX']['box%i' % box]['dimensions'] = (
                        '%8f %8f %8f' % (boxLengths, boxLengths, boxLengths)
                )
                if (boxLengths > 100.) and (rcut < 1.):
                    input_data['SIMULATION_BOX']['box%i' % box]['rcut'] = '%5e' % (boxLengths * rcut)
                if (input_data['SIMULATION_BOX']['box%i' % box]['defaults'].split()[-2] == 'T'):
                    # if ideal gas
                    input_data['SIMULATION_BOX']['box%i' % box]['rcut'] = '14.0'
            # make new run
            old_fort4 = [i for i in os.listdir(my_path) if 'fort.4.%s' % type in i]
            for fort4 in old_fort4:
                os.remove(os.path.join(my_path, fort4))
            try:
                shutil.move(os.path.join(my_path, 'fort.4'), os.path.join(my_path, 'old-fort4'))
            except FileNotFoundError:
                pass
            nextRun = 'fort.4.%s%i' % (type, findNextRun(my_path, type))
            writer.write_fort4(input_data, os.path.join(my_path, nextRun))
            input_data = None
    output_json(path, type, data)
    outputGen_json(path, type, gen_data)


if __name__ == '__main__':
    from analysis_parsers import Change

    args = vars(Change().parse_args())
    feeds = args.pop('feeds')
    main(feeds, **args)
