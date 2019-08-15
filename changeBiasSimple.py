def newBias(N, biasOld, T, impMols, impBox):
    bias_new = {mol: {box: biasOld[mol][box] for box in biasOld[mol].keys()} for mol in biasOld.keys()}
    for mol in impMols:
        n_tot_imp = sum(N[mol][i]['mean'] for i in N[mol].keys() if 'box' in i)
        n_each_box = n_tot_imp / len(impBox)
        print(mol, n_tot_imp, n_each_box)
        for box in impBox:
            bias_add = T*math.log(N[mol][box]['mean']/n_each_box)
            bias_add_stdev = T*N[mol][box]['stdev']/N[mol][box]['mean']
            if abs(bias_add_stdev) < abs(bias_add):
                print('adding %f to %s in %s'%(bias_add, mol, box))
                bias_new['mol' + mol][box] = bias_new['mol' + mol][box] + bias_add
    return bias_new

import math
from MCFlow.chem_constants import N_av
import sys
sys.path.append('MCFlow')

if __name__ == '__main__':
    from MCFlow.parser import Change
    from MCFlow.runAnalyzer import getFileData, findNextRun
    from MCFlow.file_formatting import writer, reader
    from MCFlow.getData import outputGenDB, outputDB

    my_parser = Change()
    my_parser.parser.add_argument('-I','--impurityMols',help='impurity molecules',type=str,nargs='+')
    my_parser.parser.add_argument('-ib','--impurityBoxes',help='impurity boxes',type=str, nargs='+')
    
    args = vars(my_parser.parser.parse_args())
    feeds = args.pop('feeds')

    for feed in feeds:
        args['feeds'] = [feed]
        data, gen_data = getFileData(**args)
        nbox = len(data['rho'].averages[feed].keys())
        try:
            input_data = reader.read_fort4(args['path']+'/' + feed + '/1/fort.4')
        except FileNotFoundError:
            input_data = reader.read_fort4(args['path']+'/' + feed + '/1/old-fort4')
        bias =  newBias(data['N'].averages[feed], input_data['UNIFORM_BIASING_POTENTIALS'],
                            gen_data[feed]['temperature'], args['impurityMols'], args['impurityBoxes'])
        # change data
        input_data['UNIFORM_BIASING_POTENTIALS'] = bias
        input_data['&mc_shared']['nstep'] = '%i'%args['nstep']
#       iprint, iblock = iaverage(args['nstep'])
#       input_data['&analysis']['iprint'] = '%i'%iprint
#       input_data['&analysis']['iblock'] = '%i'%iblock

        # write new files
        for sim in gen_data[feed]['indepSims']:
            input_data['&mc_shared']['seed'] = '%i'%sim
            my_path = '%s/%s/%i'%(args['path'],feed,sim)
            nextRun = findNextRun(my_path, args['type'])
            new_file = '%s/fort.4.%s%i'%(my_path, args['type'], nextRun)
            writer.write_fort4(input_data, new_file)
        outputDB(args['path'], args['feeds'],args['type'], data )
        outputGenDB(args['path'], args['feeds'],args['type'], gen_data )
