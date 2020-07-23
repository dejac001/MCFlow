from mcflow.save_data import output_json, outputGen_json
from mcflow import runAnalyzer


def main(args):
    args_to_send = args
    for key in ('rcut', 'nstep', 'time'):  #
        if key in args_to_send.keys():
            args_to_send.pop(key)

    feeds = args['feeds']
    data = {}
    gen_data = {}
    for feed in feeds:
        try:
            args['feeds'] = [feed]
            data[feed], gen_data[feed] = runAnalyzer.getFileData(**args_to_send)
        except runAnalyzer.NoFilesAnalyzed:
            print('No files to analyzed for feed %s' % feed)

    output_json(args['path'], args['type'], data, args['save_old_data'])
    outputGen_json(args['path'], args['type'], gen_data, args['save_old_data'])


def my_parser():
    from analysis_parsers import Results
    my_parser = Results()
    my_parser.parser.add_argument('-l', '--liq', help='whether or not liqid',
                                  type=bool,
                                  default=False)
    my_parser.parser.add_argument('-e', '--energies', help='whether or not to calculate dU,'
                                                           ' dH of transfer (takes extra time)',
                                  type=str, choices=['Yes', 'No'], default='No')
    my_parser.parser.add_argument('-ss', '--save_old_data', help='whether or not to save old data',
                                  type=str, choices=['Yes', 'No'], default='Yes')
    return vars(my_parser.parse_args())



if __name__ == '__main__':
    arguments = my_parser()
    main(arguments)
