def checkYoSelf(path, feeds, nindep, count_start, tag=equilName):

    for feed in feeds:
        error = False
        for i in nindep:
            os.chdir('%s/%s/%i'%(path, feed, i ))
            run = True
            count = count_start
            while run:
                for fileName in ['run.','fort12.']:
                    try:
                        with open(tag + str(count+1) + '/' +fileName + tag + str(count+1),'r') as f:
                            file_str = f.read()
                            if '\x00' in file_str:
                                print(os.getcwd(),tag,count+1,fileName,
                                      'is a binary file.You need to ignore this for analysis')
                                error = True
                                break
                            if fileName == 'run.':
                                if 'Program ended' in file_str:
                                    count += 1
                                else:
                                    print(os.getcwd(),'has not finished',tag , count+1)
                                    error = True
                                    break
                    except FileNotFoundError:
                        run = False
                        break
            if i == 1:
                count_prev = count
            elif count < count_prev:
                print('%s/%i'%(feed,i), 'is only on ',tag ,count
                      ,'while all others are on ',tag ,count_prev)
            elif count > count_prev:
                print(feed + '/1', 'is only on ', tag , count_prev,
                      'while another one is on', tag ,count)

        if not error:
            print('all runs for ',feed,'finished up through', tag, count)

import os
from MCFlow.file_organization import equilName, prodName

if __name__ == '__main__':

    from analysis_parsers import Main
    my_parser = Main()
    my_parser.other()
    my_parser.parser.description = "Determine how many runs have been successfully finished"
    args = vars(my_parser.parse_args())

    if args['guessStart'] == 1:
        # make default == 0 for checking purposes
        args['guessStart'] = 0


    checkYoSelf(args['path'], args['feeds'], args['indep'], args['guessStart'], tag=args['type'])

