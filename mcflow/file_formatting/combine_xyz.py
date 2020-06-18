def combineXYZ(**kwargs):
    assert kwargs['newFile'], 'New file needed'
    assert len(kwargs['files']) > 1, 'Need more than one file if combining'
    all_xyz_data = {'atoms':[],'coords':[]}
    for file in kwargs['files']:
        data = reader.xyz(file)
        for key, value in data.items():
            all_xyz_data[key] += value

    writer.xyz(kwargs['newFile'], all_xyz_data)



if __name__ == '__main__':
    import argparse
    from MCFlow.file_formatting import reader, writer
    parser = argparse.ArgumentParser(description='combine .xyz files')
    parser.add_argument('-f','--files',help='list of xyz files to combine',
                        type=str,nargs='+')
    parser.add_argument('-n','--newFile',help='new file name',type=str,default='None')
    args = vars(parser.parse_args())

    if args['newFile'] == 'None':
        args['newFile'] = args['files'][0][:args['files'][0].rfind('/')] + '/combined.xyz'

    combineXYZ(**args)
