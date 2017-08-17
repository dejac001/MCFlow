def write(data):
    f = open('scaled-data.dat','w')
    if len(data) == 2:
        for x,y in zip(data[0], data[1]):
            f.write('%e %e\n'%(x,y))
    elif len(data) == 4:
        for x,y, dx, dy in zip(data[0], data[1], data[2], data[3]):
            f.write('%e %e %e %e\n'%(x,y))
    f.close()

def convertUnits(my_data, xscale, yscale):
    data = np.matrix(my_data)
    data[0] = data[0]*xscale
    data[1] = data[1]*yscale
    if len(data) == 4:
        data[2] = data[2]*xscale
        data[3] = data[3]*yscale
    return data.tolist()


from plotting.read import readDat
import numpy as np

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Convert units')
    parser.add_argument('-f','--file',help='file name',type=str)
    parser.add_argument('-xs','--xscale',help='x scale factor',
                        type=float,default=1.0)
    parser.add_argument('-ys','--yscale',help='y scale factor',
                        type=float,default=1.0)
    args = vars(parser.parse_args())

    data = readDat(args['file'])
    scaled = convertUnits(data, args['xscale'], args['yscale'])
    write(scaled)
