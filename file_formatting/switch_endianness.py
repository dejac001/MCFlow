#!/usr/bin/env python

import sys
import struct

def swap(fin, fou, nbytes=8):
    bytes = fin.read(nbytes)
    fou.write(bytes[::-1])
    return bytes

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "Usage: %s /path/to/old /path/to/new\n" % sys.argv[0]
        quit()

    if sys.byteorder == 'little':  # assume file uses different endianness than the host machine
        endian = '>'
    else:
        endian = '<'

    fou = open(sys.argv[2], 'wb')
    with open(sys.argv[1]) as fin:
        for i in range(6):  # _cell_length_{x, y, z}, _cell_angle_{alpha, beta, gamma}
            swap(fin, fou, 8)
        ngridx = struct.unpack(endian+"i", swap(fin, fou, 4))[0]
        ngridy = struct.unpack(endian+"i", swap(fin, fou, 4))[0]
        ngridz = struct.unpack(endian+"i", swap(fin, fou, 4))[0]
        zntype = struct.unpack(endian+"i", swap(fin, fou, 4))[0]
        for i in range(3):  # lewald, ltailc, lshift
            swap(fin, fou, 4)
        swap(fin, fou, 8)  # rcut
        fou.write(fin.read(zntype*128))
        for i in xrange(ngridx*ngridy*3):
            for j in xrange(ngridz*zntype*3):
                swap(fin, fou, 8)
    fou.close()
