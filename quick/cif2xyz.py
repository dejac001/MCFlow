import argparse, plotter
from pymatgen.io.cifio import CifParser
out_file = '/Users/dejacor/Documents/screening'
parser = argparse.ArgumentParser(description='convert zeolite unit cell from cif to xyz format')
parser.add_argument('-c','--cifFile',help='zeolite pdb file name')
args = parser.parse_args()
cif_file = args.cifFile
cif = CifParser(cif_file)
struc = cif.get_structures(primitive=False)[0]

DATA = {}
for atom in range(len(struc)):
    xyz = struc[atom].as_dict()['xyz']
    element = struc[atom].as_dict()['label']
    if not element in DATA.keys():
        DATA[element] = [xyz]
    else:
        DATA[element].append(xyz)

# gotta correct data, for some reason they pull in my Si's as S's
corData = {}
for element in DATA.keys():
    if element == 'S':
        corData['Si'] = DATA['S']
    else:
        corData[element] = DATA[element]
name_of_zeo = list(cif.as_dict().keys())[0]
xyz_file = cif_file[:cif_file.rfind('/')]+ '/%s.xyz' % name_of_zeo
plotter.xyz(xyz_file, corData)
