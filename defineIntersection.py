def write_circles():
    r = rcut
    coords = []
    for ii, my_center in enumerate(sphere):
        h, k, l = pbc(my_center[0],a),pbc(my_center[1],b), pbc(my_center[2],c)
        for theta in np.linspace(0,2*np.pi):
            for phi in np.linspace(0,np.pi):
                x = h+r*np.cos(theta)*np.sin(phi)
                y = k+r*np.sin(theta)*np.sin(phi)
                z = l+r*np.cos(phi)
                coords.append( 
                        list(map(pbc,[x,y,z],[a,b,c]))
                        )
    data = {'atoms':['F' for i in range(len(coords))],
                   'coords':coords}
    writer.xyz('spheres.xyz',data)
#           my_ax.plot(coords[:,0],coords[:,1],'-.',color='red')
#           if len(pbc_coords) > 0:
#               pbc_coords = [[x[i],y[i]] for i in pbc_coords]
#               coords  = fold_coords(pbc_coords,[a,b])
#               my_ax.plot(coords[:,0],coords[:,1],'-.',color='red')

from MCFlow.file_formatting import reader, writer
from MCFlow.calc_tools import calculate_distance, pbc, fold_coords
import numpy as np
import pprint

rcut = 5.0/2.
a, b, c = [20.022, 19.899, 13.383]
center = [10.011,4.9748,0.]
x, y, z = center
sphere = [  
            [x,y,z],
            [a/2-x,b/2+y,c/2+z],
            [x,b/2-y,z],
            [a/2+x,y,c/2-z],
            [a-x,b-y,c-z],
            [a/2+x,b/2-y,c/2-z],
            [a-x,b/2+y,c-z],
            [a/2-x,b-y,c/2+z]
        ]

if __name__ == '__main__':    
    
    new_data  = {'atoms':[],'coords':[]}
    old_data = reader.xyz('energy_grid_114_O_merged1.5.xyz')
    for atom, xyz in zip(old_data['atoms'], old_data['coords']):
        for my_center in sphere:
            if calculate_distance(xyz,my_center,[a,b,c]) < rcut:
                atom = 'H'
        new_data['atoms'].append(atom)
        new_data['coords'].append(xyz)
    
    
    writer.xyz('new_channels.xyz',new_data)
