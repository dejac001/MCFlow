def averageG(g, rMax, dr):
    edges = np.arange(0., rMax + 1.1*dr, dr)
    num_bins = len(edges) - 1
    radii = np.zeros(num_bins)
    g_average = np.zeros(len(edges)-1)
    for i in range(len(edges)-1):
        rOuter = edges[i+1]
        rInner = edges[i]
        radii[i] = (rInner+rOuter)/2.
        g_average[i] = np.mean(g[:,i]) / (4.0/3.0*np.pi*(rOuter**3 - rInner**3))
    return g_average, radii


def pairCorrelationFunction_frame(coords1, coords2, boxLengths, rMax, dr):
    '''
    compute the three-dimensional pair correlsation function for a set of spherical particles contained
    in a cube with side length S. This simple function finds reference particles such that a sphere of
    radius rMax drawn around the particle will fit entirely within the cube, eliminating the
    need to compensate for edge effects. If no such particles exist, an error is returned.
    Try a smaller rMax...or write some code to handle edge effects! ;)
    :param coords1: a matrix of xyz coordinates of bead 1 [[x1,y1,z1]...[xN,yN,zN]]
    :param coords2: a matrix of xyz coordinates of bead 2 [[x1,y1,z1]...[xN,yN,zN]]
    :param S:  a matrix of boxlengths in each direction
    :param rMax:  outer diameter of largest spherical shell
    :param dr:  increment for increasing radius of spherical shell
    :return: returns a tuple: (g, radii, interior_indices)
        g(r)        a numpy array containing the correlation function g(r)
        radii       a numpy array containing the radii of the
                    spherical shells used to compute g(r)
        reference_indices   indices of reference particles
    '''
    def minImage(dist, boxl):
        return (
            dist - np.round( np.divide(dist,boxl) )*boxl
        )

    edges = np.arange(0., rMax + 1.1*dr, dr)
    num_bins = len(edges) -1


    n_coords1, n_coords2 = len(coords1[:,0]) , len(coords2[:,0])
    assert n_coords1 > 3, 'Wrong dimension {}'.format(n_coords1)
    num_particles = n_coords1 + n_coords2
    g = np.zeros([num_particles, num_bins])

    number_density = num_particles/(boxLengths[0]*boxLengths[1]*boxLengths[2])

    for index in range(n_coords1): #TODO: figure out if this summation is correct
        dx = coords1[index,0] - coords2[:,0]
        dy = coords1[index,1] - coords2[:,1]
        dz = coords1[index,2] - coords2[:,2]
        dx = minImage(dx, boxLengths[0])
        dy = minImage(dy, boxLengths[1])
        dz = minImage(dz, boxLengths[2])


        dist = np.sqrt(dx**2 + dy**2 + dz**2)
        (result, bins) = np.histogram(dist, bins=edges, normed=False)
        g[index,:] = result / number_density
    return g


import numpy as np
