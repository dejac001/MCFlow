def FloydWarshallWithPath(weights):
    '''
    '''
    dist = np.ones(np.shape(weights))
    prev = np.ones(np.shape(weights),dtype=int)*-9999
    N, M = np.shape(weights)
    assert N == M, 'Weights matrix not square'
    for u in range(N):
        for v in range(N):
            dist[u,v] = weights[u,v]
            prev[u,v] = u
    for k in range(N):
        for i in range(N):
            for j in range(N):
                if dist[i,j] > dist[i,k] + dist[k,j]:
                    dist[i,j] = dist[i,k] + dist[k,j]
                    prev[i,j] = prev[k,j]
    return dist, prev

def Path_new(u,v, prev, path=[]):
    if u != v:
        v = prev[u,v]
        return Path_new(u,v,prev, path=[v] + path)
    else:
        return path

def test():
    weights = np.matrix([
    [0,np.inf,-2,np.inf],
    [4, 0, 3, np.inf],
    [np.inf, np.inf, 0, 2],
    [np.inf, -1, np.inf, 0]])
    return weights

def test1():
    w = test()
    return FloydWarshallWithPath(w)

def test2():
    weights = test()
    return sps.csgraph.floyd_warshall(weights,return_predecessors=True)

import scipy.sparse as sps
import numpy as np

if __name__ == '__main__':
    d,n1 = test1()
    d,n2 = test2()
    u,v=0,1
    test_mine = Path_new(u,v,n1)
    test_yours = Path_new(u,v,n2)
    print(test_mine + [v])
    print(test_yours + [v])
