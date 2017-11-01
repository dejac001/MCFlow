def FloydWarshallWithPath(weights):
    '''
    '''
    dist = np.ones(np.shape(weights[:]))
    my_next = np.ones(np.shape(weights[:]))*-1
    N, M = np.shape(weights)
    assert N == M, 'Weights matrix not square'
    for u in range(N):
        for v in range(N):
            dist[u,v] = weights[u,v]
            my_next[u,v] = v
    for k in range(N):
        for i in range(N):
            for j in range(N):
                if dist[i,j] > dist[i,k] + dist[k,j]:
                    dist[i,j] = dist[i,k] + dist[k,j]
                    my_next[i,j] = my_next[i,k]
    return dist, my_next

def Path(u,v,nexxt):
    if nexxt[u,v] == -1:
        return []
    path = []
    while u != v:
        u = nexxt[u,v]
        path.append(u)
    return path

def test():
    weights = np.matrix([
    [0,np.inf,-2,np.inf],
    [4, 0, 3, np.inf],
    [np.inf, np.inf, 0, 2],
    [np.inf, -1, np.inf, 0]])
    return weights

import numpy as np

if __name__ == '__main__':
    w = test()
    d,n = FloydWarshallWithPath(w)
    print(d)
    print(Path(2,1,n))
    print(Path(2,0,n))
    print(Path(0,1,n))
