def addBias(n_des, vals):
    mean, std = np.mean(vals), np.std(vals)
    dG = T*np.log(mean/n_des)
    ddG = np.sqrt(pow(T*std/mean,2))
    return dG, ddG

a = '''
'''
T = 323.
import numpy as np


if __name__ == '__main__':
    print('Temperature is',T)

    x = []
    for i in a.split('\n'):
        if len(i) > 0:
            x.append(list(map(float,i.split())))

    x = np.matrix(x)
    for box in [0,1,2]:
        n_desired = np.sum(x[0,:])/3
        dG, ddG = addBias(n_desired, x[:,box])
        print('add %f +/- %f to bias in box%i'%(dG,ddG,box+1))
