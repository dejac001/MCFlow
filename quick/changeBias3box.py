def addBias(n_des, vals):
    mean, std = np.mean(vals), np.std(vals)
    dG = T*np.log(mean/n_des)
    ddG = np.sqrt(pow(T*std/mean,2))
    return dG, ddG

a = '''
1.688 0.264 0.048
1.441 0.558 0.001
1.088 0.845 0.067
1.582 0.183 0.235
1.162 0.837 0.001
1.866 0.132 0.002
1.921 0.077 0.002
1.627 0.295 0.079
'''
T = 348.
import numpy as np


if __name__ == '__main__':

    x = []
    for i in a.split('\n'):
        if len(i) > 0:
            x.append(list(map(float,i.split())))

    x = np.matrix(x)
    for box in [0,1,2]:
        n_desired = np.sum(x[0,:])/3
        dG, ddG = addBias(n_desired, x[:,box])
        print('add %f +/- %f to bias in box%i'%(dG,ddG,box+1))
