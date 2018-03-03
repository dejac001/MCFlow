a = '''
7.672 11.690 0.638
9.783 10.109 0.108
9.000 10.996 0.004
7.999 11.726 0.275
7.828 11.399 0.772
7.471 12.030 0.499
7.366 9.084 3.550
9.830 10.000 0.170
'''


if __name__ == '__main__':
    import numpy as np

    x = []
    for i in a.split('\n'):
        if len(i) > 0:
            x.append(list(map(float,i.split())))

    x = np.matrix(x)
    for box in [0,1,2]:
        dG = -323.*np.log(1/np.mean(x[:,box]))
        print('add %f to bias in box%i'%(dG,box+1))
