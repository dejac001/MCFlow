def kH(dG, MW):
    mean, stderr = dG
    b = -1/(0.008314*323.)
    a = MW/(8314*323.)
    return [
    math.exp(b*mean)*a,
    b*stderr*math.exp(b*mean)*a,
    ]

import math

dG = [-12.7,0.7]
if __name__ == '__main__':
    print(kH(dG, 102.162))
    print(kH([-39.7,0.2], 104.15))

