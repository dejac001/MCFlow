import math

a = '''
# C(g/mL)    dG(kJ/mol)    dC     dG\n
 1.228650e-04 -4.019655e+01 5.569724e-05 1.263422e+00 #4mol-8\n
 1.821449e-04 -3.890385e+01 6.153792e-05 1.321169e+00 #48mol-8\n
 1.373585e-04 -3.809177e+01 5.255584e-05 1.724485e+00 #32mol-8\n
 1.949738e-04 -3.889914e+01 7.598615e-05 1.332330e+00 #64mol-8\n
# C(g/mL)    dG(kJ/mol)    dC     dG\n
2.199527e-02 -3.947201e+01 2.358939e-03 6.456859e-01 #8mol-8\n
6.297941e-03 -3.998389e+01 9.060157e-04 5.658077e-01 #4mol-8\n
3.856303e-02 -3.988261e+01 4.006795e-03 5.228218e-01 #16mol-8\n
#7.407352e-02 -4.093807e+01 6.444124e-03 4.575015e-01 #32mol-8\n
#8.318884e-02 -4.289113e+01 9.121247e-03 5.903188e-01 #128mol-8\n
'''
means = []
stdev = []
for line in a.split('\n'):
    if (not line.startswith('#')) and (len(line.split()) == 5):
        means.append(float(line.split()[1]))
        stdev.append(float(line.split()[3]))
        if len(means) == 1:
            print(line, means, stdev)
print('%i data points'%len(means))

w = [1/pow(i,2) for i in stdev]

mean = sum([xi*wi for (xi, wi) in zip(means,w)])/sum(w)
stdev = 1/sum(w)*1.96 # 1.96 b/c large number of observation limit for weighted fit

print('dG: %f +/- %f kJ / mol'%(mean, stdev))
# calculate kH
a = 104.149/(8314*323.)
b = -1/(8.314/1000*323.)
kH = a*math.exp(b*mean)
kH_error = abs(kH*b*stdev)
print('kH: %f +/- %f g / mL / kPa'%(kH, kH_error))
