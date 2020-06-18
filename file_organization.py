import os


def rob(path, name, type, num):
    return os.path.join(path, '%s%i' % (type, num), '%s%s%i' % (name, type, num))


def tjo(path, name, type, num):
    return os.path.join(path, '%s%s%i' % (name, type, num))


# note: need to change file here based off of how you organize your files
# if os.getlogin() == 'tjo':
#     equilName='all.equil'
#     prodName='all.prod'
#     read = tjo
# else:
equilName = 'equil-'
prodName = 'prod-'
read = rob
