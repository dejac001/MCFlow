import os


def rob(path, name, type, num):
    return '%s/%s%i/%s%s%i' % (path, type, num, name, type, num)


def tjo(path, name, type, num):
    return '%s/%s%s%i' % (path, name, type, num)


# note: need to change file here based off of how you organize your files
# if os.getlogin() == 'tjo':
#     equilName='all.equil'
#     prodName='all.prod'
#     read = tjo
# else:
equilName = 'equil-'
prodName = 'prod-'
read = tjo
