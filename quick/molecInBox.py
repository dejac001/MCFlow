def getNbox(data, box):
    for key, value in data.items():
        if key.isdigit() and box in value.keys():
            if value[box]['mean'] > 0.:
                print('mol:%s in %s is %s +/- %s'%(key, box, value[box]['mean'],
                                                      value[box]['stdev']))
        else:
            if key != 'time':
                getNbox(data[key], box)


import sys, os

if __name__ == '__main__':
    assert len(sys.argv) == 3, 'usage: python3 %s box[1,2,3] feed'%sys.argv[0]
    
    box, feed_test = sys.argv[1], sys.argv[2]
    import shelve
    
    db = shelve.open('N-data.db')
    getNbox(db[feed_test],box)
