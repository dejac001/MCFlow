def getNbox(data, box, feed):
    for key, value in data.items():
        if key.isdigit() and box in value.keys():
            if value[box]['mean'] > 0.:
                print('feed: %s    mol:%s in %s is %s +/- %s'%(feed, key, box, value[box]['mean'],
                                                      value[box]['stdev']))
        else:
            if key != 'time' and len(key.split()) == 1:
                getNbox(data[key], box, feed)

def getNmolBox(data, box, mol):
    if mol in data.keys():
        return data[mol][box]['mean'], data[mol][box]['stdev']
    for key, value in data.items():
        if key != 'time':
            return getNmolBox(data[key], box, mol)


import sys, os

if __name__ == '__main__':
    assert len(sys.argv) == 3, 'usage: python3 %s box[1,2,3] feed'%sys.argv[0]
    
    box, feed_test = sys.argv[1], sys.argv[2]
    import shelve
    
    db = shelve.open('N-data.db')
    getNbox(db[feed_test],box, feed_test)
