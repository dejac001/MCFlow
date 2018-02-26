import sys, os

assert ((len(sys.argv) == 3) and
        (os.path.isfile( sys.argv[1]) )), 'usage: python3 %s file.db feed'%sys.argv[0]

file, feed_test = sys.argv[1], sys.argv[2]
import shelve, pprint

print('Here is what is in your db....')
db = shelve.open(file)
pprint.pprint(db[feed_test])
#print(db[feed_test])
