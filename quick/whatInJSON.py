import sys, os

assert ((len(sys.argv) == 3) and
        (os.path.isfile( sys.argv[1]) )), 'usage: python3 %s file.json feed'%sys.argv[0]

file, feed_test = sys.argv[1], sys.argv[2]
import json, pprint

print('Here is what is in your json file....')

with open(file, 'r') as f:
    data = json.load(f)

pprint.pprint(data[feed_test])
