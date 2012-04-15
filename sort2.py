#!/usr/bin/env python

import time
import numpy as np
import csv

start = time.time()
def elapsed():
    return time.time() - start

# count data rows, to preallocate array
f = open('links.csv', 'rb')
def count(f):
    while 1:
        block = f.read(65536)
        if not block:
             break
        yield block.count(',')

linecount = sum(count(f))
print '\n%.3fs: file has %s rows' % (elapsed(), linecount)

# pre-allocate array and load data into array
m = np.zeros(linecount, dtype=[('a', np.uint32), ('b', np.uint32)])
f.seek(0)
f = csv.reader(open('links.csv', 'rb'))
for i, row in enumerate(f):
    m[i] = int(row[0]), int(row[1])

print '%.3fs: loaded' % elapsed()
# sort in-place
m.sort(order='b')

print '%.3fs: sorted' % elapsed()
