#!/usr/bin/env python
import csv

csvfile = open("chr_sorted.csv", "rb")
print "open file"
reader=csv.reader(csvfile)
print "read file"

data = []
data.extend(reader)
csvfile.close()
print "close file"
t = filter(lambda x: x[5] != x[7], data)
print "filtered"

f = open("chr_filtered.csv", "wb")
writer = csv.writer(f)
writer.writerows(t)
f.close()
