#!/usr/bin/env python
import csv

csvfile = open("chr_filtered.csv", "rb")
print "open file"
reader=csv.reader(csvfile)
print "read file"

data = []
data.extend(reader)
csvfile.close()
print "close file"
t = filter(lambda x: x[4] != x[6], data)
print "filtered"

f = open("diff_chr_corr.csv", "wb")
writer = csv.writer(f)
writer.writerows(t)
f.close()
