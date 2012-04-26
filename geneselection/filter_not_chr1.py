#!/usr/bin/env python
import csv

csvfile = open("diff_chr_corr_no_colrows.csv", "rb")
print "open file"
reader=csv.reader(csvfile)
print "read file"

data = []
data.extend(reader)
csvfile.close()
print "close file"
t = filter(lambda x: int(x[3]) != 1, data)
print "filtered"

f = open("not_Chr1.csv", "wb")
writer = csv.writer(f)
writer.writerows(t)
f.close()
