#!/usr/bin/env python
import csv

csvfile = open("sorted_significant.txt", "rb")
print "open file"
reader=csv.reader(csvfile)
print "read file"

data = []
data.extend(reader)
csvfile.close()
print "close file"

sorted_chr = []

for i in range(1,22):
    print "before filter",i
    t = filter(lambda x: int(x[6]) == i, data)
    print "done filtered "
    s = sorted(t, key=lambda x: float(x[2]), reverse=True)
    print "sorted ", i
    sorted_chr.extend(s)
    print "extended ",i



f = open("chr_sorted.csv", "wb")
writer = csv.writer(f)
writer.writerows(sorted_chr)
f.close()
