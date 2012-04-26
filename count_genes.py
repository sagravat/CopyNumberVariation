#!/usr/bin/env python
import csv

print "opening file"
csvfile = open("chr_filtered.csv", "rb")
print "opened file"
reader=csv.reader(csvfile,delimiter=",")
print "read file"

gene_map = {}

#9214, 6064, 33.456454, 0.863746, 1, PRELP, 1, ATP2B4
i = 1
for row, col, zscore, rho, ch1, p_gene, ch2, gene in reader:
    row,col,zscore
    gene_map[gene] = gene
    
csvfile.close()
print "close file"
print "length: ", len(gene_map)

