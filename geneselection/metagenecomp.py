#!/usr/bin/env python

import csv
import getopt,sys
import operator
from scipy import *
from numpy import *
import numpy as np
from scipy.io import loadmat
from collections import OrderedDict
import argparse

def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', required=True,
                        help='correlation file name')
    parser.add_argument('-t', '--type', required=True,
                        help='amp|del')
    parser.add_argument('-at', '--avg_thresh', type=float, required=True,
                        help='average threshold')
    parser.add_argument('-et', '--extr_thresh', type=float, 
                        help='max or min threshold')
    parser.add_argument('-gr', '--gene_ranker_file', required=True,
                        help='gene ranker filename')
    parser.add_argument('-org', '--output_ranked_genes', required=False,
                         help='output file for ranked genes')
    parser.add_argument('-oug', '--output_unranked_genes', required=False,
                        help='output file for unranked genes')
    parser.add_argument('-gf', '--gistic_file', required=False,
                        help='gistic results file')
    parser.add_argument('-nf', '--not_found_file', required=False,
                        help='genes not found in gistic')

    return parser.parse_args()

def load_gene_data(filename):
    f = open(filename,'rb')
    csvreader=csv.reader(f)
    filearray = []
    filearray.extend(csvreader)
    f.close()

    return filearray

def load_gistic_data(filename):
    f = open(filename,'rb')
    csvreader=csv.reader(f, delimiter='\t')
    filearray = []
    filearray.extend(csvreader)
    f.close()

    return filearray

def read_gene_ranker(filename):
    f = open(filename,'rb')
    csvreader=csv.reader(f, delimiter=',')
    filearray = []
    filearray.extend(csvreader)
    f.close()

    return filearray

def load_corr_data(filename):
    f = open(filename,'rb')
    data = loadmat('FilteredProbes.mat')
    probes = data['FilteredProbes']
    csvreader=csv.reader(f, delimiter=',')
    filearray = []
    filearray.extend(csvreader)
    f.close()

    return [filearray,probes]

def write_ranked(filename, genes):
    f = open(filename,'wb')
    csvwriter =csv.writer(f, delimiter=',')
    csvwriter.writerow(["gene", "chr", "avg", "max", "min", "rho", "gene_rank","score","gene_lists"])
    for gene in genes:
        csvwriter.writerow(genes[gene])

    f.close()

def write_unranked(filename, genes):
    f = open(filename,'wb')
    csvwriter =csv.writer(f, delimiter=',')
    csvwriter.writerow(["gene", "chr", "avg", "max", "min", "rho"])
    for gene in genes:
        csvwriter.writerow(genes[gene])

    f.close()
"""
program to analyze amplification and deletion for a csv file
containing the following records:
    row,col,zscore,rho,ge_chr,ge_gene,cnv_chr,cnv_gene 
    args: -f <file name>, -gr <gene ranker file> -t type [amp|del], -at <avg threshold> -et <max or min threshold>
"""
if __name__ == "__main__":

    args = arg_parser()
    (filearray,probes) = load_corr_data(args.file)
    generanks = read_gene_ranker(args.gene_ranker_file)

    avg_logr = np.array([])
    max_logr = np.array([])
    min_logr = np.array([])

    for row,col,zscore,rho,ge_chr,ge_gene,cnv_chr,cnv_gene in filearray:
        avg_logr=append(avg_logr, array(mean(probes[:,col])))
        max_logr=append(max_logr, array(max(probes[:,col])))
        min_logr=append(min_logr, array(min(probes[:,col])))
    
    ranked_genes = {}
    unranked_genes = {}
    avg_avg = sort(avg_logr)[::-1]
    indices=argsort(avg_logr)[::-1]
    for (i,item) in enumerate(avg_avg):
        if args.type == "del":
            if item < args.avg_thresh and min_logr[i] < args.extr_thresh:
                row  = filearray[indices[i]]
                rho  = row[3]
                chr  = row[4]
                gene = row[5]
                if gene not in ranked_genes and gene not in unranked_genes:
                    generank = filter(lambda x: x[2] == gene, generanks)
                    if len(generank) > 0:
                        ranked_genes[gene] = [gene, int(chr), item, max_logr[i], min_logr[i], rho]
                        ranked_genes[gene].extend([int(generank[0][0]),float(generank[0][1]),int(generank[0][3])] )
                    else:
                        unranked_genes[gene] = [gene, int(chr), item, max_logr[i], min_logr[i], rho]
                        #print gene, " NOT in gene list"
        elif args.type == "amp":
            if item > args.avg_thresh and max_logr[i] > args.extr_thresh:
                row  = filearray[indices[i]]
                rho  = row[3]
                chr  = row[4]
                gene = row[5]
                if gene not in ranked_genes and gene not in unranked_genes:
                    generank = filter(lambda x: x[2] == gene, generanks)
                    if len(generank) > 0:
                        ranked_genes[gene] = [gene, int(chr), item, max_logr[i], min_logr[i], rho]
                        ranked_genes[gene].extend([int(generank[0][0]),float(generank[0][1]),int(generank[0][3])] )
                    else:
                        unranked_genes[gene] = [gene, int(chr), item, max_logr[i], min_logr[i], rho]
                        #print gene, " NOT in gene list"
        
    if args.output_unranked_genes:
        write_unranked(args.output_unranked_genes, OrderedDict(sorted(unranked_genes.items(), key=lambda x: x[1][5])))

    if args.output_ranked_genes:
        write_ranked(args.output_ranked_genes, OrderedDict(sorted(ranked_genes.items(), key=lambda x: x[1][6] )))

    if args.gistic_file:
        gistic_data = load_gistic_data(args.gistic_file)
        amplified_gistic_genes = [x[0] for x in gistic_data if x[1] == 'Amp']
        del_gistic_genes = [x[0] for x in gistic_data if x[1] == 'Del']
        if args.type == 'amp':
            for gene in ranked_genes:
                if gene not in amplified_gistic_genes:
                    print gene
        elif args.type == 'del':
            for gene in ranked_genes:
                if gene not in del_gistic_genes:
                    print gene

    if args.not_found_file:
        genes_not_found = load_gene_data(args.not_found_file)
        for gene in genes_not_found:
            generank = filter(lambda x: x[2] == gene, generanks)
            if len(generank) > 0:
                print generank
    #for key,values in sorted(ranked_genes.items(), key=lambda x: x[1][6]):
        #print values

    #print "\n"

    #for key,values in sorted(unranked_genes.items(), key=lambda x: x[1][0]):
        #print values
