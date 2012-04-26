#!/bin/sh

total=18743
for i in {1..22}
do
    count=`cut -d"," -f2 same_corr.csv | sort | grep "^${i}$" | wc -l`
    echo "$count/$total"
done
