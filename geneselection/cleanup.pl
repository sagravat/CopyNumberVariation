#!/usr/bin/env perl

open(FILE,"<tcga-generanker.csv");
$index = 1;
while (<FILE>) {
    my ($line) = $_;
    chomp($line);
    if ($index % 11 != 0) {
        print "$line\n"
    }
    $index++;
}
close(FILE);
