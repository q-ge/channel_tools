#!/bin/bash

myroot=$(dirname $0)

root=$1
chip=$2
chan=$3
cm=$4
ts=$5
count=$6

colrange=$(grep "^modulation range" $root/$chip/$chan/info | \
           sed "s/.*:\W*//" | sed "s/-/ /")
rowrange=$(grep "^result range" $root/$chip/$chan/$cm/limits | \
           sed "s/.*:\W*//" | sed "s/-/ /")

matrix="$chip.$chan.$cm.$ts.$count.cm"

find $root/$chip/$chan/$cm/TS_$ts -name "*.xz" | xargs xzcat | \
    $myroot/filter_samples $colrange $rowrange |
    ./channel_matrix $matrix $colrange $count
