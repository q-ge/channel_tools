#!/bin/bash

scriptroot=$(dirname $0)
dataroot=$1
epsilon=$2

matrices=$(ls $dataroot/matrices/*.cm)

for m in $matrices
do
    echo "Finding bandwidth for $m to precision $epsilon"
    $scriptroot/capacity $m $epsilon -q > \
        $dataroot/matrices/$(basename $m).capacity
done
