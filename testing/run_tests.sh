#!/usr/bin/env bash

binary=`pwd`/../bin/uegccdv2.09

for path in \
    N_00014_0002 \
    N_00038_0004 \
    N_00054_0006 \
    N_00066_0008 \
    N_00114_0010; do

    echo "Running: $path"

    cd $path

    eval $binary

    cd ../

done
