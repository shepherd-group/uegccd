#!/bin/bash

# Allow for -d flag to compile for debugging.
# Further flags can be added to the case statement if desired.
debug=false
while getopts d flag; do
    case "${flag}" in
        d) debug=true;;
    esac
done

if [[ "$OSTYPE" == "darwin"* ]]; then
    export F95='gfortran -framework accelerate'
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    export F95='gfortran'
fi

if $debug; then
    F95+=' -g'
    echo "NOTE: compiling with debug symbols!"
fi

# Compile the fundamentals
$F95 -O3 -c         Precision.f90
echo "compiled Precision.f90"

$F95 -O3 -c         Constants.f90
echo "compiled Constants.f90"

$F95 -O3 -c      Types.f90
echo "compiled Types.f90"

$F95 -O3 -c      Pointers.f90
echo "compiled Pointers.f90"

# These two are needed for the IO routines to function
#$F95 -O3 -c -F2003 -fbounds-check  -Wall -Wextra HEG.f90
$F95 -O3 -c      HEG.f90
echo "compiled HEG.f90"

$F95  -c         IO.f90
echo "compiled IO.f90"

# These do the CCD
$F95 -O3 -c         MP2.f90
echo "compiled MP2.f90"

$F95 -O3 -c         CCScr.f90
echo "compiled CCScr.f90"

$F95 -O3 -c         DIISCC.f90
echo "compiled DIISCC.f90"

#$F95 -O3 -c -F2003 -fbounds-check  -Wall -Wextra CCD.f90
$F95 -O3 -c         CCD.f90
echo "compiled CCD.f90"

$F95 -O3 -c         dRPA.f90
echo "compiled dRPA.f90"

# Just compiling individual functions here
$F95 -O3 -c         AllCase.f90
echo "compiled AllCase.f90"

$F95 -O3 -c         Wrappers.f90
echo "compiled Wrappers.f90"

# About time we got to main!
$F95 -O3 -c         Main.f90
echo "compiled Main.f90"

if [[ "$OSTYPE" == "darwin"* ]]; then 
    $F95 -O3 *.o -o ZCode
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    $F95 -O3 *.o -o ZCode -l openblas
fi

#mkdir -p ZRun
#mv ZCode ZRun

#rm -rf ../bin

mkdir -p ../bin
cp ZCode ../benchmarks/uegccd
mv ZCode "../bin/uegccdv2.09"
rm *.o *.mod
