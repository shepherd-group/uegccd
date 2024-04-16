#!/bin/bash 
./make
cd ZRun
./ZCode >/dev/null
grep E\( Output
grep Ec\( Output
