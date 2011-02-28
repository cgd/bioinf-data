#!/bin/bash

./fftohdf5.bash -genoinfiles ~/projects/emma-scripts/PopulationData/popdataSep09schr[[:digit:]]*.csv -aallelecol 2 -ballelecol 3 -snpcol 1 -chrcol 4 -poscol 5 -bpbuild NCBIBuild37 -firstgenocol 6 -hdf5out output.h5

