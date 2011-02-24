#!/bin/bash

./fftohdf5.bash -genoinfile ~/projects/emma-scripts/PopulationData/popdataSep09schr1.csv -aallelecol 2 -ballelecol 3 -snpcol 1 -chrcol 4 -poscol 5 -bpbuild NCBIBuild37 -firstgenocol 6 -hdf5out output1.h5

