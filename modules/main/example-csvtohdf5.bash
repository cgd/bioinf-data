#!/bin/bash

echo "creating output with alleles"
./fftohdf5.bash \
    -genoinfiles ~/projects/emma-scripts/PopulationData/popdataSep09schr[[:digit:]]*.csv \
    -aallelecol 2 -ballelecol 3 \
    -snpcol 1 -chrcol 4 -poscol 5 -bpbuild NCBIBuild37 -firstgenocol 6 -hdf5out output.h5

#echo "creating output without alleles"
#./fftohdf5.bash \
#    -genoinfiles ~/projects/emma-scripts/PopulationData/popdataSep09schr[[:digit:]]*.csv \
#    -snpcol 1 -chrcol 4 -poscol 5 -bpbuild NCBIBuild37 -firstgenocol 6 -hdf5out output-without-alleles.h5

