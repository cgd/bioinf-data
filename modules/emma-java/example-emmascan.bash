#!/bin/bash

# exit on error and don't allow the use of unset variables
set -o errexit
#set -o nounset
set +x

echo "scanning the CSV directly"
./emmascan.bash \
    -aallelecol 2 -ballelecol 3 -firstgenocol 6 \
    -genofile ~/projects/emma-scripts/PopulationData/popdataSep09schr19.csv \
    -phenofile ~/projects/emma-scripts/bone-mineral-density.txt -out chr19-csv-scanout.txt

echo "converting CSV to HDF5 format"
../main/fftohdf5.bash \
    -genoinfiles ~/projects/emma-scripts/PopulationData/popdataSep09schr19.csv \
    -aallelecol 2 -ballelecol 3 \
    -snpcol 1 -chrcol 4 -poscol 5 -bpbuild NCBIBuild37 -firstgenocol 6 -hdf5out popdata-chr19.h5

echo "scanning the HDF5"
./hdf5emmascan.bash \
    -genofile "popdata-chr19.h5" \
    -phenofile ~/projects/emma-scripts/bone-mineral-density.txt -out chr19-hdf5-scanout.txt
./hdf5emmascan.bash \
    -sex male \
    -genofile "popdata-chr19.h5" \
    -phenofile ~/projects/emma-scripts/bone-mineral-density.txt -out chr19-hdf5-scanout-male.txt
./hdf5emmascan.bash \
    -sex female \
    -genofile "popdata-chr19.h5" \
    -phenofile ~/projects/emma-scripts/bone-mineral-density.txt -out chr19-hdf5-scanout-female.txt

