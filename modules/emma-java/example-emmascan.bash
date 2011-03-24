#!/bin/bash

# exit on error and don't allow the use of unset variables
set -o errexit
#set -o nounset

./emmascan.bash -aallelecol 2 -ballelecol 3 -firstgenocol 6 -genofile ~/projects/emma-scripts/PopulationData/popdataSep09schr19.csv -phenofile ~/projects/emma-scripts/bone-mineral-density.txt -out scanout.txt

