#!/bin/bash

# exit on error and don't allow the use of unset variables
set -o errexit
#set -o nounset

SRC_DIR=`dirname $0`
CP=""
for i in `find "${SRC_DIR}/dist" -name '*.jar'`; do CP="${CP}:${i}"; done

java -enableassertions -Xmx1g -cp "${CP}" org.jax.bioinfdata.genocall.ConvertGenotypeHDF5ToFlatFileMain $@

