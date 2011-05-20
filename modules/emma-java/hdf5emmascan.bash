#!/bin/bash

# exit on error and don't allow the use of unset variables
set -o errexit
#set -o nounset

for jar_file in `find dist -name '*.jar'`; do CP="${CP}:${jar_file}"; done

java -d32 -enableassertions -Xmx1g -cp "${CP}" org.jax.emma.Hdf5EmmaMain $@

