#!/bin/bash

# exit on error and don't allow the use of unset variables
set -o errexit
#set -o nounset

for jar_file in `find dist -name '*.jar'`; do CP="${CP}:${jar_file}"; done

# the following must be used on Mac if a 64-bit JVM is the default
JAVA_HOME=/System/Library/Frameworks/JavaVM.framework/Versions/1.5/Home/

java -enableassertions -Xmx1g -cp "${CP}" org.jax.emma.CsvEmmaMain $@

