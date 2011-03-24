#!/bin/bash

# exit on error and don't allow the use of unset variables
set -o errexit
#set -o nounset

# the following must be used on Mac if a 64-bit JVM is the default
JAVA_HOME=/System/Library/Frameworks/JavaVM.framework/Versions/1.5/Home/

java -enableassertions -Xmx1g -jar `dirname $0`/dist/emma-java-1.0.jar $@

