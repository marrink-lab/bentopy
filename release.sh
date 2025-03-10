#!/bin/bash
set -xue

docker pull quay.io/pypa/manylinux2014_x86_64
docker run --rm -v `pwd`:/io quay.io/pypa/manylinux2014_x86_64 bash /io/build-wheels.sh
# Not yet an actual release script, but building the wheels is a step along 
# the way. Will be expanded.
