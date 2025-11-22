#!/bin/bash
set -xue

docker pull quay.io/pypa/manylinux2014_x86_64
docker run --rm -v `pwd`:/io quay.io/pypa/manylinux2014_x86_64 bash /io/build-wheels.sh

# Not yet an actual release script, but building the wheels is a step along 
# the way. Will be expanded.
#
# twine upload dist/*-manylinux2014_x86_64.manylinux_2_17_x86_64.whl --repository=testpypi --verbose
# twine upload dist/*-manylinux2014_x86_64.manylinux_2_17_x86_64.whl --repository=pypi     --verbose
