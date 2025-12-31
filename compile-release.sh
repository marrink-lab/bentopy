#!/bin/bash
set -xue

docker pull quay.io/pypa/manylinux_2_34
docker run --rm -v `pwd`:/io quay.io/pypa/manylinux_2_34 bash /io/build-wheels.sh
