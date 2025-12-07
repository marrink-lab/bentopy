#!/bin/bash
set -ex

# Note: The inscrutible *[!t] glob serves to not build against any
# free-threading python versions (e.g., python 3.13t), because that is
# currently giving failed builds. This probably has something to do with the
# cython support for *t python versions, but I'm not entirely sure.
# (Marieke, 2025-03-05)

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- --profile=minimal -y
export PATH="$HOME/.cargo/bin:$PATH"

# Compile wheels
for PYBIN in /opt/python/cp{312,313,314}*[!t]/bin; do
    rm -rf /io/build/
    "${PYBIN}/pip" install -U setuptools setuptools-rust
    "${PYBIN}/pip" wheel /io/ -w /io/dist/ --no-deps
done

# Bundle external shared libraries into the wheels
for whl in /io/dist/*{cp312,cp313,314}*[!t].whl; do
    auditwheel repair "$whl" -w /io/dist/
done
