#!/bin/bash

set -euo pipefail

# brew install bison boost ipopt ccache

ROOT=$(pwd)/..
# -DPAPILO=ON -DPAPILO_DIR=$ROOT/scipoptsuite-8.0.3/papilo/
cmake -DCMAKE_INSTALL_PREFIX=$ROOT/deps-installation -S$ROOT/soplex -B$ROOT/soplex/build
cmake --build $ROOT/soplex/build -j 16
cmake --install $ROOT/soplex/build

cmake -DCMAKE_INSTALL_PREFIX=$ROOT/deps-installation -S$ROOT/papilo -B$ROOT/papilo/build
cmake --build $ROOT/papilo/build -j 16
cmake --install $ROOT/papilo/build

cmake -DCMAKE_INSTALL_PREFIX=$ROOT/deps-installation -S$ROOT/scipoptsuite-8.0.3/zimpl -B$ROOT/scipoptsuite-8.0.3/zimpl/build
cmake --build $ROOT/scipoptsuite-8.0.3/zimpl/build -j 16
cmake --install $ROOT/scipoptsuite-8.0.3/zimpl/build

# Build scip like this: cmake -DCMAKE_PREFIX_PATH=$ROOT  -DLPS=pdlp -DCMAKE_CXX_COMPILER_LAUNCHER=ccache ..
