# install dependencies
apt-get update
apt-get install -y cmake g++ git clang

# build soplex
git clone https://github.com/scipopt/soplex.git
cd soplex
git checkout bugfix-60
mkdir build
cd build
cmake ..
make -j
cd ../..

# build scip in `install` directory
mkdir build
cd build
cmake .. -DAUTOBUILD=ON -DSOPLEX_DIR=soplex/build -DCMAKE_INSTALL_PREFIX=../install
make -j
make install
cd ..