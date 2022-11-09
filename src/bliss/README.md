Bliss
-----

bliss is an open-source tool for computing
canonical labelings and automorphism groups of graphs.

This is a github copy of the original software,
for more information about bliss, please see the bliss homepage at
https://users.aalto.fi/tjunttil/bliss

Compiling
---------

In Linux and macOS, one can use GNU Make to compile the "bliss" executable,
as well as the static and shared libraries, with

make -f Makefile-manual

If you are embedding bliss in a project that needs the automorphism group sizes
as GNU Multiple Precision Arithmetic Library integers,
compile with

make -f Makefile-manual gmp


On Linux and macOS with CMake installed, one can also use

cmake .
cmake --build .

To enable GNU GMP support, compile with

cmake -D USE_GMP=ON .
cmake --build .


On Windows with Visual Studio and CMake, use

cmake .
cmake --build . --config Release


Examples
--------

To see how bliss can be called from C++,
see the bliss main routine at "src/bliss.cc" and
the file in the "examples" directory.
