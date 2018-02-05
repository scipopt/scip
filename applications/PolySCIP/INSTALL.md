Building and installing the PolySCIP application               {#INSTALL}
================================================

The PolySCIP application may be built in two ways, either by using Cmake
or by using SCIP's own Makefile system. The choice of the
build system depends on the desired target platform and the availability
of those tools there, as well as your personal preferences.
For most users, we recommend to use the CMake system, especially on
non Unix platforms.


Installation information using CMake
------------------------------------

Please compile SCIP first,
see the INSTALL_CMAKE in the main SCIP directory for instructions,
or refer to the online documentation of SCIP.
The PolySCIP application can be compiled within the same build directory
as SCIP. Assuming that the build directory of SCIP was named "build",
simply execute

```
cd build
make polyscip
```

The PolySCIP application is part of the SCIP applications. To build all
applications at once, use

```
make applications
```

It is also possible to build the PolySCIP application in a stand-alone
build directory. Therefore, it is necessary to create the
stand-alone build directory first and generate the Makefile using
CMake. It might be necessary to specify the SCIP build directory
or installation directory, if SCIP has not yet been installed systemwide.

```
mkdir build
cd build
cmake .. [-DSCIP_DIR=../../]
make
```

Please refer to the [online documentation of SCIP](http://scip.zib.de/doc/html/CMAKE.php)
for a list of available
configuration options for this application, and available tests.


Installation information for SCIP's custom Makefile system on Linux
-------------------------------------------------------------------


1) Install SCIP
(see INSTALL file in SCIP directory)

2) Install PolySCIP
PolySCIP is equipped with a Makefile located in the PolySCIP
directory. Change on the command line into the PolySCIP directory and
execute 'make'

In case of compilation problems, try torun 'make' with the same command
line options that were used for the compilation of  SCIP.
For instance, if you have successfully built SCIP via 'make COMP=gnu
OPT=opt ZLIB=false', then also run 'make COMP=gnu OPT=opt ZLIB=false'
in the PolySCIP directory

After a successful build, the executable will be located in the
directory 'bin/'

3)
a) Execute 'make doc' to build the doxygen documentation in 'doc/html'
b) Execute 'cd doc; pdflatex userguide.tex' to compile the user guide





