If you want to compile SCIP yourself, it is recommended to download the SCIP Optimization Suite tarball from [scipopt.org](https://scipopt.org/index.php#download) as it contains also the LP solver [SoPlex](https://soplex.zib.de) and the presolver [PaPILO](https://github.com/scipopt/papilo).

We provide two different systems to compile the code:
- the [CMake](https://scipopt.org/doc/html/md_INSTALL.php#CMAKE) build system (recommended for new users) [(local link)](@ref CMAKE).
- the traditional [Makefile](https://scipopt.org/doc/html/md_INSTALL.php#MAKE) system [(local link)](@ref MAKE).

Be aware that generated libraries and binaries of both systems might be different and incompatible.
For further information please refer to the [online documentation of SCIP](https://scipopt.org/doc/html/INSTALL.php).

The easiest way to install SCIP is to use the SCIP Optimization Suite, which contains SCIP, SoPlex, and ZIMPL.
For that we refer to the `README.md` file of the SCIP Optimization Suite (in case of the SCIP Optimization Suite there is no need to specify any directories, the compiling process is fully automated).

Building SCIP using CMake {#CMAKE}
==================================

[CMake](https://cmake.org/) is a build system generator that can create, e.g., Makefiles for UNIX and Mac or Visual Studio project files for Windows.

CMake provides an [extensive documentation](https://cmake.org/cmake/help/latest/manual/cmake.1.html) explaining available features and use cases as well as an [FAQ section](https://gitlab.kitware.com/cmake/community/-/wikis/FAQ).
It's recommended to use the latest stable CMake version available.  `cmake --help` is also a good first step to see available options and usage information.

Windows and platform independent build instructions
---------------------------------------------------

To build SCIP you may use the CMake GUI to specify the path to SCIP and the desired location for the build.
Available options are listed and can be modified to fit your needs.
After the configuration step is done, open the generated Visual Studio solution file and compile it.
Note that compilation is tested on MSVC version >= 12.

Alternatively, you may use the command line to configure and build SCIP by creating a `build` directory and then building the configuration:

```
cmake -Bbuild -H. [-DSOPLEX_DIR=/path/to/soplex]
cmake --build build --config Release [Debug]
```

Command line instructions (Linux, macOS)
----------------------------------------

Compiling SCIP directly can be done as follows:

```
tar xvzf scip-x.y.z.tgz                                                       # unpack the tarball
cd scip-x.y.z                                                                 # change into the directory
mkdir build                                                                   # create a new directory
cd build                                                                      # change directories
cmake .. -DCMAKE_INSTALL_PREFIX=<install/dir> [-DSOPLEX_DIR=/path/to/soplex]  # configure the build
make                                                                          # start compiling SCIP
make check                                                                    # (recommended) check build
make install                                                                  # (optional) install SCIP executable, library, and headers
```

Note: For a full ctest run `ctest` instead of `make check` after compilation.

CMake checks for available third-party libraries like GMP or ZLIB and sets up the configuration accordingly.
Note that the symmetry codes [Bliss](https://users.aalto.fi/~tjunttil/bliss/) and Sassy (github.com/markusa4/sassy) are shipped with SCIP.

Note: Here is a list of apt package requirements for ubuntu or debian users that want to build the entire SCIP Optimization Suite from source tarball:
```
apt-get install wget cmake g++ m4 xz-utils libgmp-dev unzip zlib1g-dev libboost-program-options-dev libboost-serialization-dev libboost-regex-dev libboost-iostreams-dev libtbb-dev libreadline-dev pkg-config git liblapack-dev libgsl-dev flex bison libcliquer-dev gfortran file dpkg-dev libopenblas-dev rpm
```
Additionally the following dependencies need to be downloaded, compiled and installed:
 - [Hmetis](http://glaros.dtc.umn.edu/gkhome/metis/hmetis/download)
 - [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/download)
 - [Ipopt](https://github.com/coin-or/Ipopt/releases) with [Mumps](https://github.com/coin-or-tools/ThirdParty-Mumps/releases)
 - [Gmp](https://gmplib.org/#DOWNLOAD)
During the CMake configuration of the SCIP Optimization Suite the can be specified, see [CMake](https://scipopt.org/doc/html/md_INSTALL.php#CMAKE) [(local link)](@ref CMAKE) .


Troubleshooting
---------------

If you have a problem with your cmake configuration and just want to build scip with all the available dependencies, the **simplest solution is to activate the `AUTOBUILD` option**:
```
cmake .. -DAUTOBUILD=on
```
This option activates the automatic search for dependent packages like GMP, IPOPT, PaPILO, Readline, WORHP, ZIMPL, ZLIB, and deactivates the missing ones.

If you need a specific package that is not automatically found, you should try setting a hint to the installation with the specified variable.
Specific packages can also be disabled individually.

**Examples of errors and possible solutions:**

Problem:
```
-- Finding PAPILO
-- Could NOT find PAPILO (missing: PAPILO_DIR)
CMake Error at CMakeLists.txt:359 (message):
  PAPILO not found, try specifying PAPILO_DIR.

  If you have troubles configuring, you can consider setting AUTOBUILD=ON to
  try and find optional packages as available.


-- Configuring incomplete, errors occurred!
```
Solution: add `-DPAPILO_DIR=/path/to/papilo/installation` or disable PaPILO by setting `-DPAPILO=off`.

Problem:
```
-- Finding ZIMPL
-- Could NOT find ZIMPL (missing: ZIMPL_DIR)
CMake Error at CMakeLists.txt:533 (message):
  ZIMPL not found, try specifying ZIMPL_DIR.

  If you have troubles configuring, you can consider setting AUTOBUILD=ON to
  try and find optional packages as available.


-- Configuring incomplete, errors occurred!
```
Solution: add `-DZIMPL_DIR=/path/to/zimpl/installation` or disable ZIMPL by setting `-DZIMPL=off`.

Problem:
```
-- Finding IPOPT
-- Could NOT find IPOPT (missing: IPOPT_LIBRARIES) (Required is at least version "3.12.0")
CMake Error at CMakeLists.txt:564 (message):
  IPOPT not found, try specifying IPOPT_DIR.

  If you have troubles configuring, you can consider setting AUTOBUILD=ON to
  try and find optional packages as available.

-- Configuring incomplete, errors occurred!
```
Solution: add `-DIPOPT_DIR=/path/to/ipopt/installation` or disable IPOPT by setting `-DIPOPT=off`.

Problem:
```
-- Finding Solver "spx"
-- Finding Soplex
CMake Error at CMakeLists.txt:375 (find_package):
  Could not find a package configuration file provided by "SOPLEX" with any
  of the following names:

    SOPLEXConfig.cmake
    soplex-config.cmake

  Add the installation prefix of "SOPLEX" to CMAKE_PREFIX_PATH or set
  "SOPLEX_DIR" to a directory containing one of the above files.  If "SOPLEX"
  provides a separate development package or SDK, be sure it has been
  installed.


-- Configuring incomplete, errors occurred!
```
Solution: add `-DSOPLEX_DIR=/path/to/soplex/installation` or disable SOPLEX by setting `-DLPS=none` or select a different lp solver you have available by `-DLPS=grb -DGUROBI_DIR=/path/to/gurobi/installation` or `-DLPS=xprs -DXPRESS_DIR=/path/to/xpress/installation` or `-DLPS=msk -DMOSEK_DIR=/path/to/mosek/installation` or `-DLPS=cpx -DCPLEX_DIR=/path/to/cplex/installation` or `-DLPS=glob -DGLOB_DIR=/path/to/glob/installation`

Modifying a CMake configuration
-------------------------------

CMake uses an out-of-source build, i.e., compiled binaries and object files are separated from the source tree and located in another directory.
Usually this directory is called `build` or `debug` or whatever you prefer.
From within this directory, run `cmake <path/to/SCIP>` to configure your build, followed by `make` to compile the code according to the current configuration (this assumes that you chose Linux Makefiles as CMake Generator).
By default, SCIP searches for SoPlex as LP solver.
If SoPlex is not installed systemwide, the path to a CMake build directory of SoPlex must be specified (i.e. one that contains "soplex-config.cmake").
Alternatively, a different LP solver can be specified with the `LPS` variable.

Afterwards, successive calls to `make` are going to recompile modified source code, without requiring another call to `cmake`.
The initial configuration step checks your environment for available third-party libraries and packages and sets up the configuration accordingly, e.g., disabling support for GMP if not installed.

The generated executable and libraries are put in directories `bin` and `lib` respectively and will simply be named `scip` or `libscip.so`.
This is different from the naming convention of the Makefile setup that appended the configuration details like OS and third party dependencies directly to the name of the binary or library.
The CMake setup tries to follow the established Linux/UNIX compilation conventions to facilitate the use of the libraries in other applications.
The previously generated sub-libraries like `liblpi.so` or `libobjscip.so` are not created by default anymore.
They can be built using the respective targets `liblpi`, `libobjscip`, etc.
The main library `libscip.so` will contain all SCIP sources and won't have dependencies to the other sub-libraries.

There are several options that can be passed to the `cmake <path/to/SCIP>` call to modify how the code is built.
For all of these options and parameters you have to use `-D<Parameter_name>=<value>`.
Following a list of available options, for the full list run

```
cmake </path/to/SCIP/> -LH
```
and set them by running `cmake .. -D<option>=<value>`.
Options can be chained together or be specified in subsequent calls to cmake.
The existing configuration will be updated or extended.
e.g., `cmake </path/to/SCIP> -DSOPLEX_DIR=<path/to/SoPlex/build/or/install>`.

| CMake option           | Available values                   | Makefile equivalent        | Remarks                                                            |
|------------------------|------------------------------------|----------------------------|--------------------------------------------------------------------|
| `AUTOBUILD`            | `on`, `off`                        | --                         | automatically find dependencies on availability, ignores individual flags of these packages |
| `CMAKE_BUILD_TYPE`     | `Release`, `Debug`, ...            | `OPT=[opt, dbg]`           |                                                                    |
| `GMP`                  | `on`, `off`                        | `GMP=[true, false]`        | specify `GMP_DIR` if not found automatically                       |
| `IPOPT`                | `on`, `off`                        | `IPOPT=[true,false]`       | requires IPOPT version >= 3.12.0; specify `IPOPT_DIR` if not found automatically |
| `LPS`                  | `spx`, `cpx`, `grb`, `xprs`, ...   | `LPS=...`                  | specify `SOPLEX_DIR`, `CPLEX_DIR`, `MOSEK_DIR`, ... if LP solver is not found automatically |
| `SYM`                  | `bliss`, `sassy`, `none`           | `SYM=[bliss, sassy, none]` | choose symmetry handling                                           |
| `WORHP`                | `on`, `off`                        | `WORHP=[true,false]`       | should worhp be linked; specify `WORHP_DIR` if not found automatically |
| `ZIMPL`                | `on`, `off`                        | `ZIMPL=[true, false]`      | specify `ZIMPL_DIR` if not found automatically                     |
| `AMPL`                 | `on`, `off`                        | `AMPL=[true, false]`       |                                                                    |
| `READLINE`             | `on`, `off`                        | `READLINE=[true, false]`   |                                                                    |
| `..._DIR`              | `<custom/path/to/.../package>`     | --                         | e.g. `IPOPT_DIR`, `CPLEX_DIR`, `WORHP_DIR`, `Readline_DIR` ...     |
| `BOOST_ROOT`           | `<custom/path/to/.../boost>`       | --                         | hint for location of boost                                         |
| `CMAKE_INSTALL_PREFIX` | `\<path\>`                         | `INSTALLDIR=\<path\>`      |                                                                    |
| `SHARED`               | `on`, `off`                        | `SHARED=[true, false]`     |                                                                    |
| `CXXONLY`              | `on`, `off`                        | --                         | use a C++ compiler for all source files                            |
| `COVERAGE`             | `on`, `off`                        | --                         | use with gcc, lcov, gcov in **debug** mode                         |
| `COVERAGE_CTEST_ARGS`  | ctest argument string              | --                         | see `ctest --help` for arguments                                   |
| `DEBUGSOL`             | `on`, `off`                        | `DEBUGSOL=[true,false]`    | specify a debugging solution by setting the "misc/debugsol" parameter of SCIP |
| `LPSCHECK`             | `on`, `off`                        | `LPSCHECK=[true,false]`    | double check SoPlex results with CPLEX                             |
| `NOBLKMEM`             | `on`, `off`                        | `NOBLKMEM=[true,false]`    |                                                                    |
| `NOBUFMEM`             | `on`, `off`                        | `NOBUFMEM=[true,false]`    |                                                                    |
| `NOBLKBUFMEM`          | `on`, `off`                        | `NOBLKBUFMEM=[true,false]` |                                                                    |
| `MT`                   | `on`, `off`                        |                            | use static runtime libraries for Visual Studio compiler on Windows |
| `THREADSAFE`           | `on`, `off`                        | `THREADSAFE=[true,false]`  | thread safe compilation                                            |
| `SANITIZE_...`         | `on`, `off`                        | --                         | enable sanitizer in debug mode if available                        |
| `TPI`                  | `tny`, `omp`, `none`               | `TPI=[tny,omp,none]`       | enable task processing interface required for concurrent solver    |

Parameters can be set all at once or in subsequent calls to `cmake` - extending or modifying the existing configuration.

Testing with CTest
------------------

There is an extensive test suite written for [CTest]("https://cmake.org/cmake/help/latest/manual/ctest.1.html)
that may take a while to complete.
To perform a quick test to see whether the compilation was really successful you may run `make check`.
To see all available tests, run

```
ctest -N
```

and to perform a memory check, run

```
ctest -T MemCheck
```

If [Criterion](https://criterion.readthedocs.io/en/master/) is installed (set custom path with `-DCRITERION_DIR=<path>`) the target `unittests` can be used to compile and run the available unit tests.

A coverage report for the entire test suite can be generated.
This requires a modification of the compilation process.
Two variables govern the report generation, `COVERAGE` and `COVERAGE_CTEST_ARGS`.
It is recommended to use the Debug build type.

```
cmake </path/to/SCIP> -DCOVERAGE=on -DCOVERAGE_CTEST_ARGS="-R MIP -E stein -j4" -DCMAKE_BUILD_TYPE=Debug
```

In this example, coverage is enabled in combination with the build type Debug.
In addition, only the coverage for tests with "MIP" in the name are run, excluding those that have "stein" in the name.
The tests are performed in parallel using 4 cores.

Use the `coverage` target, e.g., `make coverage`, to build the coverage report.
The generated report can be found under "coverage/index.html".

Additional targets
------------------

There are several further targets available, which can be listed using `make help`.
For instance, there are some examples that can be built with `make examples` or by specifying a certain one: `make <example-name>`.
For detailed instructions see the [installation instructions for applications and examples](https://scipopt.org/doc/html/INSTALL_APPLICATIONS_EXAMPLES.php) [(local link)](@ref INSTALL_APPLICATIONS_EXAMPLES).

| CMake target    | Description                                           | Requirements                          |
|-----------------|-------------------------------------------------------|---------------------------------------|
| scip            | build SCIP executable                                 |                                       |
| applications    | build executables for all applications                |                                       |
| examples        | build executables for all examples                    |                                       |
| unittests       | build unit tests                                      | the Criterion package                 |
| all_executables | build all of the above                                |                                       |
| libscip         | build the SCIP library                                |                                       |
| install         | install SCIP                                          |                                       |
| coverage        | run the test suite and create a coverage report       | build flag `-DCOVERAGE=on`            |
| liblpi          | build the LPI library                                 |                                       |
| libobjscip      | build the ObjSCIP library for the C++ wrapper classes |                                       |


Building SCIP using the Makefile system {#MAKE}
===============================================

For Linux and Mac, reading the section "Brief installation description" below should usually be enough.
If this is not the case, you can find the "Detailed installation description" below as well as some examples.

We recommend using GCC version 4.8 or later.

Brief installation description
------------------------------

The easiest way to install SCIP is to use the SCIP Optimization Suite which contains SCIP, SoPlex, and ZIMPL.
For that we refer to the `README.md` file of the SCIP Optimization Suite (main advantage: there is no need
to specify any directories, the compiling process is fully automated).

Compiling SCIP directly can be done as follows:

```
tar xvzf scip-x.y.z.tgz                      # unpack the tarball
cd scip-x.y.z                                # change into the directory
make                                         # start compiling SCIP
make test                                    # (recommended) check your SCIP installation
make install INSTALLDIR=/path/to/install/dir # (optional) install the header, libraries, and binary
```

On your first compilation you will be asked for some soft-link targets, depending on the LP solver you want to use.
Usually, SCIP needs the following information
   a. the directory where the include files of the LP solver are located
   b. the library file(s) `lib*.a` or/and `lib*.so`

Besides that, SCIP needs similar soft-link targets for ZIMPL
   a. the directory where the include files of ZIMPL are located
   b. the library file(s) `lib*.a` or/and `lib*.so`

You will need either the `.a` or the `.so` files and can skip the others by just pressing return.

The most common compiling issue is that some libraries are missing on your system or that they are outdated.
SCIP by default requires the following packages (with usual names for Linux systems in parentheses):

- zlib     (libz-dev)
- gmp      (libgmp-dev),
- readline (libreadline-dev), and
- ncurses  (libncurses-dev)

Note that under Linux-based systems, you need to install the developer-versions of gmp/zlib/readline, in order to also have the header-files available.

If you are not able or do not want to install these packages, try compiling with:
```
make ZLIB=false READLINE=false GMP=false.
```

Detailed installation description
---------------------------------

Here is what you have to do to get SCIP running:

### 1. Compile the library and the solver program

In your SCIP main directory, enter `make [options]` with the following options:

| parameter and default | options              | description                                                                                      |
|-----------------------|----------------------|--------------------------------------------------------------------------------------------------|
| `ARCH=x86_64`         | `[x86_64, x86, sparc, mips, hppa, ppc, pwr4]` | the architecture: try to autodetect                      |
| `AMPL=true`           | `[true, false]`      | to enable or disable AMPL .nl file reader and support for using SCIP executable as solver in AMPL|
| `COMP=gnu`            | `[gnu, clang, intel]`| Use Gnu, Clang or Intel compiler.                                                                |
| `EXPRINT=cppad`       | `[cppad, none]`      | to use CppAD as expressions interpreter                                                          |
| `FILTERSQP=false`     | `[false, true]`      | to enable or disable FilterSQP interface                                                         |
| `GMP=true`            | `[true, false]`      | to enable or disable GMP library for exact counting and Zimpl support                            |
| `IPOPT=false`         | `[false, true]`      | to disable or enable IPOPT interface (needs IPOPT >= 3.12.0)                                     |
| `LPS=spx`             | `[spx1, cpx, grb, xprs, msk, clp, glop, qso, none]` | determines the LP-Solver, should be installed seperately. Options to use SoPlex (> version 2.0), SoPlex (>= version 1.4), CPLEX, Gurobi, XPRESS, MOSEK, CLP, Glop, QSopt as LP solver, no LP solver  |
| `LPSOPT=opt`          | `[opt, dbg, opt-gccold]` | Choose the debug or optimized version (or old GCC optimized) version of the LP-solver (currently only available for SoPlex and CLP). |
| `NOBLKMEM=false`      | `[false, true]`      | Turns the internal SCIP block memory off or on.                                                  |
| `NOBUFMEM=false`      | `[false, true]`      | Turns the internal SCIP buffer memory off or on.                                                 |
| `NOBLKBUFMEM=false`   | `[false, true]`      | Turns the internal SCIP block and buffer memory off or on. This way the code can be checked by valgrind or similar tools. |
| `OPT=opt`             | `[opt, dbg, perf]`   | to use optimized, debug, performance (only with Gnu compiler) analysis compilation mode. `dbg` turns on debug mode. This enables asserts and avoids macros for several function in order to ease debugging. |
| `OSTYPE`              | `[linux, darwin, cygwin, irix, windows, mingw]` | the operating system: try to autedetect                               |
| `PAPILO=false`        | `[false, true]`      | to disable or disable the MILP presolver based on the presolving library PaPILO                  |
| `READLINE=true`       | `[true, false]`      | to enable or disable readline library for interactive shell                                      |
| `SHARED=false`        | `[false, true]`      | to suppress or create shared libraries (only Gnu compiler)                                       |
| `SYM=none`            | `[none, bliss, sassy]` | to choose method for computing symmetries in mixed nonlinear integer programs                  |
| `TPI=none`            | `[none, omp, tny]`   | to disable the task processing interface or use it with the openmp or tinycthreads interface for concurrent solves |
| `VERBOSE=false`       | `[false, true]`      | to suppress or display of compiler and linker invocations                                        |
| `WORHP=false`         | `[false, true]`      | to disable or enable WORHP interface (needs WORHP >= 2.00)                                       |
| `ZIMPL=false`         | `[false, true, auto]`| to enable or disable ZIMPL file reader (needs ZIMPL and GMP to be installed)                     |
| `ZLIB=true`           | `[true, false]`      | to enable or disable zlib for reading of compressed files                                        |

For example, if you want to install SCIP on a Linux system with a x86 processor
using the gnu compiler in debug mode, using Soplex version as LP solver,
and neither an expressions interpreter nor symmetry handling techniques or multi-threading,
you would use following:

```
make OSTYPE = linux  ARCH = x86  COMP = gnu  OPT = dbg  EXPRINT = none
```

On some machines, you should use `gmake` instead of `make`.

On your first compilation you will be asked for some soft-link targets, depending on the external software you want to use.
Usually, SCIP needs the following information:
- the directory where the include files of the external software
- the library file(s) `lib*.a` or/and `lib*.so`
You will need either the `.a` or the `.so` files and can skip the others by just pressing return.

On MAC systems, GMP is often not installed in the library and include paths, e.g. in `/sw/include` and `/sw/lib`.
In this case, you have to add the paths explicitly.
In the above example add the settings:

```
USRFLAGS=-I/sw/include USRCPPFLAGS=-I/sw/include USRCFLAGS=-I/sw/include USRLDFLAGS=-L/sw/lib.
```

### 2. Installing SCIP

After compiling you can install the headers, the libraries, and the binary.
You do that by running the command:

```
make install INSTALLDIR=<directory>
```
where you substitute

- `INSTALLDIR=` to install in current directory (default)
- `INSTALLDIR=/usr/local` to install the headers (`/usr/local/include/`), the libraries (`/usr/local/lib/`), and binary (`/usr/local/bin/`) in the directory `/usr/local`

For un-installing SCIP there exist the target `uninstall` which can be used in the same way as `install`.

### 3. Instructions for manually creating the soft-links, if the query script fails:

Create necessary soft-links in the `lib/static` and `lib/include/` subdirectories of SCIP:

#### a) to use SOPLEX (Version >= 1.4.0)

For each used operating system and architecture:
```
ln -s <path to SOPLEX' *.h files> <path to SCIP>/lib/include/spxinc
ln -s <file libsoplex.[...].a> <path to SCIP>/lib/static/libsoplex.$(OSTYPE).$(ARCH).$(COMP).a
```
For example:
```
cd scip
ln -s /soplex/lib/libsoplex.linux.x86_64.gnu.opt.a lib/static/libsoplex.linux.x86_64.gnu.a
```
Warning! The `.opt` in the name of the SOPLEX library does not appear in the name of the soft-link.

#### b) to use CPLEX (Version >= 10.0)

For each used operation system and architecture:
```
ln -s <path to directory of cplex.h> <path to SCIP>/lib/include/cpxinc
ln -s <file libcplex.a> <path to SCIP>/lib/static/libcplex.$(OSTYPE).$(ARCH).$(COMP).a
```
For example:
```
cd scip
ln -s /cplex121/include/ilcplex lib/include/cpxinc
ln -s /cplex121/lib/x86-64_debian4.0_4.1/static_pic/libcplex.a lib/static/libcplex.linux.x86.gnu.a
```

#### c) to use Gurobi

For each used operation system and architecture:
```
ln -s <path to the include directory of Gurobi> <path to SCIP>/lib/include/grbinc
ln -s <file libgurobi81.so> <path to SCIP>/lib/shared/libgurobi.$(OSTYPE).$(ARCH).$(COMP).so
```
For example:
```
cd scip
ln -s /gurobi81/linux64/include lib/include/grbinc
ln -s /gurobi81/linux64/lib/libgurobi81.so lib/shared/libgurobi.linux.x86_64.gnu.so
```

#### d) to use XPRESS

For each used operation system and architecture:
```
ln -s <path to directory of xprs.h> <path to SCIP>/lib/include/xprsinc
ln -s <file libxprs.a> <path to SCIP>/lib/static/libxprs.$(OSTYPE).$(ARCH).$(COMP).a
```
For example:
```
cd scip
ln -s /xpressmp/include lib/include/xprsinc
ln -s /xpressmp/lib/libxprs.a lib/static/libxprs.linux.x86.gnu.a
```

#### e) to use MOSEK

For each used operation system and architecture:
```
ln -s <path to directory of mosek.h> <path to SCIP>/lib/include/mskincn
ln -s <file libmosek.so> <path to SCIP>/lib/shared/libmosek.$(OSTYPE).$(ARCH).$(COMP).so
```
For example:
```
cd scip
ln -s /mosek/8/tools/platform/linux64x86/h lib/include/mskinc
ln -s /mosek/8/tools/platform/linux64x86/bin/libmosek64.so lib/shared/libmosek.linux.x86_64.gnu.so
```

#### f) to use CLP

For each used operating system and architecture:
```
ln -s <path to Clp main directory> <path to SCIP>/lib/include/libclp.$(OSTYPE).$(ARCH).$(COMP).$(LPSOPT)
```
For example:
```
cd scip
ln -s /Coin-Clp lib/include/libclp.linux.x86.gnu.opt
```

#### g) to use Glop

```
ln -s <path to OR-Tools main directory> <path to SCIP>/shared/ortools
```
For example:
```
cd scip
ln -s /ortools lib/shared/ortools
```

#### h) to use ZIMPL

Use ZIMPL as additional file reader for reading *.zpl files:
```
mkdir <path to SCIP>/lib/include/zimplinc
ln -s <path to ZIMPL's *.h files> <path to SCIP>/lib/include/zimplinc/zimpl
ln -s <file libzimpl-<version>.<options>.a> <path to SCIP>/lib/static/libzimpl.$(OSTYPE).$(ARCH).$(COMP).a
```
Note that ZIMPL needs the GNU multiprecision library (GMP) to be installed on your system.

#### i) to use IPOPT as NLP solver

```
ln -s <path to IPOPT installation> <path to SCIP>/lib/ipopt.$(OSTYPE).$(ARCH).$(COMP).$(IPOPTOPT)
(e.g. `cd scip; ln -s /Ipopt lib/shared/ipopt.linux.x86.gnu.opt
```
The path to the IPOPT installation is the path under where the Ipopt build has been installed.
It should contain the directories `include/coin-or` with the Ipopt header files, the directory `lib` with the Ipopt libraries, and the file `lib/pkgconfig/ipopt.pc`.

#### j) to use WORHP as NLP solver

```
ln -s <path to WORHP installation> <path to SCIP>/lib/shared/worhp.$(OSTYPE).$(ARCH).$(COMP).$(WORHPOPT)
```
For example:
```
cd scip
ln -s /Worhp lib/shared/worhp.linux.x86.gnu.opt
```
The path to the WORHP installation is the path under where the Worhp build has been installed.
It should contain the directories `include/worhp` with the WORHP header files and the directory `lib` with the WORHP libraries.

#### k) to use FilterSQP as NLP solver

```
ln -s <path to FilterSQP library> <path to SCIP>/lib/libfiltersqp.$(OSTYPE).$(ARCH).$(COMP).a
ln -s <path to BQPD library> <path to SCIP>/lib/libbqpd.$(OSTYPE).$(ARCH).$(COMP).a
```
Make sure to replace the paths with your installation location.

#### l) to use GAMS

```
ln -s <path to GAMS system directory> <path to SCIP>/lib/shared/gams.$(OSTYPE).$(ARCH).$(COMP)
```
Make sure to replace the paths with your installation location.

### 4. Run SCIP

To run SCIP enter `bin/scip.$(OSTYPE).$(ARCH).$(COMP).$(OPT).$(LPS)`
(e.g. `bin/scip.linux.x86.gnu.opt.spx`) or just `bin/scip` for the last compiled version

### 5. Generate documentation

To generate the documentation, you need to have doxygen installed, and enter `make doc`.

### 6. Check Code with (pc)lint

To check the code with lint, you need to have flexelint installed, and enter `make lint`.
If you have pclint installed, enter `make pclint`.

### 7. Run a short test

To run a short test, enter `make [options] test` with the same options with which you compiled SCIP in step 1.
If you use `EXPRINT=none`, a few MINLP instances might be aborted.
If you use `LPS=none`, many instances will fail or take ages to be solved.

Further targets
---------------
The SCIP makefile supports several targets (used via `make ... "target"`):

| target | description|
|--|--|
| `all`    | (or no target) Build SCIP library and binary.                            |
| `links`  | Reconfigures the links in the `lib` directory.                           |
| `doc`    | Creates documentation in the `doc` directory.                            |
| `clean`  | Removes all object files.                                                |
| `depend` | Updates dependencies files. This is only needed if you add checks for preprocessor-defines `WITH_*` in source files. |
| `check`  | or `test`. Runs the check script.                                         |
| `lint`   | Statically checks the code via flexelint. The call produces the file `lint.out` which contains all the detected warnings. |
| `tags`   | Generates tags which can be used in the editor **emacs** and **xemacs**. |

The SCIP makefiles are structured as follows.

- `Makefile` This is the basic makefile in the SCIP root directory. It loads
  additional makefile information depending on the parameters set.
- `make/make.project` This file contains definitions that are useful for all codes
  that use SCIP, for instance, the examples.
- `make.\<sys\>.\<machine\>.\<compiler\>.\<dbg|opt|prf|opt-gccold\>` These file contain system/compiler specific
  definitions. If you have an unsupported compiler, you can copy one of these and modify it
  accordingly.

If your platform or compiler is not supported by SCIP you might try and copy one of the existing
makefiles in the `make` directory and modify it. If you succeed, we are always
interested in including more Makefiles into the system.

Examples
--------

###  Example 1 (defaults: SoPlex, with ZIMPL support):

Typing `make` uses SoPlex as LP solver and includes support for the modeling language ZIMPL.
You will be asked the following questions on the first call to `make` (example answers are already given):

```
make[1]: Entering directory '/sw/scip'

** creating softlinks: LPS=spx OSTYPE=linux ARCH=x86 COMP=gnu SUFFIX= ZIMPL=true ZIMPLOPT=opt IPOPT=false IPOPTOPT=opt EXPRINT=cppad

** creating directory 'lib/zimplinc'
** missing soft-link 'lib/spxinc'
** enter soft-link target file or directory for 'lib/spxinc' (return if not needed): /sw/soplex/src
-> creating softlink 'lib/spxinc' -> '/sw/soplex/src'

** missing soft-link 'lib/libsoplex.linux.x86.gnu.a'
** enter soft-link target file or directory for 'lib/libsoplex.linux.x86.gnu.a' (return if not needed): /sw/soplex/lib/libsoplex.linux.x86.gnu.opt.a
-> creating softlink 'lib/libsoplex.linux.x86.gnu.a' -> '/sw/soplex/lib/libsoplex.linux.x86.gnu.opt.a'

** missing soft-link 'lib/libsoplex.linux.x86.gnu.so'
** this soft-link is not necessarily needed since 'lib/libsoplex.linux.x86.gnu.a' already exists - press return to skip
** enter soft-link target file or directory for 'lib/libsoplex.linux.x86.gnu.so' (return if not needed):
-> skipped creation of softlink 'lib/libsoplex.linux.x86.gnu.so'. Call 'make links' if needed later.

** missing soft-link 'lib/zimplinc/zimpl'
** enter soft-link target file or directory for 'lib/zimplinc/zimpl' (return if not needed): /sw/zimpl/src
-> creating softlink 'lib/zimplinc/zimpl' -> '/sw/zimpl/src'

** missing soft-link 'lib/libzimpl.linux.x86.gnu.a'
** enter soft-link target file or directory for 'lib/libzimpl.linux.x86.gnu.a' (return if not needed): /sw/zimpl/lib/libzimpl.linux.x86.gnu.opt.a
-> creating softlink 'lib/libzimpl.linux.x86.gnu.a' -> '/sw/zimpl/lib/libzimpl.linux.x86.gnu.opt.a'

** missing soft-link 'lib/libzimpl.linux.x86.gnu.so'
** this soft-link is not necessarily needed since 'lib/libzimpl.linux.x86.gnu.a' already exists - press return to skip
** enter soft-link target file or directory for 'lib/libzimpl.linux.x86.gnu.so' (return if not needed):
-> skipped creation of softlink 'lib/libzimpl.linux.x86.gnu.so'. Call 'make links' if needed later.

make[1]: Leaving directory '/sw/scip'
```

### Example 2 (CPLEX, no ZIMPL):


Typing `make LPS=cpx ZIMPL=false` uses CPLEX as LP solver.
You will be asked the following questions on the first call to `make` (example answers are already given):

```
make[1]: Entering directory '/sw/scip'

** creating softlinks: LPS=cpx OSTYPE=linux ARCH=x86 COMP=gnu SUFFIX= ZIMPL=false

** missing soft-link 'lib/cpxinc'
** enter soft-link target file or directory for 'lib/cpxinc' (return to skip): /sw/cplex/include/ilcplex
-> creating softlink 'lib/cpxinc' -> '/sw/cplex/include/ilcplex'

** missing soft-link 'lib/libcplex.linux.x86.gnu.a'
** enter soft-link target file or directory for 'lib/libcplex.linux.x86.gnu.a' (return to skip): /sw/cplex/lib/x86_rhel4.0_3.4/static_pic/libcplex.a
-> creating softlink 'lib/libcplex.linux.x86.gnu.a' -> '/sw/cplex/lib/x86_rhel4.0_3.4/static_pic/libcplex.a'

** missing soft-link 'lib/libcplex.linux.x86.gnu.so'
** enter soft-link target file or directory for 'lib/libcplex.linux.x86.gnu.so' (return to skip):
-> skipped creation of softlink 'lib/libcplex.linux.x86.gnu.so'. Call 'make links' if needed later.

make[1]: Leaving directory '/sw/scip'
```

###   Example 3 (CLP, IPOPT, no ZIMPL):

Typing `make LPS=clp ZIMPL=false IPOPT=true` uses CLP as LP solver, and activates the interface to IPOPT.
You will be asked the following questions on the first call to `make` (example answers are already given):

```
make[1]: Entering directory '/sw/scip'

- Current settings: LPS=clp OSTYPE=linux ARCH=x86_64 COMP=gnu SUFFIX= ZIMPL=false ZIMPLOPT=opt IPOPT=true IPOPTOPT=opt EXPRINT=cppad

* SCIP needs some softlinks to external programs, in particular, LP-solvers.
* Please insert the paths to the corresponding directories/libraries below.
* The links will be installed in the 'lib' directory.
* For more information and if you experience problems see the 'INSTALL.md' file.

  -> 'clp.*' is a directory containing the Clp installation, i.e., 'clp.*/include/coin/ClpModel.hpp' should exist.
  -> 'ipopt.*' is a directory containing the ipopt installation, i.e., 'ipopt.*/include/coin/IpIpoptApplication.hpp', 'ipopt.*/lib/libipopt*', ... should exist.

- preparing missing soft-link 'lib/clp.linux.x86_64.gnu.opt':
> Enter soft-link target file or directory for 'lib/clp.linux.x86_64.gnu.opt' (return if not needed):
> /sw/Clp-1.11/build
-> creating softlink 'lib/clp.linux.x86_64.gnu.opt' -> '/sw/Clp-1.11/build'

- preparing missing soft-link 'lib/ipopt.linux.x86_64.gnu.opt':
> Enter soft-link target file or directory for 'lib/ipopt.linux.x86_64.gnu.opt' (return if not needed):
> /sw/ia64_lx26/ipopt-3.12.5/
-> creating softlink 'lib/ipopt.linux.x86_64.gnu.opt' -> '/sw/ia64_lx26/ipopt-3.12.5/'

make[1]: Leaving directory '/sw/scip'
```

### Example 4 (default: SoPlex, IPOPT, WORHP, FILTERSQP):

Typing `make IPOPT=true WORHP=true FILTERSQP=true` uses SoPlex as LP solver, and activates the interfaces to IPOPT, WORHP, and FilterSQP.
You will be asked the following questions on the first call to `make` (example answers are already given):

```
- Current settings: LPS=spx2 OSTYPE=linux ARCH=x86_64 COMP=gnu SHARED=false SUFFIX= ZIMPL=false ZIMPLOPT=opt IPOPT=true IPOPTOPT=opt FILTERSQP=true EXPRINT=cppad GAMS=false

* SCIP needs some softlinks to external programs, in particular, LP-solvers.
* Please insert the paths to the corresponding directories/libraries below.
* The links will be installed in the 'lib/include' and 'lib/static' directories.
* For more information and if you experience problems see the 'INSTALL.md' file.

  -> 'spxinc' is the path to the SoPlex 'src' directory, e.g., '<SoPlex-path>/src'.
  -> 'libsoplex.*' is the path to the SoPlex library, e.g., '<SoPlex-path>/lib/libsoplex.linux.x86.gnu.opt.a'
  -> 'ipopt.linux.x86_64.gnu.opt' is a directory containing the ipopt installation, i.e., 'ipopt.linux.x86_64.gnu.opt/include/coin/IpIpoptApplication.hpp', 'ipopt.linux.x86_64.gnu.opt/lib/libipopt*', ... should exist.
  -> 'libfiltersqp.linux.x86_64.gnu.*' is the path to the filterSQP library.
  -> 'libbqpd.linux.x86_64.gnu.*' is the path to the BQPD library.
  -> 'worhp.linux.x86_64.gnu.opt' is a directory containing the WORHP installation, i.e., 'worhp.linux.x86_64.gnu.opt/include/worhp/worhp.h' should exist.


> Enter soft-link target file or directory for 'lib/static/ipopt.linux.x86_64.gnu.opt' (return if not needed):
> /sw/ipopt-3.12.5
-> creating softlink 'lib/static/ipopt.linux.x86_64.gnu.opt' -> '/sw/ipopt-3.12.5'

> Enter soft-link target file or directory for 'lib/static/libfiltersqp.linux.x86_64.gnu.a' (return if not needed):
> /sw/libfiltersqp.a
-> creating softlink 'lib/static/libfiltersqp.linux.x86_64.gnu.a' -> '/sw/libfiltersqp.a'

> Enter soft-link target file or directory for 'lib/static/libbqpd.linux.x86_64.gnu.a' (return if not needed):
> /sw/libbqpd.a
-> creating softlink 'lib/static/libbqpd.linux.x86_64.gnu.a' -> '/sw/libbqpd.a'

> Enter soft-link target file or directory for 'lib/static/worhp.linux.x86_64.gnu.opt' (return if not needed):
> /sw/worhp-2.0
-> creating softlink 'lib/static/worhp.linux.x86_64.gnu.opt' -> '/sw/worhp-2.0'

make[1]: Leaving directory '/sw/scip'
```

Note on how to (locally) install CLP:
- create a target directory for the installation, e.g. `clp-build` (this is the directory SCIP has to link to)
- from within `clp-build`, run the `configure` script of coin-Clp, followed by `make install`

If you ever need to modify the soft-link targets, delete the soft-links in the `lib/` subdirectory and enter `make links` to generate them again.

After the soft-links have been created, the compilation of the source files should start.


Compilation problems
--------------------

If the soft-link query script does not work on your machine, read step 2 for instructions on manually creating the soft-links.

### No rule to make target lib/???
If you get an error message of the type
```
make: *** No rule to make target 'lib/???', needed by 'obj/O.linux.x86.gnu.opt/lib/scip/???.o'.  Stop.
```
the corresponding soft-link was not created or points to a wrong location.
Check the soft-link targets in the `lib/include`, `lib/static`, `lib/shared` subdirectories.
Try to delete all soft-links from those directories and call `make links` to generate them again.
If this still fails, read step 2 for instructions on manually creating the soft-links.

### No rule to make target make/make

If you get an error message of the type
```
make: *** No rule to make target 'make/make.?.?.?.?.?'.  Stop.
```
the corresponding machine dependent makefile for your architecture and compiler is missing.
Create one of the given name in the `make/` subdirectory.
You may take `make/make.linux.x86.gnu.opt` or any other file in the make subdirectory as example.

### No support for remove_history call

The readline library seems to differ slightly on different OS distributions.
Some versions do not support the `remove_history()` call.
In this case, you have to either add `-DNO_REMOVE_HISTORY` to the FLAGS in the appropriate `make/make.*` file, or to compile with `make USRFLAGS=-DNO_REMOVE_HISTORY`.
Make sure, the file `src/scip/dialog.c` is recompiled.
If this doesn't work either, disable the readline library with `make READLINE=false`.

### No support for sigaction method

On some systems, the `sigaction()` method is not available.
In this case, you have to either add `-DNO_SIGACTION` to the FLAGS in the appropriate `make/make.*` file, or to compile with `make USRFLAGS=-DNO_SIGACTION`.
Make sure, the file `src/scip/interrupt.c` is recompiled.

### No support for rand_r method

On some systems, the `rand_r()` method is not available.
In this case, you have to either add `-DNO_RAND_R` to the FLAGS in the appropriate `make/make.*` file, or to compile with `make USRFLAGS=-DNO_RAND_R`.
Make sure, the file `src/scip/misc.c` is recompiled.

### No support for strtok_r method

On some systems, the `strtok_r()` method is not available.
In this case, you have to either add `-DNO_STRTOK_R` to the FLAGS in the appropriate `make/make.*` file, or to compile with `make USRFLAGS=-DNO_STRTOK_R`.
Make sure, the file `src/scip/misc.c` is recompiled.

### No support for strerror_r method

On some systems, the `strerror_r()` method is not available.
In this case, you have to either add `-DNO_STRERROR_R` to the FLAGS in the appropriate `make/make.*` file, or to compile with `make USRFLAGS=-DNO_STRERROR_R`.
Make sure, the file `src/scip/misc.c` is recompiled.

### No support for read command

On some systems, the option [-e] is not available for the read command.
You have to compile with `READ=read`.

### Problems with Clp

In some situations, it may be necessary to adjust the flags for linking against Clp.
SCIP's makefile tries to find the file `clp_addlibs.txt`, which specifies the needed libraries.
The first thing you should check is whether `clp_addlibs.txt` is present at in path `libclp.*/share/coin/doc/Clp/` (you may have to correct this path for some Clp versions).
If this file is not present in your Clp version, SCIP tries to guess the paths and libraries: it assumes that Blas and Lapack are installed as system libraries (`libblas.a`, `liblapack.a`) and are not build into the CoinUtils library.
If that is different in your build of Clp, you may have to remove `$(LINKCXX_l)lapack$(LINKLIBSUFFIX)` from the `LPSLDFLAGS` in `Makefile` or `make.project`.
Also removing `$(LINKCXX_l)bz2$(LINKLIBSUFFIX)` may help in some cases.

### Compiler or linker errors

If you encounter other compiler or linker errors, you should recompile with `make VERBOSE=true ...` in order to get the full compiler invocation.
This might help to fix the corresponding machine dependent makefile in the make subdirectory.

Remarks on Installing under Windows using MinGW
-----------------------------------------------

To build your own Windows binaries under Windows, we recommend using the MinGW-Compiler with MSYS from mingw.org.

First install MSYS, then MinGW to the mingw folder inside the msys folder.
Now you need to install the following packages to the mingw folder:
- zlib (or use `ZLIB=false`)
- pcre (or use `ZIMPL=false` since pcre is needed for ZIMPL and ZIMPL-support in SCIP)
- gmplib (or use `ZIMPL=false` since gmplib is needed for ZIMPL and ZIMPL-support in SCIP)

(After calling `make clean` in the ZIMPL folder you will also need flex and bison to remake ZIMPL.
We recommend NOT to use `make clean` inside the ZIMPL-folder if you do not have these packages installed.)

You can download these additional packages as precompiled binaries:
- [zlib&pcre](http://gnuwin32.sourceforge.net/packages.html)
- [gmplib](http://cs.nyu.edu/exact/core/gmp/)
or compile the source on your own from the project homepages:
- [zlib](http://www.zlib.net/)
- [pcre](http://www.pcre.org/)
- [gmplib](http://www.gmplib.org/)
(The command `./configure --prefix=/mingw ; make ; make install` should succeed without problems and installs the packages into the mingw folder.)

Now `make READLINE=false` should be compiling without errors.
Please note that we do NOT support creating the doxygen documentation or readline-usage under Windows.

Since there are no real symlinks in MSYS, the include and library files of SoPlex and ZIMPL are actually copied into the SCIP-lib-folder.
When you recompile ZIMPL or SoPlex after compiling SCIP you have to copy the libraries manually into the SCIP-lib-folder and recompile SCIP afterwards.
