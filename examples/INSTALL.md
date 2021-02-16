Building and installing the examples      {#INSTALL_EXAMPLES}
====================================

Each example may be built in two ways, either by using Cmake
or by using SCIP's own Makefile system. The choice of the
build system depends on the desired target platform and the availability
of those tools there, as well as your personal preferences.
For most users, we recommend to use the CMake system, especially on
non Unix platforms.


Installation information using CMake
------------------------------------

Please compile SCIP first,
see the cmake section of INSTALL in the main SCIP directory for instructions,
or refer to the online documentation of SCIP.

In the following description, `scip_example_binary` inside its respective `SCIPExampleDir`
can be one of the following:
- `binpacking` from example `Binpacking`
- `brachistochrome`, `circle`, `circlepacking`, `gastrans`, `string` from example `CallableLibrary`
- `eventhdlr` from example `Eventhdlr`
- `gmi` from example `GMI`
- `lop` from example `LOP`
- `scipmip` from example `MIPSolver`
- `queens` from example `Queens`
- `relaxator` from example `Relaxator`
- `scflp` from example `SCFLP`
- `sciptsp` from example `TSP`
- `sudoku` from example `Sudoku`
- `vrp` from example `VRP`

The example can be compiled within the same build directory
as SCIP. Assuming that the build directory of SCIP was named "build",
simply execute

```
cmake --build build --target <scip_example_binary>
```

To build all examples at once, use

```
cmake --build build --target examples
```


It is also possible to build `scip_example_binary` in a stand-alone
build directory. Therefore, it is necessary to create the
stand-alone build directory first and configure the build using
CMake. This approach requires a systemwide installation of SCIP.
If SCIP is not installed systemwide, but in a local directory "/path/to/scip/installation",
this needs to be communicated as follows, by either specifying the `SCIP_DIR` variable or
adjusting the `CMAKE_PREFIX_PATH` variable.
The following commands need to be issued from the root directory of the example that should be built.

```
cmake -Bbuild -H. [-DSCIP_DIR=/path/to/scip/installation/lib/cmake/scip] [-DCMAKE_PREFIX_PATH=/path/to/scip/installation]
cmake --build build
```

If you are unsure what an installation directory is, "/path/to/scip/installation" should contain the directories "include" and "lib"
or the equivalents on your target operating systems.
If SCIP has been compiled into a build-directory as opposed to an installation directory, it is possible to point either of the two variables
`SCIP_DIR` or `CMAKE_PREFIX_PATH` to this build directory.
Finally, this specification should be used to give a local installation precedence over a systemwide installation of SCIP.

Please refer to the [online documentation of SCIP](http://scipopt.org/doc/html/CMAKE.php)
for a list of available configuration options and available tests.


Installation information for SCIP's custom Makefile system on Linux
-------------------------------------------------------------------

In the following, some of the names depend on your machine and your
compilation settings:

- `OSTYPE` : the operating system
             the string returned by `uname -s` in lower case with the following
             replacements:
             - "cygwin*" is replaced by only "cygwin"
             - "irix??" is replaced by only "irix"
             - "windows*" is replaced by only "windows"

- `ARCH`:   the architecture
             the string returned by `uname -m`, modified by the following
             rules to subsume some architectures:
              - "sun??" is replaced by "sparc"
              - "i?86" is replaced by "x86"
              - "IP??" is replaced by "mips"
              - "9000????" is replaced by "hppa"
              - "Power Macintosh" is replaced by "ppc"
              - "00??????????" is replaced by "pwr4"

- `COMP`:   the compiler
             `gnu`, `intel`, `compaq`, `sun`, `insure`, ... (see make/ directory)

- `OPT`:    the optimization level of compilation
             `dbg`, `opt`, or `prf`

- `LPS`:    the LP solver to use
             `spx`, `spx132`, `clp`, `cpx`, `xprs`, `msk`

Let `scip_example_binary` be one of the following inside its respective `SCIPExampleDir`:
- `binpacking` from example `Binpacking`
- `brachistochrome`, `circle`, `circlepacking`, `gastrans`, `string` from example `CallableLibrary`
- `scip` from example `Eventhdlr`
- `gmi` from example `GMI`
- `lop` from example `LOP`
- `scipmip` from example `MIPSolver`
- `queens` from example `Queens`
- `scip` from example `Relaxator`
- `scflp` from example `SCFLP`
- `sciptsp` from example `TSP`
- `sudoku` from example `Sudoku`
- `vrp` from example `VRP`

For example, if you want to install SCIP on a Linux system with a x86 processor
using the gnu compiler in debug mode, and using Soplex version >= 1.4.0
as LP solver, you would have the following names:
```
make OSTYPE=linux ARCH=x86 COMP=gnu OPT=dbg LPS=spx
```

Here is what you have to do to compile and run the example project using SCIP as a library:

1. Install and compile SCIP as described in the make section of the INSTALL file of SCIP's main
   directory, and make sure to create the necessary softlinks in SCIP's lib
   directory

2. In the directory SCIPExampleDir edit the variable SCIPDIR if necessary - it should
   point to the directory that contains SCIP.

3. Compile the scip_example_binary example project:
   In the main directory SCIPExampleDir, enter `make OPT=<...> LPS=<...> COMP=<...>`
   with the following options:
   - `OPT=opt`       to use optimized compilation mode (default)
   - `OPT=dbg`       to use debug compilation mode
   - `OPT=prf`       to use performance analysis compilation mode
   - `LPS=spx`       to use SOPLEX Version >= 1.4.0 as LP solver (default)
   - `LPS=cpx`       to use CPLEX as LP solver
   - `LPS=xprs`      to use XPRESS as LP solver
   - `LPS=msk`       to use MOSEK as LP solver
   - `LPS=clp`       to use CLP as LP solver
   - `COMP=gnu`      to use GNU c/c++ compiler (default)
   - other compilers are available (see make/ directory)

  For CallableLibrary:
   - `IPOPT=true`    to enable using Ipopt as NLP solver

4. To run the program enter `bin/scip_example_binary.$(OSTYPE).$(ARCH).$(COMP).$(OPT).$(LPS)`
   (e.g. `bin/scip_example_binary.linux.x86.gnu.opt.spx`) or `bin/scip_example_binary` which is a link
   to last compiled version

5. For the examples Binpacking, CallableLibrary, GMI, LOP, Queens, Relaxator, SCFLP, TSP, VRP
   you can generate a documentation. Make sure you have doxygen installed, and enter `make doc`

On some machines, you should use gmake instead of make.
