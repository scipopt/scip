# SCIP: Solving Constraint Integer Programs

Welcome to what is currently one of the fastest academically developed solvers
for mixed integer programming (MIP) and mixed integer nonlinear programming
(MINLP). In addition, SCIP provides a highly flexible framework for constraint
integer programming and branch-cut-and-price. It allows for total control of the
solution process and the access of detailed information down to the guts of the
solver.

The original instance of this repository is hosted at
[git.zib.de](https://git.zib.de) and a read-only
mirror is available at
[github.com/scipopt/scip](https://github.com/scipopt/scip).

Further information and resources are available through the official website at
[scipopt.org](https://scipopt.org):

- [online documentation](https://scipopt.org/doc) of the code with information how to get started,
- downloads of official releases and binaries for various platforms,
- release reports and scientific articles describing the algorithms behind SCIP,
- information how to cite SCIP when you use it in scientific publications,
- ...

For installation instructions have a look [here](INSTALL.md) or in the [online
documentation](https://scipopt.org/doc/html/INSTALL.php).

# Exact SCIP
This branch is a development version of the *numerically exact* variant of MIP solver SCIP.

The exact solving mode is based on the framework described in
> William J. Cook, Thorsten Koch, Daniel E. Steffy, Kati Wolter: [A hybrid branch-and-bound approach for exact rational mixed-integer programming.](https://doi.org/10.1007/s12532-013-0055-6) Math. Program. Comput. 5(3): 305-344 (2013)

That framework has been revised and extended by symbolic presolving, using the parallel presolving library [PaPILO](https://github.com/lgottwald/PaPILO), as well as an exact repair heuristic.
The results of this revision are described in
> Leon Eifler, Ambros Gleixner: [A computational status update for exact rational mixed integer programming](https://doi.org/10.1007/s10107-021-01749-5)  Math. Program. (2022).

Furthermore, the most recent version features separation and verification of safe Gomory mixed integer cuts.

Certificate printing is supported (although presolving cannot yet be verified), and can be checked using [VIPR](https://github.com/ambros-gleixner/VIPR).

The input formats that are currently supported for the exact solving mode are [ZIMPL](https://zimpl.zib.de/), MPS, and LP.

# Installation instructions for exact rational SCIP

This development version of exact SCIP can only be built using CMake, the Makefile system is not yet functional.

Dependencies which are optional for the floating point version but mandatory for exact SCIP:
* [GMP](https://gmplib.org/) for rational arithmetic in ZIMPL, SoPlex, SCIP, and PaPILO
* [Boost multiprecision library](https://www.boost.org/) (Version should be 1.70 or newer) for rationals in SCIP and PaPILO
* [MPFR](https://www.mpfr.org/) for approximating rationals with floating-point numbers in SCIP.

## Step by step guide:

### 1. Build SoPlex, PaPILO (optional), Zimpl (optional) :

```
# go into the directory
cd <one of zimpl/soplex/papilo>
mkdir build
cmake ..
make
cd ../..
```

### 2. Build SCIP:

```
cd scip
mkdir build
cmake .. -DSOPLEX_DIR=<PATH_TO_SOPLEX_BUILD_DIRECTORY> -DPAPILO_DIR=<PATH_TO_PAPILO_BUILD_DIRECTORY> -DZIMPL_DIR=<PATH_TO_ZIMPL_BUILD_DIRECTORY>
make
```

Test if everything works correctly with the MIPEX ctest (viprchk needs to be installed for the certificate checking to work)

```
ctest -R MIPEX
```

# USAGE

Enabling/disabling the exact solving mode is done with parameter `exact/enabled`, note that this has to be done before reading a problem instance. Further advanced parameters for exact solving can be set in the `exact` submenu.

Certificate printing can be enabled by setting `certificate/filename` to a non-default value. If cutting planes are enabled, then the certificate file needs to be completed using the `viprcomp` script prior to verifiaction.
