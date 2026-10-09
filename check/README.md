# Testsuite

Please make sure you have read the section about testruns in the scip documentation https://scipopt.org/doc/html/TEST.php before you continue reading.

In general the execution of the testruns is the following:

- To start a testrun call `make test` for a local one, `make testcluster` for one on the cluster.
  (These are for running scip tests but variants exist for other solvers, see below)
- The `Makefile` will call the corresponding `check_<solver>.sh` or `check_cluster_<solver>.sh` script
  with the correct variables. These scripts will
  + configure the environment variables for local and cluster runs `configuration_set.sh`, `configuration_cluster.sh`
  + configure the test output files such as the .eval, the .tmp and the .set files `configuration_logfiles.sh`
  + cluster run scripts have a `waitcluster.sh` script in the loop that makes the scripts wait
    instead of overloading the cluster queue with jobs.
  + run `evalcheck_cluster.sh` to evaluate and save any old testruns that lay around
  + `configuration_tmpfile_setup_{cbc,cplex,gurobi,scip,xpress}.sh`
    reset and fill a tmpfile to run the solver with; tmpfile will set correct limits, read in settings, and control
    display of the solving process
  + with `run.sh` (or `run.*.sh`) run or submit the jobs for the current testrun
  + some local `check_*.sh` scripts run `evalcheck.sh` or `evalcheck_cluster.sh` automatically
    after the testrun, for the others and for clusterruns the user has to do it themselves;
    the `evalcheck*.sh` scripts evaluate the testrun and concatenate the individual
    logfiles; the naming here are a bit misleading, the main difference between the two files is that
    the `_cluster` script is more involved and cares about cleaning up after itself; both of them can
    be used to evaluate a `check.*.eval` file; Recommended is to use the command
    `evalcheck_cluster.sh results/check.*.eval` from the `check` folder

## Local make targets

### Run a local testrun with SCIP

make test
  - `check.sh`
    + `configuration_set.sh`
    + `configuration_logfiles.sh`
    + `evalcheck_cluster.sh`
    + `configuration_tmpfile_setup_{cbc,cplex,gurobi,scip,xpress}.sh`
    + `run.sh`

### Run a local testrun with FiberSCIP

make testfscip
  - `check_fscip.sh`
    + `configuration_set.sh`
    + `configuration_logfiles.sh`
    + `evalcheck_cluster.sh`
    + `run_fscip.sh`

### Run a local testrun with GLPK, Mosek, Symphony

make test{glpk,mosek,symphony}
  - `check_{glpk,mosek,symphony}.sh`
    +`getlastprob.awk`
    +`evalcheck.sh`

### Generate a coverage report

make coverage
  - `check_coverage.sh`
    + `evalcheck.sh`

### Count feasible solutions

make testcount
  - `check_count.sh`
    + `getlastprob.awk`
    + `evalcheck_count.sh`
       - `configuration_solufile.sh`
       - `check_count.awk`

## Cluster make targets

### Run a clusterrun with SCIP, CBC, CPLEX, Gurobi, Xpress

make testcluster{,cbc,cpx,gurobi,xpress}
  - `check_cluster.sh`
    + `configuration_cluster.sh`
    + `configuration_set.sh`
    + `waitcluster.sh`
    + `configuration_logfiles.sh`
    + `configuration_tmpfile_setup_{cbc,cplex,gurobi,scip,xpress}.sh`
    + `run.sh`

### Run a clusterrun with FiberSCIP

make testclusterfscip
  - `check_cluster_fscip.sh`
    + `configuration_cluster.sh`
    + `configuration_set.sh`
    + `waitcluster.sh`
    + `configuration_logfiles.sh`
    + `run_fscip.sh`

### Run a clusterrun with Mosek

make testclustermosek
  - `check_cluster_mosek.sh`
    + `configuration_cluster.sh`
    + `waitcluster.sh`
    + `run.sh`

### Start a GAMS testrun on the cluster

make testgamscluster
  - `check_gamscluster.sh`
    + `schulz.sh`
    + `configuration_solufile.sh`
    + `waitcluster.sh`
    + `rungamscluster.sh`
    + `finishgamscluster.sh`
      - evalcheck_gamscluster.sh
        + `configuration_solufile.sh`
        + `check_count.awk`

## Other scripts

### Generate comparison of two testruns

- `allcmpres.sh`
  + `cmpres.awk`

### Compares averages of several SCIP result files

- `average.sh`
  + `average.awk`

### Compute averages over instances for different permutations

- `permaverage.sh`
  + `permaverage.awk`

### Compare different versions of runs with permuations

- `permcmpresall.sh`
  + `permcmpresall.awk`

# Files

## AWK files

- `check_*.awk` parses and checks the check files of a testrun and output a .res file table

- `check_count.awk` counts feasible solutions

- `getlastprob.awk` gets the last problem of a check.*-outfile

- `cmpres.awk` check comparison report generator - compare two res files

- `average.awk` computes averages of several SCIP result files

- `permaverage.awk` compute averages over instances for different permutations

- `permcmpresall.awk` compares different versions of runs with permutations

## Bash scripts

### Configuration

- `configuration_tmpfile_setup_{cbc,cplex,gurobi,scip,xpress}.sh`
  resets and fills a batch file TMPFILE to run the solver with;
  tmpfile will set correct limits, read in settings, and control display of the solving process

- `configuration_logfiles.sh` configures the right test output files such as the .eval, the .tmp and the .set files to run a test on

- `configuration_solufile.sh` configures SOLUFILE env variable from name of testset

- `configuration_cluster.sh` configures environment variables for cluster runs;
  it is to be invoked inside a `check_cluster*.sh` script

- `configuration_set.sh` configures environment variables that are needed for test runs both on the cluster and locally;
  it is to be invoked inside a `check(_cluster)*.sh` script
  + `configuration_solufile.sh`

### Running

- `run.sh` executes EXECNAME on one instance and produces the logfiles;
  can be executed either locally or on a cluster node;
  is to be invoked inside a `check(_cluster)*.sh` script

- `run_fscip.sh` executes FiberSCIP on one instance and produces the logfiles;
  can be executed either locally or on a cluster node;
  is to be invoked inside a `check(_cluster)*.sh` script

- `rungamscluster.sh` executes GAMS on one instance and producing the logfiles;
  can be executed either locally or on a cluster node;
  is to be invoked inside a `check(_cluster)*.sh` script

### Evaluation

- `evalcheck.sh` evaluates one or more testrun by each checking the check*-outfile;
  is to be invoked inside a `check*.sh` script
  + `evaluate.sh`
    + `configuration_solufile.sh`
    + `check_*.awk`

- `evalcheck_cluster.sh` evaluates a testrun and concatenates the individual logfiles, possibly uploads to Rubberband
  + `evaluate.sh`
    + `configuration_solufile.sh`
    + `check_*.awk`

- `evaluate.sh` calls check_*.awk, depending on the solver used, on the testrun files and writes the output in a .res file
  + `configuration_solufile.sh`
  + `check_*.awk`

- `finishgamscluster.sh` cleans up after GAMS testrun
  + `evalcheck_gamscluster.sh`
    - `configuration_solufile.sh`
    - `check_count.awk`

### Helper

- `waitcluster.sh`:
  in order to not overload the cluster, no jobs are submitted if the queue is too full;
  instead, this script waits until the queue load falls under a threshold and returns
  for the calling script to continue submitting jobs

- `schulz.sh`: Supervises processes and sends given signals via kill when elapsed time exceeds given thresholds

## Directories and other files

- `CMakeLists.txt` is part of the cmake system, this file's main purpose is to add tests

- the directory `mipstarts/` contains some files used by `CMakeLists.txt` to generate tests

- in the directory `coverage/` are all necessary files located that are used in producing the coverage report for SCIP

- the directory `interactiveshell/` contains files for `CMakeLists.txt` that get configured to batchfiles to be executed by SCIP

- in the `testset/` directory reside the `*.test` and `*.solu` files to be specified via `TEST=short`;
  a `.test` file lists problem files, one file per line, absolute and relative paths to its location;
  a `.solu` file with the same basename as the `.test` file contains information about feasibility and best known objective value;
  it is optional for a testrun

- in the `instances` directory are some example instances that can be solved with SCIP
