# Testsuite

Please make sure you have read the section about testruns in the scip documentation https://scipopt.org/doc/html/TEST.php before you continue reading.

In general the execution of the testruns is the following:

- To start a testrun call `make test` for a local one, `make testcluster` for one on the cluster.
  (These are for running scip tests but variants exist for other solvers, see below)
- The `Makefile` will call the corresponding `check_<solver>.sh` or `check_cluster_<solver>.sh` script
  with the correct variables. These scripts will
  + configure the environment variables for local and cluster runs `. configuration_set.sh`, `. configuration_cluster.sh`
  + configure the test output files such as the .eval, the .tmp and the .set files `. configuration_logfiles.sh`
  + cluster run scripts have a `waitcluster.sh` script in the loop that makes the scripts wait
    instead of overloading the cluster queue with jobs.
  + run `evalcheck_cluster.sh` to evaluate and save any old testruns that lay around.
  + `. configuration_tmpfile_setup_{cbc,cplex,gurobi,scip,xpress}.sh`
    reset and fill a tmpfile to run the solver with. Tmpfile will set correct limits, read in settings, and control
    display of the solving process.
  + with `run.sh` (or `run.*.sh`) run or submit the jobs for the current testrun
  + some local `check_*.sh` scripts run `evalcheck.sh` or `evalcheck_cluster.sh` automatically
    after the testrun, for the others and for clusterruns the user has to do it themselves.
    The `evalcheck*.sh` scripts evaluate the testrun and concatenate the individual
    logfiles. The naming here are a bit misleading, the main difference between the two files is that
    the `_cluster` script is more involved and cares about cleaning up after itself. Both of them can
    be used to evaluate a `check.*.eval` file. Recommended is to use the command
    `evalcheck_cluster.sh results/check.*.eval` from the `check` folder.

## local make targets

### run a local testrun with scip

make test
  - `check.sh`
    + `. configuration_set.sh`
    + `. configuration_logfiles.sh` END
    + `evalcheck_cluster.sh` END
    + `. configuration_tmpfile_setup_{cbc,cplex,gurobi,scip,xpress}.sh` END
    + `run.sh` END

### run a local testrun with fscip

make testfscip
  - `check_fscip.sh`
    + `. configuration_set.sh`
    + `. configuration_logfiles.sh` END
    + `evalcheck_cluster.sh` END
    + `run_fscip.sh` END

### run a local testrun with glpk, mosek, symphony

make test{glpk,mosek,symphony}
  - `check_{glpk,mosek,symphony}.sh`
    +`getlastprob.awk` END
    +`evalcheck.sh`

### generate a coverage report

make coverage
  - `check_coverage.sh`
    + `evalcheck.sh`

### count feasible solutions

make testcount
  - `check_count.sh`
    + `getlastprob.awk` END
    + `evalcheck_count.sh`
       - `. configuration_solufile.sh` END
       - `check_count.awk` END

## cluster make targets

### run a clusterrun with scip, cbc, cpx, gurobi, xpress

make testcluster{,cbc,cpx,gurobi,xpress}
  - `check_cluster.sh`
    + `. configuration_cluster.sh`
    + `. configuration_set.sh`
    + `waitcluster.sh` END
    + `. configuration_logfiles.sh` END
    + `. configuration_tmpfile_setup_{cbc,cplex,gurobi,scip,xpress}.sh` END
    + `run.sh` END

### run a clusterrun with fscip

make testclusterfscip
  - `check_cluster_fscip.sh`
    + `. configuration_cluster.sh`
    + `. configuration_set.sh`
    + `waitcluster.sh` END
    + `. configuration_logfiles.sh` END
    + `run_fscip.sh` END

### run a clusterrun with mosek

make testclustermosek
  - `check_cluster_mosek.sh`
    + `. configuration_cluster.sh`
    + `waitcluster.sh` END
    + `run.sh` END

### start a gams testrun on the cluster

make testgamscluster
  - `check_gamscluster.sh`
    + `schulz.sh` END
    + `. configuration_solufile.sh` END
    + `waitcluster.sh` END
    + `rungamscluster.sh` END
    + `finishgamscluster.sh`
      - evalcheck_gamscluster.sh
        + `. configuration_solufile.sh` END
        + `check_count.awk` END

## other Scripts

### Generate comparison of two testruns.

- `allcmpres.sh`
  + `cmpres.awk` END

### compares averages of several SCIP result files

- `average.sh`
  + `average.awk` END

### compute averages over instances for different permuations

- `permaverage.sh`
  + `permaverage.awk` END

### compare different versions of runs with permuations

- `permcmpresall.sh`
  + `permcmpresall.awk` END

# Files

## AWK files

- `check_*.awk` END
  Parse and check the check files of a testrun and output a .res file table.

- `check_count.awk`
  Count feasible solutions

- `getlastprob.awk` END
  Get the last problem of a check.*-outfile

- `cmpres.awk` END
  Check Comparison Report Generator - Compare two res files.

- `average.awk` END
  compute averages of several SCIP result files

- `permaverage.awk` END
  compute averages over instances for different permuations

- `permcmpresall.awk` END
  compare different versions of runs with permuations

## Makefiles

- `wakeup-slurm` END
  A Makefile that can wake up sleeping slurm nodes and queues before submitting tests to avoid errors.

## Bash Scripts

### Configuration

- `. configuration_tmpfile_setup_{cbc,cplex,gurobi,scip,xpress}.sh` END
  Resets and fills a batch file TMPFILE to run the solver with.
  Tmpfile will set correct limits, read in settings, and control display of the solving process.

- `. configuration_logfiles.sh` END
  Configures the right test output files such as the .eval, the .tmp and the .set files to run a test on.

- `. configuration_solufile.sh` END
  Configures SOLUFILE env variable from name of testset.

- `. configuration_cluster.sh`
  Configures environment variables for cluster runs.
  It is to be invoked inside a `check_cluster*.sh` script.
  It calls `wakeup-slurm` to the respective queue.
  + `wakeup-slurm` END

- `. configuration_set.sh`
  Configures environment variables that are needed for test runs both on the cluster and locally.
  It is to be invoked inside a `check(_cluster)*.sh` script.
  + `. configuration_solufile.sh` END

### Running

- `run.sh` END
  The script executing EXECNAME on one instance and producing the logfiles.
  Can be executed either locally or on a cluster node.
  Is to be invoked inside a `check(_cluster)*.sh` script.

- `run_fscip.sh` END
  The script executing fscip on one instance and producing the logfiles.
  Can be executed either locally or on a cluster node.
  Is to be invoked inside a `check(_cluster)*.sh` script.

- `rungamscluster.sh` END
  The script executing gams on one instance and producing the logfiles.
  Can be executed either locally or on a cluster node.
  Is to be invoked inside a `check(_cluster)*.sh` script.

### Evaluation

- `evalcheck.sh`
  Evaluates one or more testrun by each checking the check*-outfile.
  Is to be invoked inside a `check*.sh` script.
  + `. evaluate.sh`
    + `. configuration_solufile.sh` END
    + `check_*.awk` END

- `evalcheck_cluster.sh` END
  Evaluates a testrun and concatenates the individual logfiles, possibly uploads to rubberband.
  + `. evaluate.sh`
    + `. configuration_solufile.sh` END
    + `check_*.awk` END

- `. evaluate.sh`
  Depending on the solver used for the testrun calls corresponding check_*.awk on the testrun files and writes the output in a .res file.
  + `. configuration_solufile.sh` END
  + `check_*.awk` END

- `finishgamscluster.sh`
  Cleans up after gams testrun.
  + evalcheck_gamscluster.sh
    - `. configuration_solufile.sh` END
    - `check_count.awk` END

### helper

- `waitcluster.sh` END
  In order to not overload the cluster, no jobs are submitted if the queue is too full
  instead, this script waits until the queue load falls under a threshold and returns
  for the calling script to continue submitting jobs.

### uncategorized

- `schulz.sh` END
  Supervises processes and sends given signals via kill when elapsed time exceeds given thresholds.

## folders and other files

- `CMakeLists.txt` is part of the cmake system, this file's main purpose is to add tests.

- The directory `mipstarts/` contains some files used by `CMakeLists.txt` to generate tests.

- In the directory `coverage/` are all necessary files located that are used in producing the coverate report for SCIP.

- The directory `interactiveshell/` contains files for `CMakeLists.txt` that get configured to batchfiles to be executed by SCIP.

- In the `testset/` directory reside the *.test and *.solu files to be specified via `TEST=short`.
  A `.test` file lists problem files, one file per line, absolute and relative paths to its location.
  A `.solu` file with the same basename as the `.test` file contains information about feasibility and best known objective value.
  It is optional for a testrun.

- In the `instances` direcotry are some example instances that can be solved with SCIP.
