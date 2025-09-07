### Description

%{all_commits}

### Code review

* [ ] The code change is correct.
* [ ] The naming and place of new methods is clear and consistent **or** no new methods have been added.
* [ ] Emphasis settings are up to date.
* [ ] Ensure that new code and features use randomization where necessary to account for performance variability.
* [ ] The code is sufficiently documented.
* [ ] The coding style is OK, see https://scipopt.org/doc/html/CODE.php.

### Documentation and building

* [ ] The CHANGELOG is up to date (including API changes if present in this MR).
* [ ] The user documentation is up to date (doc/xternal.c, doc/inc/faq/, installation instructions, ...).
* [ ] Both build systems and makedist.sh are up to date. Especially, newly added, renamed or removed source files have been added to, renamed in or removed from src/CMakeLists.txt.
* [ ] If new CMake options have been added, non-default values are written to the log (see `# write log file of non-default parameter settings` in `CMakeLists.txt`).

### Testing

* [ ] ctest passes without errors (type some of `jenkins ctest {soplex master,soplex bugfix,cplex,gurobi,mosek,xpress,highs}`).
* [ ] The performance impact has been checked (type some of `jenkins performance {mip,minlp,pb} (quick|continue|)`), **or** the changed code will not be executed by default.
* [ ] The new code is sufficiently covered by tests (perhaps, new coverage settings or new unit tests have been added).
* Consider a debug run (type some of `jenkins debug {short,minlp,mip,pb}`).

### Does this merge request introduce an API change? :warning:

* [ ] No, **or** as far as possible, the code ensures backwards compatibility.
* [ ] No, **or** the `SCIP_APIVERSION` will be updated (use label 'public SCIP API'; in particular if new `SCIP_EXPORT` methods have been added).
  To update `SCIP_APIVERSION`, run `scripts/updateversion.py -a` on branch **master** after the changes of this MR have arrived in the master branch (if this is a merge into vX-minor, then first merge vX-minor into master or wait for the Jenkins job that does that).
* [ ] No parameter was added/deleted/changed **or** the MR is tagged with the label 'default parameter'. (This includes the addition of parameters by adding/deleting a plugin.)
