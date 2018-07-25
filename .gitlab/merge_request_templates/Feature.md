**Description:**
[Add a description of the introduced changes here.]

**Code review:**
* [ ] Is the code sufficiently documented? Are CHANGELOG entries added? Is the coding style OK?
* [ ] Is the naming and place of new methods clear and consistent?
* [ ] Is additional user documentation in doc/xternal.c, doc/inc/faq/, INSTALL file, etc. added (if necessary)?
* [ ] Are coverage settings added to test new code (if necessary)?
* [ ] Have unit tests been added (if necessary)?
* [ ] Are new files added to makedist.sh and both build systems?
* [ ] Does ctest pass (jenkins ctest ...)?
* [ ] Has performance impact been checked on mi(nl)pdev-solvable (if activated by default)?

**:warning: If this merge introduces an API change:**
* [ ] Look for a satisfactory solution that ensures backwards compatibility.
* [ ] Increase SCIP_APIVERSION after the merge.
* [ ] Interface changes must be documented in doc/xternal.c and the CHANGELOG.
* [ ] Tag this MR with the label 'default parameter' and inform one of the developers responsible for SAP (default: Jakob) if a parameter was added/deleted/changed.