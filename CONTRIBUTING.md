# Contributing
Before making a contribution to SCIP, you'll need to agree to the terms outlined
in the [Contributor License Agreement](CLA)

## Workflow
Please see the [wiki](https://wiki.zib.de/confluence/display/SCIPOPTSUITE/Development+Workflow)
for information about our contribution workflow with git. If you don't have
access to this page, request it from @bzfgleix.

Additionally, you should install git hooks that check your commit messages, e.g. for correct formatting:
[localhooks](https://git.zib.de/optimization/localhooks)

## Bug Reports
Use the following sample template in the Issue description to ensure consistency
across bug reports.

```
|||
| --- | --- |
| Affected Versions | 3.2.2, 3.2.1 |
| SCIP Commit | 8eec1de1206 |
| LP Solver | [c155d8895e3](https://git.zib.de/integer/soplex/commit/c155d8895e3)|
| OS | Linux |
| Architecture | x86 64-bit |
| Compile command | `make OPT=dbg LPS=spx` |
| Compiler | gcc-4.8.4 |
| Attachments | |

## Steps to Reproduce

## Extended Description
```

When creating SCIP issues in the SCIP repository, all or part of the commit hash
will auto-link to the commit. Commit hashes from other repositories, however,
won't autolink, so be sure to include an explicit link. Also, please apply the
~bug label when using the above template.

## Linting

Currently, `flexelint` and `cppcheck` are the preferred light-weight static analysis tooling that we use. `flexelint` supressions are done inline. `cppcheck` supresssions can also be performed inline, but it's typically preferably to disable a category of checks in `supressions.txt`
