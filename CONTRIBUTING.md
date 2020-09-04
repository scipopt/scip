# Contributing

Before making a contribution to SCIP, you'll need to agree to the terms outlined
in the [Contributor License Agreement](CLA).

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
| LP Solver | SoPlex [c155d8895e3](https://git.zib.de/integer/soplex/commit/c155d8895e3)|
| OS | Linux |
| Architecture | x86 64-bit |
| Compile command | `make OPT=dbg LPS=spx` or `cmake -DCMAKE_BUILD_TYPE=Debug`|
| Compiler | gcc-4.8.4 |
| Attachments | |

## Steps to Reproduce

## Extended Description
```

When creating SCIP issues in the SCIP repository, all or part of the SCIP commit hash
will auto-link to the commit. Commit hashes from other repositories, however,
won't autolink, so be sure to include an explicit link. Also, please apply the
~bug label when using the above template.

# Code Analysis

Currently, we use (at least) four tools to perform a static analysis of the code.
These tools have helped to detect many errors that would not have been discovered otherwise.
The code analysis tools `pclint` (older versions: `flexelint`), `cppcheck`, and the Clang Static Analyzer (`scan-build`)
investigate the source code and possibly output warnings.
Moreover, from time to time, `coverity` is used to analyze different pathways of the code.
In addition, one can use `valgrind` to detect memory leaks and access problems.
These tools (except `coverity`) are also automatically run in the gitlab CI for merge requests.

## Linting

In the following, we will give some hints on how to deal with warnings from `pclint`/`flexelint` (`lint` in the following).
First, it should be noted that most warnings have a good reason, but the static analysis of code has its limits
and `lint` does not always understand the code we write.
Many of the `lint` warnings are globally disabled in the file `lint/lint.lnt` and `pclint/pclint.lnt`.
One should take care when adding suppressions to these files, since the warnings might be of value.

One example is warning ``715 - named parameter symbol of function symbol not subsequently referenced``.
In a framework like SCIP, in which callback functions play an important role, this warning is triggered quite often,
since some function parameters might not be needed. However, we do not want to disable this warning globally, since it
often helps to identify function parameters that are useless, e.g., by a copy-paste error.

`Lint` warnings can be disabled by adding a comment in the particular line indicated by the warning, like ``/*lint !e527*/``.
Note that `lint` parses asserts, so adding ``assert( v != NULL );``, for example,  can convince `lint`` that ``v`` is not a NULL-pointer.
Currently, `lint` does not always understand more complicated asserts, so it might be useful to add redundant asserts
at particular places.

Some other warnings that often occur are the following:

- Warnings connected to type conversions: for example, ``unsigned int`` is converted to ``int``
  or ``size_t`` to ``int``. In principle, this comes with a reduction of precision, but this might be o.k. in many
  contexts, especially, since ``int`` is the default type for counting variables. In general, adding an explicit type
  cast will silence `lint`, but should be used carefully.

- One further warning that often occurs is ``777 - testing floating point values for equality``. This often
  indicates that one should have used something like ``SCIPisEq()`` instead.

- The warning ``534 - ignoring return value of function symbol`` often indicates that one has forgotten ``SCIP_CALL()`` around
  a function. If this is not the case often something like ``(void) SCIPsnprintf(...)`` helps.


## cppcheck

`cppcheck` suppressions can be performed inline, but sometimes it is preferable to disable a category of checks in `suppressions.txt`.
