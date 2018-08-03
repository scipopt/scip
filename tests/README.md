# SCIP Unit Tests

Write and run unit tests for SCIP.

- [Overview](#overview)
- [Write](#write)
  - [Examples](#examples)
- [Compile](#compile)
- [Run](#run)
- [Debug](#debug)

## Overview

A unit test is an automated piece of code that invokes a unit of work in the system and then checks a single assumption about the behavior of that unit of work. The SCIP Unit Test Suite leverages [Criterion](http://criterion.readthedocs.io/en/master/) as the testing framework and [ctest](https://cmake.org/cmake/help/v2.8.8/ctest.html) as the runner. The SCIP Unit Test Suite is very much in a state of development. Check out the [unit test suite milestone](https://git.zib.de/integer/scip/milestones/2) for more information.

## Write

Tests are organized into topic-specific directories in `src`. When writing new tests, find the directory that best suites your test, or create one if it doesn't already exist. For example, if a test is meant to illustrate a bug, place is in `src/bugs/`. Use `#include "include/scip_test.h"` to access Criterion and the `SCIP_CALL` macro. Ensure that this is the **last** included header.

**NOTE** If your test needs `SCIP` code (eg, you are implementing a constraint handler in your test, see `src/cons/cons.c`), place `#include "include/scip_test.h"` after the SCIP code.

Criterion comes with [fixtures](http://criterion.readthedocs.io/en/master/starter.html?highlight=fixture#fixtures) and [asserts](http://criterion.readthedocs.io/en/master/assert.html) built-in, and also supports [parameterized tests](http://criterion.readthedocs.io/en/master/parameterized.html).

### Examples

Here are some test examples that can help you get started writing unit tests.

| Example Type| Location |
| ------ | ------ |
| catch a signal | `src/bugs/depthlevel.c` |
| parameterized test | unittest_framework_tmp branch, `src/cons/expr/simplify.c` |
| check stdout | unittest_framework_tmp branch, `src/cons/expr/walk.c` |

## Compile

Smart test discovery is already built into the `Makefile`, so anything in `src` (at any level of nesting) will be detected and compiled into the equivalent path in the `bin` directory. Also, the `Makefile` generates the test "makefile" for `ctest`. There should never be a reason to directly modify any Makefile unless you are hacking on the SCIP Unit Test Suite.

The easiest way to compile and run the tests is:

```
make
```

**NOTE** `make` will read the options used for building `SCIP` from the binary. It uses `scip --version` to find out the options.
If `SHARED=true`, or `OPT=dbg`  were not used when compiling `SCIP`, `make` will end with a proper error. If the binary is not found, it will
assume `SHARED=true` and `OPT=dbg`.

This command will check for [Criterion](http://criterion.readthedocs.io/en/master/) in ./Criterion, download and install it if not found, and compile and run all tests in `src/`.
If you already have installed Criterion on you system, execute `touch Criterion` or `mkdir Criterion` before calling make.

**NOTE** Some tests might need to include c files from SCIP. For tests to be recompilied the included c file gets recompiled, run `make depend`.

## Run

See above for the easiest way to compile and run tests. For simply running tests:

```
make
```

This creates `CTestTestfile.cmake` with a list of the test to run and then calls `ctest --output-on-failure`. By default, tests in `src/bugs/` are not compiled or run since they take a long time. To compile and run them:

```
make BUGS=true
```

You can also run a single test, e.g. `
```
 >> ./bin/cons/quadratic/gauge.linux.x86_64.gnu.dbg.spx2
```

Note, that parameterized tests will not work on systems that have address
space layout randomization (ASLR) enabled. One can disable ASLR for a
specific process (and its children) by calling it in a modified environment, e.g.,
```
 >> setarch `uname -m` -R ./bin/cons/quadratic/gauge.linux.x86_64.gnu.dbg.spx2
```

This is the approach that is also followed by the Makefile when running
the whole test suite.

Alternatively, one can disable ASLR system-wide (requires root access):
```
 >> sudo echo 0 > /proc/sys/kernel/randomize_va_space
```

TODO: Define a policy for moving/removing tests in `src/bugs` once the bugs are fixed.

## Debug (up to Criterion 2.2.2)

If a test fails, use `gdb` to debug. For example:

```
 >> ./bin/cons/quadratic/gauge.linux.x86_64.gnu.dbg.spx2
         [----] src/cons/quadratic/gauge.c:112: Assertion failed: gauge unavailable, pointless to continue
         [FAIL] separation::gauge: (0.00s)
         [====] Synthesis: Tested: 1 | Passing: 0 | Failing: 1 | Crashing: 0
```

The test suite is `separation` and the test name is `gauge`. To debug:

```
>> gdb --args bin/cons/quadratic/gauge.linux.x86_64.gnu.dbg.spx2 --single separation::gauge
(gdb) br src/cons/quadratic/gauge.c:112
```

Criterion by default prints all of the critical debugging information (test_suite::test_name, file and line number were to break). When a test crashes, there is no need to `break` in `gdb`.

## Debug (with Criterion 2.3 and later)

If a test fails, one can use `gdb` or `undodb-gdb` to debug. For example:

```
 >> ./bin/cons/quadratic/gauge.linux.x86_64.gnu.dbg.spx2
         [----] src/cons/quadratic/gauge.c:112: Assertion failed: gauge unavailable, pointless to continue
         [FAIL] separation::gauge: (0.00s)
         [====] Synthesis: Tested: 1 | Passing: 0 | Failing: 1 | Crashing: 0
```

The test suite is `separation` and the test name is `gauge`. To debug with `gdb` write in one terminal:

```
>> bin/cons/quadratic/gauge.linux.x86_64.gnu.dbg.spx2 --filter *gauge* --debug
```
This will start a `gdbserver`. To connect, in another terminal use
```
>> gdb bin/cons/quadratic/gauge.linux.x86_64.gnu.dbg.spx2 -ex "target remote localhost:1234"
```

If one doesn't want to use a `gdbserver` use:
```
>> bin/cons/quadratic/gauge.linux.x86_64.gnu.dbg.spx2 --filter *gauge* --debug=idle
```
This will give the PID of the process which can then be attached to a `undodb-gdb` or `gdb` session with
```
>> gdb --pid <pid-number>
```

After this, execute `continue` twice in gdb.
