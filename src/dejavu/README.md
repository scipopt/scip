# Compilation, Library, Tests
dejavu is a solver and C++ library for the fast detection and manipulation of [combinatorial symmetry](https://automorphisms.org/quick_start/symmetry/). 
Below, you can find some information on how to build the solver and include the library.
More detailed information can be found in our [get started guide](https://automorphisms.org/) or in the [full documentation](https://automorphisms.org/documentation/), which can also be built from the source code using [doxygen](https://www.doxygen.nl/).


## Compilation
Using *cmake*, the project should compile without any further dependencies:
```text
cmake . -DTEST=0
make
```
Compilation produces a binary *dejavu*. It accepts a DIMACS graph as input, and computes the automorphism group of the graph. For available options and more descriptions, please refer to our [guide](https://automorphisms.org/quick_start/standalone/).

## Use dejavu as a library
dejavu is a header-only library. You can simply add dejavu to your C++ project by including the respective header file: 
```cpp
#include "dejavu.h"
```

Note that currently, dejavu requires to be *compiled with C++ version 14*. For a more thorough description, please refer to our [guide](https://automorphisms.org/quick_start/cpp_api/).

By default, dejavu is compiled without assertions. We recommend activating assertions for debugging purposes (by adding the definition `DEJDEBUG`). Assertions do however slow the code considerably.

## Running the tests
Using *cmake*, a test target `dejavu_test` can be produced. The test target depends on [googletest](https://github.com/google/googletest) and a set of [test graphs](https://automorphisms.org/graphs/graphs.zip) . The following set of commands should download all of the dependencies and make the test target:
```text
cmake . -DTEST=1
make dejavu_test
```

Then, you may run the tests:
```text
./dejavu_test
```

Note that when running cmake with `-DTEST=1`, the `dejavu` target will also be compiled with assertions on (and hence, may run slower).
