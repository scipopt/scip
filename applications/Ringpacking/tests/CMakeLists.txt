#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      *
#*                                                                           *
#*  Licensed under the Apache License, Version 2.0 (the "License");          *
#*  you may not use this file except in compliance with the License.         *
#*  You may obtain a copy of the License at                                  *
#*                                                                           *
#*      http://www.apache.org/licenses/LICENSE-2.0                           *
#*                                                                           *
#*  Unless required by applicable law or agreed to in writing, software      *
#*  distributed under the License is distributed on an "AS IS" BASIS,        *
#*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
#*  See the License for the specific language governing permissions and      *
#*  limitations under the License.                                           *
#*                                                                           *
#*  You should have received a copy of the Apache-2.0 license                *
#*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#@file    tests/CMakeLists.txt
#@brief   CMake lists file for including ringpacking unit tests
#@author  Benjamin Mueller

include(CTest)

#
# define the C99 standard for older compilers (gcc < 5.0)
#
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED on)

# path to CMake modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../../cmake/Modules)

# find Criterion
find_package(CRITERION)

#
# unit tests depend on the presence of the Criterion framework
#
if( CRITERION_FOUND )

  #
  # glob recurse into the src/ subdirectory and collect all C-files
  #
  file(GLOB_RECURSE TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)

  #
  # besides Criterion, we also include the current source directory for test specific implementations,
  # and the source directory of the parent SCIP directory to be able to include C-sources of SCIP
  # into our unittests
  #
  include_directories(${CRITERION_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../../../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../../../tests)

  # link to math library if it is available
  find_library(LIBM m)
  if(NOT LIBM)
    set(LIBM "")
  endif()

  #
  # Prepare unit test executable for every source file
  #
  foreach(testSrc ${TEST_SRCS})
    #
    # extract the filename without an extension and directory (NAME_WE)
    # and extract the directory name (DIRECTORY) in dedicated variables.
    #
    get_filename_component(testName ${testSrc} NAME_WE)
    get_filename_component(dir ${testSrc} DIRECTORY)

    #
    # replace forward slashes by hyphens
    #
    string(REPLACE "/" "-" dirhyphen ${dir})
    set(testName "unittest-${dirhyphen}-${testName}")

    #
    # add compile target
    #
    add_executable(${testName} src/${testSrc}
      ../src/reader_rpa.c
      ../src/probdata_rpa.c
      ../src/pricer_rpa.c
      ../src/cons_rpa.c
      ../src/pattern.c
    )

    #
    # link the executable against both criterion and the SCIP library target
    #
    target_link_libraries(${testName} ${CRITERION_LIBRARIES} ${SCIP_LIBRARIES} ${LIBM})

    #
    # use the directory name (relative to src/) to build the executable
    #
    set_target_properties(${testName} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${dir})

    #
    # add a test to build this unit test
    #
    add_test(NAME ${testName}-build
      COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${testName}
      )

    #
    # avoid that several build jobs try to concurrently build the SCIP library
    # note that this ressource lock name is not the actual libscip target
    #
    set_tests_properties(${testName}-build
      PROPERTIES
      RESOURCE_LOCK libscip
      )

    #
    # Finally, add the executable as a test to the test framework
    #
    if( COVERAGE )
      set( ADDITIONAL_FLAGS "--no-early-exit" )
    else()
      set( ADDITIONAL_FLAGS "" )
    endif()

    add_test(NAME ${testName}
      COMMAND $<TARGET_FILE:${testName}> ${ADDITIONAL_FLAGS}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      )
    set_tests_properties(${testName}
      PROPERTIES
      DEPENDS ${testName}-build
      )
  endforeach(testSrc)
endif()
