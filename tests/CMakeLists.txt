#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#@file    tests/CMakeLists.txt
#@brief   CMake lists file for including unit tests
#@author  Gregor Hendel

include(CTest)

#
# this directory contains unit tests for certain units of SCIP
#
# the build file generation recursively scans all src-directories for
# C-files and creates an executable that respects the given directory
# structure (relative to ./src/) in the build directory under tests/.
# The executable is dynamically linked to (and automatically depends on)
# the generated SCIP library.
#
# Finally, a test is added for every executable.
#

#
# define the C99 standard for older compilers (gcc < 5.0)
#
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED on)

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
    include_directories(${CRITERION_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../src)
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
        add_executable(${testName} src/${testSrc})

        #
        # link the executable against both criterion and the SCIP library target
        #
        target_link_libraries(${testName} ${CRITERION_LIBRARIES} libscip m)

        #
        # use the directory name (relative to src/) to build the executable
        #
        set_target_properties(${testName} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${dir})
        #
        # add the newly created executable to the dependencies of the unittest target
        #
        add_dependencies(unittests ${testName})

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
