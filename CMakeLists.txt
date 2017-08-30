cmake_minimum_required(VERSION 3.0)

set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}")

project(SCIP)

set(SCIP_VERSION_MAJOR 4)
set(SCIP_VERSION_MINOR 0)
set(SCIP_VERSION_PATCH 1)
set(SCIP_VERSION_SUB 0)
set(SCIP_VERSION_API 0)

set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "${SCIP_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${SCIP_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${SCIP_VERSION_PATCH}")
set(CPACK_PACKAGE_VENDOR "Zuse Institute Berlin")
include(CPack)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# path to e.g. findGMP module
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

option(SHARED "Build shared libraries" ON)
set(BUILD_SHARED_LIBS ${SHARED})
message(STATUS "Build shared libraries: " ${SHARED})

# make 'Release' the default build type
if(CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

option(ZLIB "should zlib be linked" ON)
option(READLINE "should readline be linked" ON)
option(GMP "should gmp be linked" ON)
option(ZIMPL "should zimpl be linked" ON)
option(IPOPT "should ipopt be linked" ON)
option(WORHP "should worhp be linked" OFF)
option(PARASCIP "should scip be compiled thread safe" OFF)
option(LPSCHECK "double check SoPlex results with CPLEX" OFF)
option(LEGACY "build SCIP/SoPlex in legacy mode" OFF)
option(NOBLKBUFMEM "should block and buffer memory be disabled" OFF)
option(NOBLKMEM "should block and buffer memory be disabled" OFF)
option(NOBUFMEM "should block and buffer memory be disabled" OFF)
option(DEBUGSOL "should the debug solution mechanism be enabled" OFF)
option(COVERAGE "enable coverage support" OFF)
SET(COVERAGE_CTEST_ARGS "" CACHE STRING "additional ctest arguments for coverage")
option(MT "use static runtime libraries for Visual Studio compiler" OFF)
option(CXXONLY "use a c++ compiler for all source files" OFF)

set(TPI none CACHE STRING "options for thread support library")  #create the variable
set_property(CACHE TPI PROPERTY STRINGS none tny omp)  #define list of values GUI will offer for the variable

set(EXPRINT cppad CACHE STRING "options for expression interpreter")  #create the variable
set_property(CACHE EXPRINT PROPERTY STRINGS none cppad )  #define list of values GUI will offer for the variable

set(LPS spx CACHE STRING "options for LP solver")  #create the variable
set_property(CACHE LPS PROPERTY STRINGS spx cpx grb glop msk )  #define list of values GUI will offer for the variable

#set options for memory management
if( NOBLKBUFMEM )
   set(NOBLKMEM ON FORCE)
   set(NOBUFMEM ON FORCE)
endif()

if(NOBLKMEM)
   add_definitions(-DBMS_NOBLOCKMEM)
endif()

if(NOBUFMEM)
   add_definitions(-DSCIP_NOBUFFERMEM)
endif()

if(DEBUGSOL)
   add_definitions(-DWITH_DEBUG_SOLUTION)
endif()

#set the correct rpath for OS X
set(CMAKE_MACOSX_RPATH ON)

#set defines for Windows
if(WIN32)
    add_definitions(-DNO_SIGACTION)
    add_definitions(-DNO_STRTOK_R)
endif()
if(MSVC)
    add_definitions(/W4)
    add_definitions(/wd4100)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

# Visual Studio compiler with static runtime libraries
if(MSVC AND MT)
    set(CMAKE_CXX_FLAGS_RELEASE "/MT")
    set(CMAKE_CXX_FLAGS_DEBUG "/MTd")
    set(CMAKE_C_FLAGS_RELEASE "/MT")
    set(CMAKE_C_FLAGS_DEBUG "/MTd")
endif()

#set expression interpreter file that should be used
if(EXPRINT STREQUAL "cppad")
    set(exprinterpret nlpi/exprinterpret_cppad.cpp)
elseif(EXPRINT STREQUAL "none")
    set(exprinterpret nlpi/exprinterpret_none.c)
else()
    message(FATAL_ERROR "EXPRINT option has wrong value")
endif()

#if tpi is not none force parascip to on
if(NOT (TPI STREQUAL "none"))
   set(PARASCIP ON FORCE)
endif()

#set sources files, libraries and defines for tpi
if(TPI STREQUAL "none")
    set(tpisources tpi/tpi_none.c)
    set(THREAD_LIBRARIES "")
    add_definitions(-DTPI_NONE)
elseif(TPI STREQUAL "tny")
    add_definitions(-DTPI_TNYC)
    set(tpisources tpi/tpi_tnycthrd.c tinycthread/tinycthread.c)
    find_package(Threads REQUIRED)
    set(THREAD_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
elseif(TPI STREQUAL "omp")
    add_definitions(-DTPI_OMP)
    set(tpisources tpi/tpi_openmp.c)
    find_package(OpenMP REQUIRED)
    set(THREAD_LIBRARIES "")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(FATAL_ERROR "TPI option has wrong value")
endif()

#set nparascip define if parascip is off
if(NOT PARASCIP)
    add_definitions(-DNPARASCIP)
endif()

#check if legacy mode must be used
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.7")
        set(LEGACY ON)
    endif()
endif()
message(STATUS "LEGACY mode for old compilers: ${LEGACY}")

# create a target for updating the current git hash
file(WRITE ${CMAKE_BINARY_DIR}/scip_update_githash.cmake "
find_program(GIT git)
if(GIT)
   execute_process(
      COMMAND \${GIT} describe --always --dirty
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE GITHASH OUTPUT_STRIP_TRAILING_WHITESPACE)
   if(GITHASH STREQUAL \"\")
      set(GITHASH \"------\")
   endif()
   string(REGEX REPLACE \"^.*-g\" \"\" GITHASH \${GITHASH})
else()
   set(GITHASH \"------\")
endif()
if(GITHASH STREQUAL \"------\")
   # read stored hash from file
   file(STRINGS \${DST} GITHASH)
   string(REGEX REPLACE \"#define SCIP_GITHASH \" \"\" GITHASH \${GITHASH})
else()
   file(WRITE \${DST} \"#define SCIP_GITHASH \\\"\${GITHASH}\\\"\")
endif()
message(STATUS \"Git hash: \" \${GITHASH})
")
add_custom_target(scip_update_githash
                  COMMAND ${CMAKE_COMMAND} -DDST=${PROJECT_SOURCE_DIR}/src/scip/githash.c
                                           -P ${CMAKE_BINARY_DIR}/scip_update_githash.cmake)

# create a target for updating the current build configuration
set(NEWLINE "\\\\n")
file(WRITE ${CMAKE_BINARY_DIR}/scip_update_buildflags.cmake "
file(WRITE \${DST} \"#define SCIP_BUILDFLAGS \\\" ARCH=\${CMAKE_SYSTEM}${NEWLINE} COMP=\${CMAKE_SYSTEM_PROCESSOR}${NEWLINE} BUILD=\${CMAKE_BUILD_TYPE}\\\"\")
")
add_custom_target(scip_update_buildflags
                  COMMAND ${CMAKE_COMMAND} -DDST=${PROJECT_SOURCE_DIR}/src/scip/buildflags.c
                                           -P ${CMAKE_BINARY_DIR}/scip_update_buildflags.cmake)

add_definitions(-DWITH_SCIPDEF)



# use C++11 standard
set(CMAKE_CXX_STANDARD 11)

# ZLIB
if(ZLIB)
    find_package(ZLIB)
endif()
if(ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    add_definitions(-DWITH_ZLIB)
else()
    set(ZLIB_LIBRARIES "")
endif()

# Readline
if(READLINE)
    find_package(Readline)
endif()
if(READLINE_FOUND)
    include_directories(${Readline_INCLUDE_DIRS})
    add_definitions(-DWITH_READLINE)
else()
    set(Readline_LIBRARY "")
endif()

# GMP
if(GMP)
    find_package(GMP)
endif()
if(GMP_FOUND)
    include_directories(${GMP_INCLUDE_DIRS})
    add_definitions(-DWITH_GMP)
else()
    set(GMP_LIBRARIES "")
endif()

#search the selected LP solver library
if(LPS STREQUAL "spx")
    find_package(SOPLEX REQUIRED CONFIG)
    if(LPSCHECK)
        find_package(CPLEX REQUIRED)
        add_definitions(-DWITH_LPSCHECK)
    endif()
elseif(LPS STREQUAL "cpx")
    find_package(CPLEX REQUIRED)
elseif(LPS STREQUAL "glop")
    find_package(GLOP REQUIRED)
elseif(LPS STREQUAL "grb")
    find_package(GUROBI REQUIRED)
elseif(LPS STREQUAL "clp")
    find_package(CLP REQUIRED)
elseif(LPS STREQUAL "xprs")
    find_package(XPRESS REQUIRED)
elseif(LPS STREQUAL "msk")
    find_package(MOSEK REQUIRED)
elseif(LPS STREQUAL "none")
    set(lpi lpi/lpi_none.c)
else()
    message(FATAL_ERROR "option LPS has wrong value")
endif()

#setup the proper lpi file for the selected LP solver
if(SOPLEX_FOUND)
    # SoPlex headers can be directly included
    include_directories(${SOPLEX_INCLUDE_DIRS})
    set(LPS_LIBRARIES ${SOPLEX_LIBRARIES})
    if(SHARED)
        set(LPS_PIC_LIBRARIES ${SOPLEX_PIC_LIBRARIES})
    else()
        set(LPS_PIC_LIBRARIES ${LPS_LIBRARIES})
    endif()
    if(LEGACY)
        set(lpi lpi/lpi_spx1.cpp)
    else()
        set(lpi lpi/lpi_spx2.cpp)
    endif()
endif()

if(CLP_FOUND)
    include_directories(${CLP_INCLUDE_DIRS})
    set(lpi lpi/lpi_clp.cpp)
    set(LPS_LIBRARIES ${CLP_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_LIBRARIES})
endif()

if(CPLEX_FOUND)
    include_directories(${CPLEX_INCLUDE_DIRS})
    # only use lpi_cpx.c if LPSCHECK is not enabled
    if(LPS STREQUAL "cpx")
        set(lpi lpi/lpi_cpx.c)
    endif()
    set(LPS_LIBRARIES ${LPS_LIBRARIES} ${CPLEX_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_PIC_LIBRARIES} ${LPS_LIBRARIES})
endif()

if(GLOP_FOUND)
    include_directories(${GLOP_INCLUDE_DIRS})
    set(lpi lpi/lpi_glop.cpp)
    set(LPS_LIBRARIES ${GLOP_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_LIBRARIES})
endif()

if(GUROBI_FOUND)
    include_directories(${GUROBI_INCLUDE_DIRS})
    set(lpi lpi/lpi_grb.c)
    set(LPS_LIBRARIES ${GUROBI_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_LIBRARIES})
endif()

if(XPRESS_FOUND)
    include_directories(${XPRESS_INCLUDE_DIRS})
    set(lpi lpi/lpi_xprs.c)
    set(LPS_LIBRARIES ${XPRESS_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_LIBRARIES})
endif()

if(MOSEK_FOUND)
    include_directories(${MOSEK_INCLUDE_DIRS})
    set(lpi lpi/lpi_msk.c)
    set(LPS_LIBRARIES ${MOSEK_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_LIBRARIES})
endif()

#
# enable coverage support
#
# it is very important to execute this prior to adding subdirectories because of
# the compiler flag changes
#
if( COVERAGE )
    include(CodeCoverage)
    APPEND_COVERAGE_COMPILER_FLAGS()
    set(COVERAGE_EXCLUDES '/usr*')
    SETUP_TARGET_FOR_COVERAGE(NAME coverage
                          EXECUTABLE bash -c "${CMAKE_CTEST_COMMAND} ${COVERAGE_CTEST_ARGS}"
                          DEPENDENCIES all_executables
                          )
endif()

# ZIMPL headers need to be copied to have the "zimpl/*.h" prefix
if(ZIMPL)
    find_package(ZIMPL CONFIG)
endif()
if(ZIMPL_FOUND)
    add_definitions(-DWITH_ZIMPL)
    include_directories(${ZIMPL_INCLUDE_DIRS})
    if(NOT SHARED)
        set(ZIMPL_PIC_LIBRARIES ${ZIMPL_LIBRARIES})
    endif()
else()
    set(ZIMPL_LIBRARIES "")
    set(ZIMPL_PIC_LIBRARIES "")
endif()

# IPOPT
if(IPOPT)
    find_package(IPOPT)
endif()
if(IPOPT_FOUND)
    include_directories(${IPOPT_INCLUDE_DIRS})
    #on debian IPOPT package needs this definition to work
    add_definitions( -DHAVE_CSTDDEF )
    set(NLPI_LIBRARIES ${IPOPT_LIBRARIES})
    set(nlpi nlpi/nlpi_ipopt.cpp)
else()
    set(nlpi nlpi/nlpi_ipopt_dummy.c)
endif()

# WORHP
if(WORHP)
    find_package(WORHP)
endif()
if(WORHP_FOUND)
    include_directories(${WORHP_INCLUDE_DIRS})
    set(nlpi ${nlpi} nlpi/nlpi_worhp.c)
    set(NLPI_LIBRARIES ${NLPI_LIBRARIES} ${WORHP_LIBRARIES})
else()
    set(nlpi ${nlpi} nlpi/nlpi_worhp_dummy.c)
endif()

# FilterSQP (with CMake, nlpi_filtersqp doesn't build anyway)
set(nlpi ${nlpi} nlpi/nlpi_filtersqp_dummy.c)

# run checks to figure out how the rounding mode can be set
include(CheckSymbolExists)
check_symbol_exists(FE_DOWNWARD "fenv.h" LINUX_ROUNDING)
if(LINUX_ROUNDING)
    add_definitions(-DROUNDING_FE)
else()
    check_symbol_exists(FP_RND_RM "float.h" OSF_ROUNDING)
    if(OSF_ROUNDING)
        add_definitions(-DROUNDING_FP)
    else()
        check_symbol_exists(RC_DOWN "float.h" MS_ROUNDING)
        if(MS_ROUNDING)
                add_definitions(-DROUNDING_MS)
        else()
                warning("cannot figure out how to set rounding mode")
        endif()
    endif()
endif()

find_package(Criterion)

# go to src/ and compile the code
add_subdirectory(src)

#
# we set the SCIP_DIR variable explicitly for the following examples/applications and unittests that depend on SCIP.
#
set(SCIP_DIR ${CMAKE_BINARY_DIR})
#
# add SCIP tests
#
add_subdirectory(check)

#
# add unit tests as a single target. Including tests will add the unit tests as single executables
#
add_custom_target(unittests)
add_subdirectory(tests EXCLUDE_FROM_ALL)
add_custom_target(all_executables DEPENDS scip unittests examples applications)
#
# add examples
#
# use sub directory bin/examples for executables of examples
#
set(OLD_CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples)
add_subdirectory(examples)

#
# add applications
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OLD_CMAKE_RUNTIME_OUTPUT_DIRECTORY}/applications)
add_subdirectory(applications)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OLD_CMAKE_RUNTIME_OUTPUT_DIRECTORY})

enable_testing()


include(FeatureSummary)
feature_summary(WHAT ALL)
