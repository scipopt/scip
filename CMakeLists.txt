cmake_minimum_required(VERSION 3.3)

# option() honors normal variables.
if(POLICY CMP0077)
   cmake_policy(SET CMP0077 NEW)
endif()

# MSVC runtime library flags are selected by an abstraction.
if(POLICY CMP0091)
   cmake_policy(SET CMP0091 NEW)
endif()

set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}")

project(SCIP)

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
   # for intel compilers we need to set the fp-model precise flag to disable incorrect floating point optimizations
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fp-model precise")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model precise")
   # to prevent symbols from svml exported twice we use the -shared-intel flag, see https://git.zib.de/integer/scip/issues/2872
   # this happens with the cplex static lib
   if(LPS MATCHES "cpx")
     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -shared-intel ")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -shared-intel ")
   endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-psabi")
   # require at least gcc 4.8
   if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
     message(WARNING "GCC version not supported, should be at least 4.8!")
   endif()
endif()


set(SCIP_VERSION_MAJOR 7)
set(SCIP_VERSION_MINOR 0)
set(SCIP_VERSION_PATCH 2)
set(SCIP_VERSION_SUB 4)
set(SCIP_VERSION_API 77)

get_directory_property(hasParent PARENT_DIRECTORY)
if(hasParent)
  set(SCIPOptSuite_VERSION_MAJOR ${SCIP_VERSION_MAJOR} PARENT_SCOPE)
  set(SCIPOptSuite_VERSION_MINOR ${SCIP_VERSION_MINOR} PARENT_SCOPE)
  set(SCIPOptSuite_VERSION_PATCH ${SCIP_VERSION_PATCH} PARENT_SCOPE)
  set(SCIP_DIR ${PROJECT_BINARY_DIR} PARENT_SCOPE)
endif()

set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "${SCIP_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${SCIP_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${SCIP_VERSION_PATCH}")
set(CPACK_PACKAGE_VENDOR "Zuse Institute Berlin")
set(CPACK_PACKAGE_CONTACT "http://scipopt.org")
include(CPack)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

if(SCIPOptSuite_BINARY_DIR)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SCIPOptSuite_BINARY_DIR}/bin)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${SCIPOptSuite_BINARY_DIR}/lib)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${SCIPOptSuite_BINARY_DIR}/lib)
  #set(SCIP_DIR ${SCIPOptSuite_BINARY_DIR} PARENT_SCOPE)
endif()

# path to e.g. findGMP module
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

option(SHARED "Build shared libraries" ON)
set(BUILD_SHARED_LIBS ${SHARED})
message(STATUS "Build shared libraries: " ${SHARED})

# make 'Release' the default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

option(ZLIB "should zlib be linked" ON)
option(READLINE "should readline be linked" ON)
option(GMP "should gmp be linked" ON)
option(STATIC_GMP "prefer static gmp library" OFF)
option(PAPILO "should papilo library be linked" ON)
option(ZIMPL "should zimpl be linked" ON)
option(IPOPT "should ipopt be linked" ON)
option(WORHP "should worhp be linked" OFF)
# for backward compability
option(PARASCIP "should scip be compiled thread safe (deprecated, use THREADSAFE)" OFF)
option(THREADSAFE "should scip be compiled thread safe" ON)
option(LPSCHECK "double check SoPlex results with CPLEX" OFF)
option(LEGACY "build SCIP/SoPlex in legacy mode (spx1 interface)" OFF)
option(NOBLKBUFMEM "should block and buffer memory be disabled" OFF)
option(NOBLKMEM "should block memory be disabled" OFF)
option(NOBUFMEM "should buffer memory be disabled" OFF)
option(DEBUGSOL "should the debug solution mechanism be enabled" OFF)
option(SANITIZE_ADDRESS "should the address sanitizer be enabled in debug mode if available" OFF)
option(SANITIZE_MEMORY "should the memory sanitizer be enabled in debug mode if available" OFF)
option(SANITIZE_UNDEFINED "should the undefined behavior sanitizer be enabled in debug mode if available" ON)
option(SANITIZE_THREAD "should the thread sanitizer be enabled in debug mode if available" OFF)
option(COVERAGE "enable coverage support" OFF)
SET(COVERAGE_CTEST_ARGS "" CACHE STRING "additional ctest arguments for coverage")
option(MT "use static runtime libraries for Visual Studio compiler" OFF)
option(CXXONLY "use a c++ compiler for all source files" OFF)
option(AMPL "Enable AMPL support" OFF)

set(TPI none CACHE STRING "options for thread support library")  #create the variable
set_property(CACHE TPI PROPERTY STRINGS none tny omp)  #define list of values GUI will offer for the variable

set(EXPRINT cppad CACHE STRING "options for expression interpreter")  #create the variable
set_property(CACHE EXPRINT PROPERTY STRINGS none cppad )  #define list of values GUI will offer for the variable

set(LPS spx CACHE STRING "options for LP solver")  #create the variable
set_property(CACHE LPS PROPERTY STRINGS spx cpx grb xprs clp glop msk qso none )  #define list of values GUI will offer for the variable

set(SYM bliss CACHE STRING "options for symmetry computation")  #create the variable
set_property(CACHE SYM PROPERTY STRINGS bliss none )  #define list of values GUI will offer for the variable

#set options for memory management
if( NOBLKBUFMEM )
   set(NOBLKMEM ON FORCE)
   set(NOBUFMEM ON FORCE)
endif()

if(NOBLKMEM)
   set(BMS_NOBLOCKMEM on)
endif()

if(NOBUFMEM)
   set(SCIP_NOBUFFERMEM on)
endif()

if(DEBUGSOL)
   set(WITH_DEBUG_SOLUTION on)
endif()

#set the correct rpath for OS X
set(CMAKE_MACOSX_RPATH ON)

#set defines for Windows
if(WIN32)
    set(SCIP_NO_SIGACTION on)
    set(SCIP_NO_STRTOK_R on)
endif()
if(MSVC)
#    add_definitions(/W4)
    add_definitions(/wd4100)
    add_definitions(/wd4244)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

# Visual Studio compiler with static runtime libraries
if(MSVC AND MT)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
endif()

#set expression interpreter file that should be used
if(EXPRINT STREQUAL "cppad")
    set(exprinterpret nlpi/exprinterpret_cppad.cpp)
elseif(EXPRINT STREQUAL "none")
    set(exprinterpret nlpi/exprinterpret_none.c)
else()
    message(FATAL_ERROR "EXPRINT option has wrong value")
endif()

# if tpi is not none force THREADSAFE to on
if(NOT (TPI STREQUAL "none"))
   set(THREADSAFE ON FORCE)
endif()

#set sources files, libraries and defines for tpi
if(TPI STREQUAL "none")
    set(tpisources tpi/tpi_none.c)
    set(THREAD_LIBRARIES "")
    set(TPI_NONE on)
elseif(TPI STREQUAL "tny")
    set(TPI_TNYC on)
    set(tpisources tpi/tpi_tnycthrd.c tinycthread/tinycthread.c)
    find_package(Threads REQUIRED)
    set(THREAD_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
elseif(TPI STREQUAL "omp")
    set(TPI_OMP on)
    set(tpisources tpi/tpi_openmp.c)
    find_package(OpenMP REQUIRED)
    set(THREAD_LIBRARIES "")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(FATAL_ERROR "TPI option has wrong value")
endif()

# PARASCIP is deprecated, use THREADSAFE
if(PARASCIP)
  set(THREADSAFE ON FORCE)
endif()

#set SCIP_THREADSAFE define
if(THREADSAFE)
    set(SCIP_THREADSAFE on)
endif()

set(NEWLINE "\\\\n")

# create a target for updating the current git hash
file(WRITE ${PROJECT_BINARY_DIR}/scip_update_githash.cmake "
find_program(GIT git)
if(EXISTS \${DST})
   file(STRINGS \${DST} GITHASH_OLD)
   string(REGEX REPLACE \"#define SCIP_GITHASH \\\"(.*)\\\"\" \"\\\\1\" GITHASH_OLD \${GITHASH_OLD})
endif()
if((GIT) AND (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git))
   execute_process(
      COMMAND \${GIT} describe --always --dirty
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE GITHASH OUTPUT_STRIP_TRAILING_WHITESPACE)
   string(REGEX REPLACE \"^.*-g\" \"\" GITHASH \${GITHASH})
   if(NOT \${GITHASH} STREQUAL \"\${GITHASH_OLD}\")
      file(WRITE \${DST} \"#define SCIP_GITHASH \\\"\${GITHASH}\\\"\n\")
   endif()
else()
   set(GITHASH \${GITHASH_OLD})
   if(NOT GITHASH)
      message(STATUS \"Compiling without git information\")
      set(GITHASH \"NoGitInfo\")
   endif()
   file(WRITE \${DST} \"#define SCIP_GITHASH \\\"\${GITHASH}\\\"\n\")
endif()
message(STATUS \"Git hash: \" \${GITHASH})
")
add_custom_target(scip_update_githash
   COMMAND ${CMAKE_COMMAND}
      -DDST=${PROJECT_SOURCE_DIR}/src/scip/githash.c
      -P ${PROJECT_BINARY_DIR}/scip_update_githash.cmake)

set(WITH_SCIPDEF on)

# use C++14 standard
set(CMAKE_CXX_STANDARD 14)

# set function visibility default to hidden
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)

# ZLIB
if(ZLIB)
    message(STATUS "Finding ZLIB")
    if(NOT TARGET ZLIB::ZLIB)
        find_package(ZLIB REQUIRED)
    endif()
    message(STATUS "Finding ZLIB - found")
    set(SCIP_WITH_ZLIB ON)
else()
    message(STATUS "Support ZLIB: OFF")
    set(SCIP_WITH_ZLIB OFF)
endif()

# Readline
if(READLINE)
    message(STATUS "Finding Readline")
    find_package(Readline REQUIRED)
    message(STATUS "Finding Readline - found")
    include_directories(${Readline_INCLUDE_DIRS})
    set(SCIP_WITH_READLINE ON)
else()
    message(STATUS "Support Readline: OFF")
    set(Readline_LIBRARY "")
    set(SCIP_WITH_READLINE OFF)
endif()

# GMP
if(GMP)
    message(STATUS "Finding GMP")
    find_package(GMP REQUIRED)
    message(STATUS "Finding GMP - found")
    include_directories(${GMP_INCLUDE_DIRS})
    set(SCIP_WITH_GMP ON)
else()
    message(STATUS "Support GMP: OFF")
    set(GMP_LIBRARIES "")
    set(SCIP_WITH_GMP OFF)
endif()

# look for presolvelib
if(PAPILO)
    message(STATUS "Finding PAPILO")
    find_package(PAPILO CONFIG REQUIRED)
    message(STATUS "Finding PAPILO - found")
    set(SCIP_WITH_PAPILO ON)
else()
    message(STATUS "Support PAPILO: OFF")
    set(PAPILO_IMPORTED_TARGETS "")
    set(SCIP_WITH_PAPILO OFF)
endif()

#search the selected LP solver library
message(STATUS "Finding Solver \"${LPS}\"")
if(LPS STREQUAL "spx")
    message(STATUS "Finding Soplex")
    find_package(SOPLEX REQUIRED CONFIG)
    if(LPSCHECK)
        find_package(CPLEX REQUIRED)
        set(SCIP_WITH_LPSCHECK on)
    endif()
elseif(LPS STREQUAL "cpx")
    find_package(CPLEX REQUIRED)
elseif(LPS STREQUAL "glop")
    find_package(GLOP REQUIRED)
elseif(LPS STREQUAL "grb")
    find_package(GUROBI REQUIRED)
elseif(LPS STREQUAL "qso")
    find_package(QSO REQUIRED)
    set(BUILD_SHARED_LIBS off)
    message(STATUS "Turning off shared libraries.")
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
    message(STATUS "Finding SOPLEX - found")
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
else()
    message(STATUS "Support SOPLEX: OFF")
endif()

if(CLP_FOUND)
    message(STATUS "Finding CLP - found")
    include_directories(${CLP_INCLUDE_DIRS})
    set(lpi lpi/lpi_clp.cpp)
    set(LPS_LIBRARIES ${CLP_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_LIBRARIES})
else()
    message(STATUS "Support CLP: OFF")
endif()

if(CPLEX_FOUND)
    message(STATUS "Finding CPLEX - found")
    include_directories(${CPLEX_INCLUDE_DIRS})
    # only use lpi_cpx.c if LPSCHECK is not enabled
    if(LPS STREQUAL "cpx")
        set(lpi lpi/lpi_cpx.c)
    endif()
    set(LPS_LIBRARIES ${LPS_LIBRARIES} ${CPLEX_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_PIC_LIBRARIES} ${LPS_LIBRARIES})
else()
    message(STATUS "Support CPLEX: OFF")
endif()

if(GLOP_FOUND)
    message(STATUS "Finding GLOP - found")
    include_directories(${GLOP_INCLUDE_DIRS})
    set(lpi lpi/lpi_glop.cpp)
    set(LPS_LIBRARIES ${GLOP_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_LIBRARIES})
    set(CMAKE_CXX_STANDARD 17)
else()
    message(STATUS "Support GLOP: OFF")
endif()

if(GUROBI_FOUND)
    message(STATUS "Finding GUROBI - found")
    include_directories(${GUROBI_INCLUDE_DIRS})
    set(lpi lpi/lpi_grb.c)
    set(LPS_LIBRARIES ${GUROBI_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_LIBRARIES})
else()
    message(STATUS "Support GUROBI: OFF")
endif()

if(XPRESS_FOUND)
    message(STATUS "Finding XPRESS - found")
    include_directories(${XPRESS_INCLUDE_DIRS})
    set(lpi lpi/lpi_xprs.c)
    set(LPS_LIBRARIES ${XPRESS_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_LIBRARIES})
else()
    message(STATUS "Support XPRESS: OFF")
endif()

if(MOSEK_FOUND)
    message(STATUS "Finding MOSEK - found")
    include_directories(${MOSEK_INCLUDE_DIRS})
    set(lpi lpi/lpi_msk.c)
    set(LPS_LIBRARIES ${MOSEK_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_LIBRARIES})
else()
    message(STATUS "Support MOSEK: OFF")
endif()

if(QSO_FOUND)
    message(STATUS "Finding QSO - found")
    include_directories(${QSO_INCLUDE_DIRS})
    set(lpi lpi/lpi_qso.c)
    set(LPS_LIBRARIES ${QSO_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_LIBRARIES})
else()
    message(STATUS "Support QSO: OFF")
endif()

#search the selected symmetry computation program
message(STATUS "Finding symmetry computation program \"${SYM}\"")
if(SYM STREQUAL "bliss")
    message(STATUS "Finding BLISS")
    find_package(BLISS)
    if(BLISS_FOUND)
        message(STATUS "Finding BLISS - found")
        include_directories(${BLISS_INCLUDE_DIRS})
        set(sym symmetry/compute_symmetry_bliss.cpp)
        set(SYM_LIBRARIES ${BLISS_LIBRARIES})
        set(SYM_PIC_LIBRARIES ${BLISS_LIBRARIES})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BLISS_DEFINITIONS}")
    else()
        message(STATUS "Finding BLISS - not found")
        set(sym symmetry/compute_symmetry_none.cpp)
    endif()
elseif(SYM STREQUAL "none")
    message(STATUS "Support SYM: OFF")
    set(sym symmetry/compute_symmetry_none.cpp)
else()
    message(FATAL_ERROR "option SYM has wrong value")
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

    #
    # create a CMake script file that is executed to run the coverage test.
    # With a script file, the return code of the tests is simply ignored,
    # and a coverage report is generated even if some tests fail currently.
    #
    file(WRITE ${PROJECT_BINARY_DIR}/RunCoverage.cmake "execute_process(COMMAND ctest ${CMAKE_CTEST_COMMAND} ${COVERAGE_CTEST_ARGS})")

    #
    # setup the coverage target to execute the RunCoverage script
    #
    SETUP_TARGET_FOR_COVERAGE(NAME coverage
                          EXECUTABLE ${CMAKE_COMMAND} -P ${PROJECT_BINARY_DIR}/RunCoverage.cmake
                          DEPENDENCIES all_executables
                          )
else()
  message(STATUS "CodeCoverage: OFF")
endif()

# ZIMPL headers need to be copied to have the "zimpl/*.h" prefix
if(ZIMPL)
    message(STATUS "Finding ZIMPL")
    find_package(ZIMPL CONFIG REQUIRED)
    message(STATUS "Finding ZIMPL - found")
    set(SCIP_WITH_ZIMPL ON)
    include_directories(${ZIMPL_INCLUDE_DIRS})
    if(NOT SHARED)
        set(ZIMPL_PIC_LIBRARIES ${ZIMPL_LIBRARIES})
    endif()
else()
    message(STATUS "Support ZIMPL: OFF")
    set(ZIMPL_LIBRARIES "")
    set(ZIMPL_PIC_LIBRARIES "")
    set(SCIP_WITH_ZIMPL OFF)
endif()

# IPOPT
if(IPOPT)
    message(STATUS "Finding IPOPT")
    find_package(IPOPT 3.12.0 REQUIRED)
    message(STATUS "Finding IPOPT - found")
    include_directories(${IPOPT_INCLUDE_DIRS})
    #on debian IPOPT package needs this definition to work
    set(HAVE_CSTDDEF on)
    set(NLPI_LIBRARIES ${IPOPT_LIBRARIES})
    set(nlpi nlpi/nlpi_ipopt.cpp)
else()
    message(STATUS "Support IPOPT: OFF")
    set(nlpi nlpi/nlpi_ipopt_dummy.c)
endif()

# WORHP
if(WORHP)
    message(STATUS "Finding WORHP")
    find_package(WORHP REQUIRED)
    message(STATUS "Finding WORHP - found")
    include_directories(${WORHP_INCLUDE_DIRS})
    set(nlpi ${nlpi} nlpi/nlpi_worhp.c)
    set(NLPI_LIBRARIES ${NLPI_LIBRARIES} ${WORHP_LIBRARIES})
else()
    message(STATUS "Support WORHP: OFF")
    set(nlpi ${nlpi} nlpi/nlpi_worhp_dummy.c)
endif()

# FilterSQP (with CMake, nlpi_filtersqp doesn't build anyway)
set(nlpi ${nlpi} nlpi/nlpi_filtersqp_dummy.c)

# run checks to figure out how the rounding mode can be set
include(CheckSymbolExists)
check_symbol_exists(FE_DOWNWARD "fenv.h" LINUX_ROUNDING)
if(LINUX_ROUNDING)
    set(SCIP_ROUNDING_FE on)
else()
    check_symbol_exists(FP_RND_RM "float.h" OSF_ROUNDING)
    if(OSF_ROUNDING)
        set(SCIP_ROUNDING_FP on)
    else()
        check_symbol_exists(RC_DOWN "float.h" MS_ROUNDING)
        if(MS_ROUNDING)
            set(SCIP_ROUNDING_MS on)
        else()
            message(WARNING "cannot figure out how to set rounding mode")
        endif()
    endif()
endif()

message(STATUS "Finding CRITERION")
find_package(CRITERION)
if(CRITERION_FOUND)
    message(STATUS "Finding CRITERION - found")
else()
    message(STATUS "Finding CRITERION - not found")
endif()


# export compilation settings to header file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/scip/config.h.in ${PROJECT_BINARY_DIR}/scip/config.h @ONLY)

# go to src/ and compile the code
add_subdirectory(src)

# store directory
set(OLD_CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

include(CTest)
if(BUILD_TESTING)
    #
    # add SCIP tests
    #
    add_subdirectory(check)

    #
    # add unit tests as a single target. Including tests will add the unit tests as single executables
    #
    add_custom_target(unittests)
    add_subdirectory(tests EXCLUDE_FROM_ALL)
    add_subdirectory(doc EXCLUDE_FROM_ALL)

    #
    # add examples
    #
    # use sub directory bin/examples for executables of examples
    #
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OLD_CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples)
    add_subdirectory(examples)

    #
    # add applications
    #
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OLD_CMAKE_RUNTIME_OUTPUT_DIRECTORY}/applications)
    add_subdirectory(applications)

    add_custom_target(all_executables DEPENDS scip unittests examples applications)
endif()

if(AMPL)
    #
    # add AMPL interface
    #
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OLD_CMAKE_RUNTIME_OUTPUT_DIRECTORY}/interfaces/ampl)
    add_subdirectory(interfaces/ampl EXCLUDE_FROM_ALL)
else()
    message(STATUS "Support AMPL: OFF")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OLD_CMAKE_RUNTIME_OUTPUT_DIRECTORY})

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    include(FeatureSummary)
    feature_summary(WHAT ALL)
endif()
