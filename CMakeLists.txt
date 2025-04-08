cmake_minimum_required(VERSION 3.11)

# option() honors normal variables.
if(POLICY CMP0077)
   cmake_policy(SET CMP0077 NEW)
endif()

# MSVC runtime library flags are selected by an abstraction.
if(POLICY CMP0091)
   cmake_policy(SET CMP0091 NEW)
endif()

# Avoid warnings about removed FindBoost.cmake until our minimum required version is >= 3.30.
if(POLICY CMP0167)
   cmake_policy(SET CMP0167 NEW)
endif()

set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}")

set(SCIP_VERSION_MAJOR 10)
set(SCIP_VERSION_MINOR 0)
set(SCIP_VERSION_PATCH 0)
set(SCIP_VERSION_API 141)

project(SCIP
    VERSION ${SCIP_VERSION_MAJOR}.${SCIP_VERSION_MINOR}.${SCIP_VERSION_PATCH}
    LANGUAGES CXX C
    DESCRIPTION "SCIP is currently one of the fastest non-commercial solvers for mixed integer programming (MIP) and mixed integer nonlinear programming (MINLP)"
    )

get_directory_property(hasParent PARENT_DIRECTORY)
if(hasParent)
  set(SCIPOptSuite_VERSION_MAJOR ${SCIP_VERSION_MAJOR} PARENT_SCOPE)
  set(SCIPOptSuite_VERSION_MINOR ${SCIP_VERSION_MINOR} PARENT_SCOPE)
  set(SCIPOptSuite_VERSION_PATCH ${SCIP_VERSION_PATCH} PARENT_SCOPE)
  set(SCIP_DIR ${PROJECT_BINARY_DIR} PARENT_SCOPE)
endif()

set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
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
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

option(SHARED "Build shared libraries" ON)
set(BUILD_SHARED_LIBS ${SHARED})
message(STATUS "Build shared libraries: " ${SHARED})

# make 'Release' the default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

option(AUTOBUILD "find and use dependencies as available" OFF)
set(AUTOBUILD_MSG "\nIf you have troubles configuring, you can consider setting AUTOBUILD=ON to try and find optional packages as available.")
option(ZLIB "should zlib be linked" ON)
option(READLINE "should readline be linked" ON)
option(GMP "should gmp be linked" ON)
option(STATIC_GMP "prefer static gmp library" OFF)
option(PAPILO "should papilo library be linked" ON)
option(ZIMPL "should zimpl be linked" ON)
option(AMPL "should ampl interface be compiled" ON)
option(IPOPT "should ipopt be linked" ON)
option(LAPACK "should lapack be linked" OFF)
option(WORHP "should worhp be linked" OFF)
# for backward compatibility
option(PARASCIP "should scip be compiled thread safe (deprecated, use THREADSAFE)" OFF)
option(THREADSAFE "should scip be compiled thread safe" ON)
option(LPSCHECK "double check SoPlex results with CPLEX" OFF)
option(LEGACY "build SCIP/SoPlex in legacy mode (spx1 interface)" OFF)
option(NOBLKBUFMEM "should block and buffer memory be disabled" OFF)
option(NOBLKMEM "should block memory be disabled" OFF)
option(NOBUFMEM "should buffer memory be disabled" OFF)
option(DEBUGSOL "should the debug solution mechanism be enabled" OFF)
option(SANITIZE "should sanitizers be enabled if available" OFF)
option(COVERAGE "enable coverage support" OFF)
SET(COVERAGE_CTEST_ARGS "" CACHE STRING "additional ctest arguments for coverage")
option(MT "use static runtime libraries for Visual Studio compiler" OFF)
option(CXXONLY "use a c++ compiler for all source files" OFF)

set(MPFR AUTO CACHE STRING "should MPFR be linked")
set_property(CACHE MPFR PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER ${MPFR} MPFR)

set(EXACTSOLVE AUTO CACHE STRING "compile with support for exact rational solve (requires GMP, MPFR and Boost)")
set_property(CACHE EXACTSOLVE PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER ${EXACTSOLVE} EXACTSOLVE_UPPER) # Also allow "auto"

# Initialize SCIP_WITH_EXACTSOLVE.
if(EXACTSOLVE)
    set(SCIP_WITH_EXACTSOLVE ON) # We disable it later if dependencies are missing.
else()
    set(SCIP_WITH_EXACTSOLVE OFF)
endif()

# Initialize SCIP_WITH_MPFR.
if(MPFR)
    set(SCIP_WITH_MPFR ON) # We disable it later if exact solving mode is disabled.
else()
    set(SCIP_WITH_MPFR OFF)
    if(SCIP_WITH_EXACTSOLVE)
        message(WARNING "EXACTSOLVE is disabled because this requires MPFR.")
        set(SCIP_WITH_EXACTSOLVE OFF)
    endif()
endif()

# Initialize SCIP_WITH_GMP.
if(GMP)
    set(SCIP_WITH_GMP ON)
else()
    set(SCIP_WITH_GMP OFF)
    if(SCIP_WITH_EXACTSOLVE)
        message(WARNING "EXACTSOLVE is disabled because this requires GMP. ${AUTOBUILD_MSG}")
        set(SCIP_WITH_EXACTSOLVE OFF)
    endif()
endif()

set(TPI tny CACHE STRING "options for thread support library")  #create the variable
set_property(CACHE TPI PROPERTY STRINGS none tny omp)  #define list of values GUI will offer for the variable

set(EXPRINT cppad CACHE STRING "options for expression interpreter")  #create the variable
set_property(CACHE EXPRINT PROPERTY STRINGS none cppad )  #define list of values GUI will offer for the variable

set(LPS spx CACHE STRING "options for LP solver")  #create the variable
set_property(CACHE LPS PROPERTY STRINGS spx cpx grb xprs clp glop msk qso highs none)  #define list of values GUI will offer for the variable
set(LPSEXACT AUTO CACHE STRING "options for exact LP solver") # Create the variable
set_property(CACHE LPSEXACT PROPERTY STRINGS spx qsoex none AUTO) # Define list of values GUI will offer for the variable.
string(TOUPPER ${LPSEXACT} LPSEXACT_UPPER) # To also allow "auto".

set(SYM snauty CACHE STRING "options for symmetry computation")  #create the variable
set_property(CACHE SYM PROPERTY STRINGS bliss sbliss nauty snauty dejavu none )  #define list of values GUI will offer for the variable

if(NOT (CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" OR CMAKE_BUILD_TYPE STREQUAL "MinSizeRel"))
    message(FATAL_ERROR "You have not selected a valid CMAKE_BUILD_TYPE. Please select 'Release' or 'Debug'.")
endif()

if (NOT (LPSEXACT_UPPER STREQUAL "AUTO") AND NOT (LPSEXACT_UPPER STREQUAL "NONE") AND NOT (LPSEXACT_UPPER STREQUAL "SPX") AND NOT (LPSEXACT_UPPER STREQUAL "QSOEX"))
    message(FATAL_ERROR "Invalid value \"${LPSEXACT}\" for LPSEXACT. Please select 'auto', 'none', 'spx' or 'qsoex'.")
endif()

if(AUTOBUILD)
    message(WARNING "You are using the automatic build option AUTOBUILD=ON to find and use dependencies automatically.")
else()
    message(WARNING "You are not using the automatic build option. ${AUTOBUILD_MSG}")
endif() # AUTOBUILD

# Use C++14 standard; may be increased further down.
set(CMAKE_CXX_STANDARD 14)

#search the selected symmetry computation program
message(STATUS "Finding symmetry computation program \"${SYM}\"")

if(SYM STREQUAL "bliss" OR SYM STREQUAL "sbliss")
    if(SYM STREQUAL "bliss")
        message(STATUS "Support SYM: bliss")
        set(sym symmetry/compute_symmetry_bliss.cpp)
    else()
        message(STATUS "Support SYM: sassy+bliss")
        set(sym symmetry/compute_symmetry_sassy_bliss.cpp symmetry/build_dejavu_graph.cpp)
        # sassy needs C++-14
        set(CMAKE_CXX_STANDARD 14)
    endif()

    find_package(Bliss CONFIG HINTS ${BLISS_DIR})
    if(Bliss_FOUND)
        set(SYM_LIBRARIES Bliss::libbliss)
        set(SYM_PIC_LIBRARIES Bliss::libbliss)
        message(STATUS "Found Bliss: ${Bliss_DIR}")
    else()
        message(FATAL_ERROR "Bliss not found, try installing https://github.com/scipopt/bliss and specifying BLISS_DIR.")
    endif()
elseif(SYM STREQUAL "nauty" OR SYM STREQUAL "snauty")
    if(NAUTY_DIR)
        include_directories(${NAUTY_DIR})
        set(SYM_LIBRARIES ${NAUTY_DIR}/nauty/nauty.a)
        set(SYM_PIC_LIBRARIES ${NAUTY_DIR}/nauty/nauty.a)
    else()
        include_directories(src/nauty)
        set(sym nauty/nauty.c nauty/nautil.c nauty/nausparse.c nauty/schreier.c nauty/naurng.c)
    endif()
    if(SYM STREQUAL "nauty")
        message(STATUS "Support SYM: nauty")
        set(sym ${sym} symmetry/compute_symmetry_nauty.c)
    elseif(SYM STREQUAL "snauty")
        message(STATUS "Support SYM: sassy+nauty")
        set(sym ${sym} symmetry/compute_symmetry_sassy_nauty.cpp symmetry/build_dejavu_graph.cpp)
        # sassy needs C++-14
        set(CMAKE_CXX_STANDARD 14)
    endif()
elseif(SYM STREQUAL "dejavu")
    include_directories(src/dejavu)
    message(STATUS "Support SYM: dejavu")
    set(sym ${sym} symmetry/compute_symmetry_dejavu.cpp symmetry/build_dejavu_graph.cpp)
    # dejavu needs C++-14
    set(CMAKE_CXX_STANDARD 14)
elseif(SYM STREQUAL "none")
    message(STATUS "Support SYM: OFF")
    set(sym symmetry/compute_symmetry_none.cpp)
else()
    message(FATAL_ERROR "option SYM has wrong value")
endif()

# if changing these flags, also update GCCWARN/GXXWARN in make/make.project
set(ADD_C_AND_CXX_FLAGS -pedantic -Wall -W -Wpointer-arith -Wcast-align -Wwrite-strings -Wshadow  -Wredundant-decls -Wdisabled-optimization -Wno-long-long -Wno-unknown-pragmas -Wno-unused-parameter)
set(ADD_C_FLAGS -Wsign-compare -Wstrict-prototypes -Wmissing-declarations -Wmissing-prototypes -Wno-override-init -Wno-uninitialized -Wno-maybe-uninitialized)
set(ADD_CXX_FLAGS -Wnon-virtual-dtor -Wreorder -Woverloaded-virtual -Wsign-promo -Wsynth -Wcast-qual -Wno-strict-overflow -Wno-psabi)

# disable fused floating point contraction to enhance reproducibility across compilers and architectures
if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_compile_options(/fp:precise)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_compile_options(-fp-model=precise)
    # to prevent symbols from svml exported twice we use the -shared-intel flag, see https://git.zib.de/integer/scip/issues/2872
    # this happens with the cplex static lib
    if(LPS MATCHES "cpx")
        add_compile_options(-shared-intel) # for c and cxx
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # require at least gcc 4.8
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
        message(WARNING "GCC version not supported, should be at least 4.8!")
    endif()
    add_compile_options(-ffp-contract=off ${ADD_C_AND_CXX_FLAGS})
    add_compile_options("$<$<COMPILE_LANGUAGE:C>:${ADD_C_FLAGS}>$<$<COMPILE_LANGUAGE:CXX>:${ADD_CXX_FLAGS}>")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # clang specific suppressions
    set(ADD_C_AND_CXX_FLAGS ${ADD_C_AND_CXX_FLAGS} -Wno-unknown-warning-option)
    set(ADD_C_FLAGS ${ADD_C_FLAGS} -Wno-initializer-overrides)
    add_compile_options(-ffp-contract=off ${ADD_C_AND_CXX_FLAGS})
    add_compile_options("$<$<COMPILE_LANGUAGE:C>:${ADD_C_FLAGS}>$<$<COMPILE_LANGUAGE:CXX>:${ADD_CXX_FLAGS}>")
endif()

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
    add_definitions(/Zc:__cplusplus)
    add_definitions(/EHsc)
endif()

# Visual Studio compiler with static runtime libraries
if(MSVC AND MT)
    add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>:/MTd>$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELEASE>>:/MT>$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:DEBUG>>:/MTd>$<$<AND:$<COMPILE_LANGUAGE:C>,$<CONFIG:RELEASE>>:/MT>")
endif()

#set expression interpreter file that should be used
if(EXPRINT STREQUAL "cppad")
    set(exprinterpret scip/exprinterpret_cppad.cpp)
elseif(EXPRINT STREQUAL "none")
    set(exprinterpret scip/exprinterpret_none.c)
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
    set(TPI_TNY on)
    set(tpisources tpi/tpi_tnycthrd.c tinycthread/tinycthread.c)
    find_package(Threads REQUIRED)
    set(THREAD_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
elseif(TPI STREQUAL "omp")
    set(TPI_OMP on)
    set(tpisources tpi/tpi_openmp.c)
    find_package(OpenMP REQUIRED)
    set(THREAD_LIBRARIES "")
    add_compile_options("$<$<COMPILE_LANGUAGE:C>:${OpenMP_C_FLAGS}>$<$<COMPILE_LANGUAGE:CXX>:${OpenMP_CXX_FLAGS}>")
    add_link_options("${OpenMP_CXX_FLAGS}")
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
   string(REGEX REPLACE \"#define SCIP_GITHASH \\\"(.*)\\\"\" \"\\\\1\" GITHASH_OLD \"\${GITHASH_OLD}\")
endif()
if((GIT) AND (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git))
   execute_process(
      COMMAND \${GIT} describe --always --dirty
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE GITHASH OUTPUT_STRIP_TRAILING_WHITESPACE)
   string(REGEX REPLACE \"^.*-g\" \"\" GITHASH \"\${GITHASH}\")
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

# ZLIB
if(ZLIB OR AUTOBUILD)
    message(STATUS "Finding ZLIB")
    if(NOT TARGET ZLIB::ZLIB)
        find_package(ZLIB)
        if((NOT ZLIB_FOUND) AND (NOT AUTOBUILD))
            message(FATAL_ERROR "ZLIB not found, try specifying ZLIB_DIR. ${AUTOBUILD_MSG}")
        endif()
    endif()
endif()
if(ZLIB_FOUND)
    message(STATUS "Finding ZLIB - found")
    set(SCIP_WITH_ZLIB ON)
else()
    message(STATUS "Support ZLIB: OFF")
    set(SCIP_WITH_ZLIB OFF)
endif()

# Readline
if(READLINE OR AUTOBUILD)
    message(STATUS "Finding Readline and Curses")
    find_package(Readline)
    if((NOT Readline_FOUND) AND (NOT AUTOBUILD))
        message(FATAL_ERROR "READLINE not found, try specifying READLINE_DIR. ${AUTOBUILD_MSG}")
    endif()
    find_package(Curses)
endif()
if(Readline_FOUND AND CURSES_FOUND)
    message(STATUS "Finding Readline - found")
    include_directories(${Readline_INCLUDE_DIRS})
    set(SCIP_WITH_READLINE ON)
    set(Readline_LIBRARY ${Readline_LIBRARY} ${CURSES_LIBRARIES})
elseif(Readline_FOUND)
    message(STATUS "Finding Readline - found, but without Curses library.")
    include_directories(${Readline_INCLUDE_DIRS})
    set(SCIP_WITH_READLINE ON)
else()
    message(STATUS "Support Readline: OFF")
    set(Readline_LIBRARY "")
    set(SCIP_WITH_READLINE OFF)
endif()

# GMP
if(SCIP_WITH_GMP OR AUTOBUILD)
    message(STATUS "Finding GMP")
    find_package(GMP)
    if(NOT GMP_FOUND)
        set(SCIP_WITH_EXACTSOLVE OFF)
        set(SCIP_WITH_GMP OFF)
        if(NOT AUTOBUILD) # GMP was explicitly requested.
            message(FATAL_ERROR "GMP not found, try specifying GMP_DIR. ${AUTOBUILD_MSG}")
        elseif(EXACTSOLVE AND NOT (EXACTSOLVE_UPPER STREQUAL "AUTO")) # EXACTSOLVE was explicitly requested.
          message(FATAL_ERROR "GMP not found, which is required for EXACTSOLVE; try specifying GMP_DIR. ${AUTOBUILD_MSG}")
        elseif(EXACTSOLVE) # Both are automatic.
          message(WARNING "GMP not found, which is required for EXACTSOLVE; try specifying GMP_DIR. ${AUTOBUILD_MSG}")
        else() # EXACTSOLVE disabled. TODO: check if there is any reason for linking to it at all.
          message(STATUS "GMP not found; try specifying GMP_DIR. ${AUTOBUILD_MSG}")
        endif()
    endif()
else()
    message(STATUS "Skipping GMP.")
endif()
if(SCIP_WITH_GMP)
    # TODO: Avoid global usage.
    include_directories(${GMP_INCLUDE_DIRS})
else()
    set(GMP_LIBRARIES "")
endif()

# MPFR
if(SCIP_WITH_EXACTSOLVE AND SCIP_WITH_MPFR)
    message(STATUS "Finding MPFR")
    find_package(MPFR)
    if(NOT MPFR_FOUND)
        set(SCIP_WITH_MPFR OFF)
        set(SCIP_WITH_EXACTSOLVE OFF)
        if(MPFR AND NOT (MPFR STREQUAL "AUTO")) # MPFR was explicitly requested.
            message(FATAL_ERROR "MPFR not found; try specifying MPFR_DIR.")
        elseif(EXACTSOLVE AND NOT (EXACTSOLVE_UPPER STREQUAL "AUTO")) # EXACTSOLVE was explicitly requested.
            message(FATAL_ERROR "MPFR not found, which is required for EXACTSOLVE; try specifying MPFR_DIR.")
        else() # Both are automatic.
            message(WARNING "MPFR not found, which is required for EXACTSOLVE; try specifying MPFR_DIR.")
        endif()
    endif()
else()
   message(STATUS "Skipping MPFR.")
endif()

# look for presolvelib
if(PAPILO OR AUTOBUILD)
    message(STATUS "Finding PAPILO")
    find_package(PAPILO CONFIG HINTS ${PAPILO_DIR})
    if((NOT PAPILO_FOUND) AND (NOT AUTOBUILD))
        message(FATAL_ERROR "PAPILO not found, try specifying PAPILO_DIR. ${AUTOBUILD_MSG}")
    endif()
endif()
if(PAPILO_FOUND)
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
    find_package(SOPLEX REQUIRED CONFIG HINTS ${SOPLEX_DIR})
    if(NOT SOPLEX_FOUND)
        message(FATAL_ERROR "Requested LP solver SoPlex not found.")
    endif()
    if (DEFINED SOPLEX_WITH_PAPILO)
        message(STATUS "SOPLEX links PAPILO")
        if((NOT SCIP_WITH_PAPILO)) # TODO not sure how to handle AUTOBUILD
            message(FATAL_ERROR "Please link a SOPLEX built without PAPILO to build SCIP without PAPILO.")
        endif ()
    endif ()
    if(LPSCHECK)
        find_package(CPLEX REQUIRED)
        if(NOT CPLEX_FOUND)
            message(FATAL_ERROR "Requested LP check solver CPLEX not found.")
        endif()
        set(SCIP_WITH_LPSCHECK on)
    endif()
elseif(LPS STREQUAL "cpx")
    find_package(CPLEX REQUIRED)
    if(NOT CPLEX_FOUND)
        message(FATAL_ERROR "Requested LP solver CPLEX not found.")
    endif()
elseif(LPS STREQUAL "glop")
    if(NOT TPI STREQUAL "tny")
        message(FATAL_ERROR "When using GLOP as an lp solver, TPI needs to be set to 'tny'.")
    endif()
    find_package(GLOP REQUIRED)
    if(NOT GLOP_FOUND)
        message(FATAL_ERROR "Requested LP solver Glop not found.")
    endif()
elseif(LPS STREQUAL "grb")
    find_package(GUROBI REQUIRED)
    if(NOT GUROBI_FOUND)
        message(FATAL_ERROR "Requested LP solver Gurobi not found.")
    endif()
elseif(LPS STREQUAL "qso")
    find_package(QSO REQUIRED)
    if(NOT QSO_FOUND)
        message(FATAL_ERROR "Requested LP solver Qsopt not found.")
    endif()
    if(BUILD_SHARED_LIBS)
        set(BUILD_SHARED_LIBS off)
        message(STATUS "LP solver QSopt turns off shared libraries.")
    endif()
elseif(LPS STREQUAL "clp")
    find_package(CLP REQUIRED)
    if(NOT CLP_FOUND)
        message(FATAL_ERROR "Requested LP solver CLP not found.")
    endif()
elseif(LPS STREQUAL "xprs")
    find_package(XPRESS REQUIRED)
    if(NOT XPRESS_FOUND)
        message(FATAL_ERROR "Requested LP solver XPRESS not found.")
    endif()
elseif(LPS STREQUAL "msk")
    find_package(MOSEK REQUIRED)
    if(NOT MOSEK_FOUND)
        message(FATAL_ERROR "Requested LP solver MOSEK not found.")
    endif()
elseif(LPS STREQUAL "highs")
    find_package(HIGHS REQUIRED CONFIG HINTS ${HIGHS_DIR})
    if(NOT HIGHS_FOUND)
        message(FATAL_ERROR "Requested LP solver HiGHS not found.")
    endif()
    find_package(Threads REQUIRED)
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
    set(LPS_PIC_LIBRARIES ${SOPLEX_PIC_LIBRARIES})
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
    set(LPS_PIC_LIBRARIES ${CLP_LIBRARIES})
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
    set(LPS_PIC_LIBRARIES ${LPS_PIC_LIBRARIES} ${CPLEX_LIBRARIES})
else()
    message(STATUS "Support CPLEX: OFF")
endif()

if(GLOP_FOUND)
    message(STATUS "Finding GLOP - found")
    include_directories(${GLOP_INCLUDE_DIRS})
    set(lpi lpi/lpi_glop.cpp)
    set(LPS_LIBRARIES ${LPS_LIBRARIES} ${GLOP_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_PIC_LIBRARIES} ${GLOP_LIBRARIES})
    set(CMAKE_CXX_STANDARD 17)
else()
    message(STATUS "Support GLOP: OFF")
endif()

if(GUROBI_FOUND)
    message(STATUS "Finding GUROBI - found")
    include_directories(${GUROBI_INCLUDE_DIRS})
    set(lpi lpi/lpi_grb.c)
    set(LPS_LIBRARIES ${LPS_LIBRARIES} ${GUROBI_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_PIC_LIBRARIES} ${GUROBI_LIBRARIES})
else()
    message(STATUS "Support GUROBI: OFF")
endif()

if(XPRESS_FOUND)
    message(STATUS "Finding XPRESS - found")
    include_directories(${XPRESS_INCLUDE_DIRS})
    set(lpi lpi/lpi_xprs.c)
    set(LPS_LIBRARIES ${LPS_LIBRARIES} ${XPRESS_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_PIC_LIBRARIES} ${XPRESS_LIBRARIES})
else()
    message(STATUS "Support XPRESS: OFF")
endif()

if(MOSEK_FOUND)
    message(STATUS "Finding MOSEK - found")
    include_directories(${MOSEK_INCLUDE_DIRS})
    set(lpi lpi/lpi_msk.c)
    set(LPS_LIBRARIES ${LPS_LIBRARIES} ${MOSEK_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_PIC_LIBRARIES} ${MOSEK_LIBRARIES})
else()
    message(STATUS "Support MOSEK: OFF")
endif()

if(HIGHS_FOUND)
    if(LPS STREQUAL "highs")
        set(lpi lpi/lpi_highs.cpp)
    endif()
    set(LPS_LIBRARIES ${LPS_LIBRARIES} highs::highs)
    set(LPS_PIC_LIBRARIES ${LPS_PIC_LIBRARIES} highs::highs)
    message(WARNING "EXPERIMENTAL: You are using the experimental HiGHS interface. Please report any issues on https://github.com/scipopt/scip.")
endif()

if(QSO_FOUND)
    message(STATUS "Finding QSO - found")
    include_directories(${QSO_INCLUDE_DIRS})
    set(lpi lpi/lpi_qso.c)
    set(LPS_LIBRARIES ${LPS_LIBRARIES} ${QSO_LIBRARIES})
    set(LPS_PIC_LIBRARIES ${LPS_PIC_LIBRARIES} ${QSO_LIBRARIES})
else()
    message(STATUS "Support QSO: OFF")
endif()

# Boost
if(SCIP_WITH_EXACTSOLVE)
    set(BOOST_REQUIRED_VERSION 1.68.0)
    if(Boost_FOUND AND (Boost_VERSION VERSION_GREATER_EQUAL ${BOOST_REQUIRED_VERSION}))
        message(STATUS "Already found Boost version \"${Boost_VERSION_STRING}\", minimum required is \"${BOOST_REQUIRED_VERSION}\".")
    else()
        find_package(Boost ${BOOST_REQUIRED_VERSION})
    endif()
    if(Boost_FOUND AND (Boost_VERSION VERSION_GREATER_EQUAL ${BOOST_REQUIRED_VERSION}))
        set(SCIP_WITH_BOOST TRUE)
    else()
        set(SCIP_WITH_EXACTSOLVE OFF)
        if(EXACTSOLVE AND NOT (EXACTSOLVE_UPPER STREQUAL "AUTO")) # EXACTSOLVE was explicitly requested.
            message(FATAL_ERROR "Boost version ${BOOST_REQUIRED_VERSION} or higher not found, which is required for EXACTSOLVE; try specifying Boost_DIR.")
        else() # EXACTSOLVE is automatic.
        endif()
    endif()
endif()

# Find exact LP solvers.
set(TEST_EXAXTSOLVE OFF)
set(SCIP_WITH_LPSEXACT OFF) # Human-readable name of exact LP solver
if(SCIP_WITH_EXACTSOLVE)
    # SoPlex
    if((LPSEXACT_UPPER STREQUAL "SPX") OR ((LPSEXACT_UPPER STREQUAL "AUTO" AND NOT SCIP_WITH_LPSEXACT)))
        set(BAD_SOPLEX_REASON 0)
        if(SOPLEX_FOUND)
            message(STATUS "Already found LP solver SoPlex.")
            if(NOT SOPLEX_WITH_MPFR)
                set(BAD_SOPLEX_REASON "does not support MPFR")
            endif()
        else()
            message(STATUS "Finding exact LP solver SoPlex.")
            find_package(SOPLEX CONFIG HINTS ${SOPLEX_DIR})
        endif()
        if(NOT SOPLEX_FOUND)
            set(BAD_SOPLEX_REASON "not found")
        elseif(NOT SOPLEX_WITH_MPFR)
            set(BAD_SOPLEX_REASON "does not support MPFR")
        endif()
        if(SOPLEX_FOUND AND SOPLEX_WITH_MPFR)
            set(SCIP_WITH_LPSEXACT "SoPlex")
            set(TEST_EXAXTSOLVE ON)
        elseif(LPSEXACT_UPPER STREQUAL "SPX" AND NOT (EXACTSOLVE_UPPER STREQUAL "AUTO")) # SoPlex and EXACTSOLVE were explicitly requested.
            message(FATAL_ERROR "LP solver SoPlex ${BAD_SOPLEX_REASON}, which is required for EXACTSOLVE; try specifying SOPLEX_DIR.")
        elseif(LPSEXACT_UPPER STREQUAL "SPX") # SoPlex was explicitly requested, but EXACTSOLVE is AUTO.
            message(WARNING "LP solver SoPlex ${BAD_SOPLEX_REASON}, which is required for EXACTSOLVE; try specifying SOPLEX_DIR.")
            set(SCIP_WITH_EXACTSOLVE OFF)
        else()
            message(STATUS "LP solver SoPlex ${BAD_SOPLEX_REASON}.")
        endif()
    endif()
    # QSopt_ex
    if((LPSEXACT_UPPER STREQUAL "QSOEX") OR ((LPSEXACT_UPPER STREQUAL "AUTO" AND NOT SCIP_WITH_LPSEXACT)))
        message(STATUS "Finding exact LP solver QSopt_ex.")
        find_package(QSOEX)
        if(QSOEX_FOUND)
            message(STATUS "Finding QSopt_ex dependency EGLIB.")
            find_package(EGLIB)
        endif()
        if(QSOEX_FOUND AND EGLIB_FOUND)
            set(SCIP_WITH_LPSEXACT "QSopt_ex")
            set(SCIP_WITH_EGLIB ON)
            set(TEST_EXAXTSOLVE ON)
        elseif(LPSEXACT_UPPER STREQUAL "QSOEX" AND NOT (EXACTSOLVE_UPPER STREQUAL "AUTO")) # QSopt_ex and EXACTSOLVE were explicitly requested.
            message(FATAL_ERROR "Exact LP solver QSopt_ex not found, which is required for EXACTSOLVE; try specifying QSOEX_DIR.")
        elseif(LPSEXACT_UPPER STREQUAL "QSOEX") # QSopt_ex was explicitly requested, but EXACTSOLVE is AUTO.
            message(WARNING "Exact LP solver QSopt_ex not found, which is required for EXACTSOLVE; try specifying QSOEX_DIR.")
            set(SCIP_WITH_EXACTSOLVE OFF)
        endif()
    endif()
    # No exact LP solver found.
    if(NOT SCIP_WITH_LPSEXACT AND (EXACTSOLVE_UPPER STREQUAL "AUTO"))
        message(STATUS "No exact LP solver found, which is required for EXACTSOLVE.")
        set(SCIP_WITH_EXACTSOLVE OFF)
        set(SCIP_WITH_BOOST FALSE)
    elseif(NOT SCIP_WITH_LPSEXACT AND EXACTSOLVE)
        message(FATAL_ERROR "No exact LP solver found, which is required for EXACTSOLVE.")
    else()
        set(SCIP_WITH_BOOST TRUE)
    endif()
else()
    set(SCIP_WITH_LPSEXACT "none")
    set(SCIP_WITH_BOOST FALSE)
endif()

# If EXACTSOLVE was disabled, we also do not need MPFR.
if(NOT SCIP_WITH_EXACTSOLVE)
    set(SCIP_WITH_MPFR OFF)
endif()

# Now we are sure whether we link to MPFR.
if(SCIP_WITH_MPFR)
    # TODO: Avoid global usage.
    include_directories(${MPFR_INCLUDE_DIRS})
else()
    set(MPFR_LIBRARIES "")
endif()

# Now we are sure about EXACTSOLVE.
if(SCIP_WITH_EXACTSOLVE)
    message(STATUS "Building SCIP with support for exact solving mode using exact LP solver ${SCIP_WITH_LPSEXACT}.")

    message(STATUS "Finding VIPR checker")
    find_program(VIPRCOMP "viprcomp" HINTS ${VIPR_DIR} $ENV{VIPR_DIR})
    find_program(VIPRCHK "viprchk" HINTS ${VIPR_DIR} $ENV{VIPR_DIR})
    if(VIPRCOMP AND VIPRCHK)
        message(STATUS "Finding VIPR checker - found")
    else()
        message(STATUS "Finding VIPR checker - not found: verification is not enabled for exact ctest. Make sure to install VIPR (https://github.com/scipopt/vipr) to enable certificate checking.")
    endif()

    # TODO: Avoid global usage.
    include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

    if(SCIP_WITH_LPSEXACT STREQUAL "QSopt_ex")
        if(BUILD_SHARED_LIBS)
            set(BUILD_SHARED_LIBS off)
            message(STATUS "Exact LP solver QSopt_ex turns off shared libraries.")
        endif()
        include_directories(${QSOEX_INCLUDE_DIRS})
        include_directories(${EGLIB_INCLUDE_DIRS})
        set(LPS_LIBRARIES ${LPS_LIBRARIES} ${QSOEX_LIBRARIES} ${EGLIB_LIBRARIES})
        set(LPS_PIC_LIBRARIES ${LPS_PIC_LIBRARIES} ${QSOEX_LIBRARIES} ${EGLIB_LIBRARIES})
        set(SCIP_WITH_QSOEX ON)
        set(lpiexact lpiexact/lpiexact_qsoex.c)
    endif()

    if(SCIP_WITH_LPSEXACT STREQUAL "SoPlex")
        include_directories(${SOPLEX_INCLUDE_DIRS})
        set(lpiexact lpiexact/lpiexact_spx.cpp)
        set(LPS_LIBRARIES ${LPS_LIBRARIES} ${SOPLEX_LIBRARIES})
        set(LPS_PIC_LIBRARIES ${LPS_PIC_LIBRARIES} ${SOPLEX_PIC_LIBRARIES})
    endif()

    if(SCIP_WITH_LPSEXACT STREQUAL "none")
        set(lpiexact lpiexact/lpiexact_none.c)
    endif()
else()
    message(STATUS "Building SCIP without support for exact solving mode.")
    set(lpiexact lpiexact/lpiexact_none.c)
endif()

# Now we now whether we actually build shared libraries.
set(SHARED ${BUILD_SHARED_LIBS})
if(SHARED)
  message(STATUS "Building shared libraries.")
else()
  message(STATUS "Building static libraries.")
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
if(ZIMPL OR AUTOBUILD)
    message(STATUS "Finding ZIMPL")
    find_package(ZIMPL 3.5.0 CONFIG HINTS ${ZIMPL_DIR})
    if((NOT ZIMPL_FOUND) AND (NOT AUTOBUILD))
        message(FATAL_ERROR "ZIMPL not found, try specifying ZIMPL_DIR. ${AUTOBUILD_MSG}")
    endif()
endif()
if(ZIMPL_FOUND)
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

if(AMPL)
    set(SCIP_WITH_AMPL ON)
    include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/src/amplmp/include)
    set(amplsrc scip/reader_nl.cpp amplmp/src/format.cpp amplmp/src/expr-info.cpp amplmp/src/nl-reader.cpp amplmp/src/os.cpp amplmp/src/posix.cpp)
else()
    message(STATUS "Support AMPL: OFF")
    set(SCIP_WITH_AMPL OFF)
endif()

# IPOPT
if(IPOPT OR AUTOBUILD)
    message(STATUS "Finding IPOPT")
    find_package(IPOPT 3.12.0)
    if((NOT IPOPT_FOUND) AND (NOT AUTOBUILD))
        message(FATAL_ERROR "IPOPT not found, try specifying IPOPT_DIR. ${AUTOBUILD_MSG}")
    endif()
endif()
if(IPOPT_FOUND)
    message(STATUS "Finding IPOPT - found")
    include_directories(${IPOPT_INCLUDE_DIRS})
    #on debian IPOPT package needs this definition to work
    set(HAVE_CSTDDEF on)
    set(NLPI_LIBRARIES ${IPOPT_LIBRARIES})
    set(nlpi scip/nlpi_ipopt.cpp)
    set(LAPACK off)
else()
    message(STATUS "Support IPOPT: OFF")
    set(nlpi scip/nlpi_ipopt_dummy.c)
endif()

# LAPACK
if(LAPACK)
    message(STATUS "Finding LAPACK")
    find_package(LAPACK REQUIRED)
    message(STATUS "Finding LAPACK - found")
    set(SCIP_WITH_LAPACK on)
else()
    message(STATUS "Support LAPACK: OFF")
    set(LAPACK_LIBRARIES "")
    set(SCIP_WITH_LAPACK off)
endif()

# WORHP
if(WORHP OR AUTOBUILD)
    message(STATUS "Finding WORHP")
    find_package(WORHP)
    if((NOT WORHP_FOUND) AND (NOT AUTOBUILD))
        message(FATAL_ERROR "WORHP not found, try specifying WORPH_DIR. ${AUTOBUILD_MSG}")
    endif()
endif()
if(WORHP_FOUND)
    message(STATUS "Finding WORHP - found")
    include_directories(${WORHP_INCLUDE_DIRS})
    set(nlpi ${nlpi} scip/nlpi_worhp.c)
    set(NLPI_LIBRARIES ${NLPI_LIBRARIES} ${WORHP_LIBRARIES})
else()
    message(STATUS "Support WORHP: OFF")
    set(nlpi ${nlpi} scip/nlpi_worhp_dummy.c)
endif()

# FilterSQP (with CMake, nlpi_filtersqp doesn't build anyway)
set(nlpi ${nlpi} scip/nlpi_filtersqp_dummy.c)

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

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OLD_CMAKE_RUNTIME_OUTPUT_DIRECTORY})

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    include(FeatureSummary)
    feature_summary(WHAT ALL)
endif()

# --------------------------------------------------------------
# write log file of non-default parameter settings
file(WRITE ${CMAKE_BINARY_DIR}/cmake.log "")

macro(log_option option)
    string(TOUPPER ${${option}} UPPER)
    if(NOT (${UPPER} STREQUAL "AUTO"))
        file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -D${option}=${${option}}")
    endif()
    if(${option}_DIR)
        file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -D${option}_DIR=${${option}_DIR}")
    endif()
endmacro()

log_option(EXACTSOLVE)
log_option(LPSEXACT)
log_option(MPFR)

if(AUTOBUILD)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log "-DAUTOBUILD=on")
endif()
if(NOT SHARED)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DSHARED=off")
endif()
if(NOT ZLIB)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DZLIB=off")
endif()
if(NOT READLINE)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DREADLINE=off")
endif()
if(NOT GMP)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DGMP=off")
endif()
if(STATIC_GMP)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DSTATIC_GMP=on")
endif()
if(NOT PAPILO)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DPAPILO=off")
else()
  if(PAPILO_DIR)
    file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DPAPILO_DIR=${PAPILO_DIR}")
  endif()
endif()
if(NOT ZIMPL)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DZIMPL=off")
else()
  if(ZIMPL_DIR)
    file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DZIMPL_DIR=${ZIMPL_DIR}")
  endif()
endif()
if(NOT AMPL)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DAMPL=off")
endif()
if(NOT IPOPT)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DIPOPT=off")
else()
  if(IPOPT_DIR)
    file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DIPOPT_DIR=${IPOPT_DIR}")
  endif()
endif()
if(LAPACK)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DLAPACK=on")
endif()
if(WORHP)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DWORHP=on")
  if(WORHP_DIR)
    file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DWORHP_DIR=${WORHP_DIR}")
  endif()
endif()
if(PARASCIP)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DPARASCIP=on")
endif()
if(NOT THREADSAFE)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DTHREADSAFE=off")
endif()
if(LPSCHECK)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DLPSCHECK=on")
endif()
if(LEGACY)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DLEGACY=on")
endif()
if(NOBLKBUFMEM)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DNOBLKBUFMEM=on")
endif()
if(NOBLKMEM)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DNOBLKMEM=on")
endif()
if(NOBUFMEM)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DNOBUFMEM=on")
endif()
if(DEBUGSOL)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DDEBUGSOL=on")
endif()
if(SANITIZE)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DSANITIZE=${SANITIZE}")
endif()
if(COVERAGE)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DCOVERAGE=on")
endif()
if(MT)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DMT=on")
endif()
if(CXXONLY)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DCXXONLY=on")
endif()
if(TPI STREQUAL "tny")
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DTPI=tny")
elseif(TPI STREQUAL "omp")
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DTPI=omp")
endif()
if(EXPRINT STREQUAL "none")
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DEXPRINT=none")
endif()
if(LPS STREQUAL "spx")
  if(SOPLEX_DIR)
    file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DSOPLEX_DIR=${SOPLEX_DIR}")
  endif()
elseif(LPS STREQUAL "cpx")
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DLPS=${LPS}")
  if(CPLEX_DIR)
    file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DCPLEX_DIR=${CPLEX_DIR}")
  endif()
elseif(LPS STREQUAL "glop")
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DLPS=${LPS}")
  if(GLOP_DIR)
    file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DGLOP_DIR=${GLOP_DIR}")
  endif()
elseif(LPS STREQUAL "grb")
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DLPS=${LPS}")
  if(GUROBI_DIR)
    file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DGUROBI_DIR=${GUROBI_DIR}")
  endif()
elseif(LPS STREQUAL "qso")
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DLPS=${LPS}")
  if(QSO_DIR)
    file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DQSO_DIR=${QSOI_DIR}")
  endif()
elseif(LPS STREQUAL "clp")
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DLPS=${LPS}")
  if(CLP_DIR)
    file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DCLP_DIR=${CLP_DIR}")
  endif()
elseif(LPS STREQUAL "xprs")
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DLPS=${LPS}")
  if(XPRESS_DIR)
    file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DXPRESS_DIR=${XPRESS_DIR}")
  endif()
elseif(LPS STREQUAL "msk")
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DLPS=${LPS}")
  if(MOSEK_DIR)
    file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DMOSEK_DIR=${MOSEK_DIR}")
  endif()
elseif(LPS STREQUAL "highs")
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DLPS=${LPS}")
elseif(LPS STREQUAL "none")
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DLPS=${LPS}")
endif()
if(NOT SYM STREQUAL "sbliss")
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DSYM=${SYM}")
endif()
if((SYM STREQUAL "nauty" OR SYM STREQUAL "snauty") AND NAUTY_DIR)
  file(APPEND ${CMAKE_BINARY_DIR}/cmake.log " -DNAUTY_DIR=${NAUTY_DIR}")
endif()
file(APPEND ${CMAKE_BINARY_DIR}/cmake.log "\n")
message(STATUS "Written file with changed options to ${CMAKE_BINARY_DIR}/cmake.log.")
# --------------------------------------------------------------
