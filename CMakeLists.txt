cmake_minimum_required(VERSION 3.9)

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

set(SCIP_VERSION_MAJOR 8)
set(SCIP_VERSION_MINOR 1)
set(SCIP_VERSION_PATCH 0)
set(SCIP_VERSION_SUB 0)
set(SCIP_VERSION_API 104)

project(SCIP
    VERSION ${SCIP_VERSION_MAJOR}.${SCIP_VERSION_MINOR}.${SCIP_VERSION_PATCH}.${SCIP_VERSION_SUB}
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
option(SANITIZE_UNDEFINED "should the undefined behavior sanitizer be enabled in debug mode if available" OFF)
option(SANITIZE_THREAD "should the thread sanitizer be enabled in debug mode if available" OFF)
option(COVERAGE "enable coverage support" OFF)
SET(COVERAGE_CTEST_ARGS "" CACHE STRING "additional ctest arguments for coverage")
option(MT "use static runtime libraries for Visual Studio compiler" OFF)
option(CXXONLY "use a c++ compiler for all source files" OFF)

set(TPI none CACHE STRING "options for thread support library")  #create the variable
set_property(CACHE TPI PROPERTY STRINGS none tny omp)  #define list of values GUI will offer for the variable

set(EXPRINT cppad CACHE STRING "options for expression interpreter")  #create the variable
set_property(CACHE EXPRINT PROPERTY STRINGS none cppad )  #define list of values GUI will offer for the variable

set(LPS spx CACHE STRING "options for LP solver")  #create the variable
set_property(CACHE LPS PROPERTY STRINGS spx cpx grb xprs clp glop msk qso none )  #define list of values GUI will offer for the variable

set(SYM sbliss CACHE STRING "options for symmetry computation")  #create the variable
set_property(CACHE SYM PROPERTY STRINGS bliss sassy nauty none )  #define list of values GUI will offer for the variable

if(NOT (CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" OR CMAKE_BUILD_TYPE STREQUAL "MinSizeRel"))
    message(FATAL_ERROR "You have not selected a valid CMAKE_BUILD_TYPE. Please select 'Release' or 'Debug'.")
endif()

if(AUTOBUILD)
    message(WARNING "You are using the automatic build option AUTOBUILD=ON to find and use dependencies automatically.")
else()
    message(WARNING "You are not using the automatic build option. ${AUTOBUILD_MSG}")
endif() # AUTOBUILD

# use C++14 standard
set(CMAKE_CXX_STANDARD 14)

#search the selected symmetry computation program
message(STATUS "Finding symmetry computation program \"${SYM}\"")

unset(BLISS_TARGET)
unset(NAUTY_TARGET)
if(SYM STREQUAL "bliss")
    message(STATUS "Support SYM: bliss")
    set(sym symmetry/compute_symmetry_bliss.cpp)

    # modify configuration for bliss
    set(BUILD_SHARED_LIBS OFF)
    set(TMPFLAGS ${CMAKE_C_FLAGS})
    set(TMXFLAGS ${CMAKE_CXX_FLAGS})
    set(CMAKE_C_FLAGS -w)
    set(CMAKE_CXX_FLAGS -w)

    add_subdirectory(src/bliss EXCLUDE_FROM_ALL)

    # undo modification
    set(CMAKE_C_FLAGS ${TMPFLAGS})
    set(CMAKE_CXX_FLAGS ${TMXFLAGS})
    set(BUILD_SHARED_LIBS ${SHARED})

    set(BLISS_TARGET libbliss)
    set(SYM_LIBRARIES libbliss)
    set(SYM_PIC_LIBRARIES libbliss)
elseif(SYM STREQUAL "sbliss")
    message(STATUS "Support SYM: sbliss")
    set(sym symmetry/compute_symmetry_sassy.cpp)

    # sassy is currently based on bliss, so configure bliss
    set(BUILD_SHARED_LIBS OFF)
    set(TMPFLAGS ${CMAKE_C_FLAGS})
    set(TMXFLAGS ${CMAKE_CXX_FLAGS})
    set(CMAKE_C_FLAGS -w)
    set(CMAKE_CXX_FLAGS -w)

    add_subdirectory(src/bliss)

    # undo modification
    set(CMAKE_C_FLAGS ${TMPFLAGS})
    set(CMAKE_CXX_FLAGS ${TMXFLAGS})
    set(BUILD_SHARED_LIBS ${SHARED})

    set(BLISS_TARGET libbliss)
    set(SYM_LIBRARIES libbliss)
    set(SYM_PIC_LIBRARIES libbliss)

    # sassy needs C++-17
    set(CMAKE_CXX_STANDARD 17)

    if(NOT BLISS_INCLUDE_DIR OR NOT BLISS_LIBRARY)
      find_path(BLISS_INCLUDE_DIR
	NAMES bliss/graph.hh
	PATH_SUFFIXES include src)
    endif()
elseif(SYM STREQUAL "nauty")
    message(STATUS "Support SYM: nauty")
    set(sym symmetry/compute_symmetry_nauty.c)

    if(SHARED)
        set(PIC_FLAG "CFLAGS='-fPIC'")
    else()
        set(PIC_FLAG "")
    endif()
    file(GLOB files ${CMAKE_CURRENT_SOURCE_DIR}/src/nauty/*)
    # add custom command to build the library using configure and make
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lib/libnauty.a
        COMMAND 
        mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/src/nauty &&
        cd ${CMAKE_CURRENT_BINARY_DIR}/src/nauty &&
        cp ${CMAKE_CURRENT_SOURCE_DIR}/src/nauty/* . &&
        ./configure &&
        make ${PIC_FLAG} nauty.a &&
        cp nauty.a ${CMAKE_CURRENT_BINARY_DIR}/lib/libnauty.a
        DEPENDS ${files})
    
    # add custom target and library target
    add_custom_target(libnauty_target DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/lib/libnauty.a)
    add_library(libnauty STATIC IMPORTED GLOBAL)
    add_dependencies(libnauty libnauty_target)
    set_target_properties(libnauty PROPERTIES 
        IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/lib/libnauty.a)

    set(NAUTY_TARGET libnauty)
    set(SYM_LIBRARIES libnauty)
    set(SYM_PIC_LIBRARIES libnauty)

elseif(SYM STREQUAL "none")
    message(STATUS "Support SYM: OFF")
    set(sym symmetry/compute_symmetry_none.cpp)
else()
    message(FATAL_ERROR "option SYM has wrong value")
endif()

# if changing these flags, also update GCCWARN/GXXWARN in make/make.project
set(ADD_C_AND_CXX_FLAGS -pedantic -Wall -W -Wpointer-arith -Wcast-align -Wwrite-strings -Wshadow  -Wredundant-decls -Wdisabled-optimization -Wno-long-long -Wno-unknown-pragmas -Wno-unused-parameter)
set(ADD_C_FLAGS -Wsign-compare -Wstrict-prototypes -Wmissing-declarations -Wmissing-prototypes -Wdeclaration-after-statement -Wno-override-init -Wno-uninitialized -Wno-maybe-uninitialized)
set(ADD_CXX_FLAGS -Wctor-dtor-privacy -Wnon-virtual-dtor -Wreorder -Woverloaded-virtual -Wsign-promo -Wsynth -Wcast-qual -Wno-strict-overflow -Wno-psabi)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
# clang specific suppressions
set(ADD_C_AND_CXX_FLAGS ${ADD_C_AND_CXX_FLAGS} -Wno-unknown-warning-option)
set(ADD_C_FLAGS ${ADD_C_FLAGS} -Wno-initializer-overrides)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
   # for intel compilers we need to set the fp-model precise flag to disable incorrect floating point optimizations
   add_compile_options(-fp-model precise)
   # to prevent symbols from svml exported twice we use the -shared-intel flag, see https://git.zib.de/integer/scip/issues/2872
   # this happens with the cplex static lib
   if(LPS MATCHES "cpx")
       add_compile_options(-shared-intel) # for c and cxx
   endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # require at least gcc 4.8
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
        message(WARNING "GCC version not supported, should be at least 4.8!")
    endif()
    add_compile_options(${ADD_C_AND_CXX_FLAGS})
    add_compile_options("$<$<COMPILE_LANGUAGE:C>:${ADD_C_FLAGS}>$<$<COMPILE_LANGUAGE:CXX>:${ADD_CXX_FLAGS}>")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(${ADD_C_AND_CXX_FLAGS})
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
    message(STATUS "Finding Readline")
    find_package(Readline)
    if((NOT Readline_FOUND) AND (NOT AUTOBUILD))
        message(FATAL_ERROR "READLINE not found, try specifying READLINE_DIR. ${AUTOBUILD_MSG}")
    endif()
endif()
if(Readline_FOUND)
    message(STATUS "Finding Readline - found")
    include_directories(${Readline_INCLUDE_DIRS})
    set(SCIP_WITH_READLINE ON)
else()
    message(STATUS "Support Readline: OFF")
    set(Readline_LIBRARY "")
    set(SCIP_WITH_READLINE OFF)
endif()

# GMP
if(GMP OR AUTOBUILD)
    message(STATUS "Finding GMP")
    find_package(GMP)
    if((NOT GMP_FOUND) AND (NOT AUTOBUILD))
        message(FATAL_ERROR "GMP not found, try specifying GMP_DIR. ${AUTOBUILD_MSG}")
    endif()
endif()
if(GMP_FOUND)
    message(STATUS "Finding GMP - found")
    include_directories(${GMP_INCLUDE_DIRS})
    set(SCIP_WITH_GMP ON)
else()
    message(STATUS "Support GMP: OFF")
    set(GMP_LIBRARIES "")
    set(SCIP_WITH_GMP OFF)
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
    if (DEFINED SOPLEX_WITH_PAPILO)
        message(STATUS "SOPLEX links PAPILO")
        if((NOT SCIP_WITH_PAPILO)) # TODO not sure how to handle AUTOBUILD
            message(FATAL_ERROR "Please link a SOPLEX built without PAPILO to build SCIP without PAPILO.")
        endif ()
    endif ()
    if(LPSCHECK)
        find_package(CPLEX REQUIRED)
        set(SCIP_WITH_LPSCHECK on)
    endif()
elseif(LPS STREQUAL "cpx")
    find_package(CPLEX REQUIRED)
elseif(LPS STREQUAL "glop")
    if(NOT TPI STREQUAL "tny")
        message(FATAL_ERROR "When using GLOP as an lp solver, TPI needs to be set to 'tny'.")
    endif()
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
else()
    message(STATUS "Support IPOPT: OFF")
    set(nlpi scip/nlpi_ipopt_dummy.c)
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
