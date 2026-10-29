# Overall structure
#
# POLICIES
# PROJECT DEFINITIONS
#    Version
#    Project declaration
#    Minimum dependency versions
#    Directories and search paths
#    Packaging
# OPTIONS
#    Declarations with preprocessing/printing
# BASIC SETTINGS
#    Shared/static libraries
#    Build type
#    LTO
#    Githash
# DEPENDENCIES
#    Search and process
#    Check for position-independent code
# TOOLS
#    Code coverage.
#    Global compilation settings
# SUBDIRECTORIES
# INSTALLATION
#    Licenses
#    Config(-version) files
#    Build-tree targets
# LOGGING
#    Non-default parameters

cmake_minimum_required(VERSION 3.11)

##### POLICIES #####

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

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

# Use rpaths on Mac OS, i.e., CMAKE_MACOSX_RPATH is ON by default.
if(POLICY CMP0042)
   cmake_policy(SET CMP0042 NEW)
endif()

##### PROJECT DEFINITIONS #####

# SCIP version - it is cached in order to be readable for super project (SCIPOptSuite).
set(SCIP_VERSION_MAJOR 11
   CACHE STRING "Major component of SCIP version")
set(SCIP_VERSION_MINOR 0
   CACHE STRING "Minor component of SCIP version")
set(SCIP_VERSION_PATCH 0
   CACHE STRING "Patch component of SCIP version")
set(SCIP_VERSION_API 167
   CACHE STRING "API component of SCIP version")

# Project
project(SCIP
   VERSION ${SCIP_VERSION_MAJOR}.${SCIP_VERSION_MINOR}.${SCIP_VERSION_PATCH}
   LANGUAGES CXX C
   DESCRIPTION "SCIP is currently one of the fastest non-commercial solvers for mixed integer programming (MIP) and mixed integer nonlinear programming (MINLP)"
   )

set(BOOST_MIN_VERSION 1.68) # Minimum Boost version if we use it at all.
set(IPOPT_MIN_VERSION 3.12) # Minimum Ipopt version if we use it at all.
set(ZIMPL_MIN_VERSION 3.5.0) # Minimum ZIMPL version if we use it at all.

set(NAUTY_SHIPPED_VERSION 2.8.8) # Version of nauty that is shipped with SCIP.
set(DEJAVU_SHIPPED_VERSION 2.1) # Version of dejavu that is shipped with SCIP.
set(CPPAD_SHIPPED_VERSION 20180000.0) # Version of CppAD that is shipped with SCIP.

# Set installation directories according to GNU coding standards.
include(GNUInstallDirs)

# Path to CMake modules such as FindGMP.cmake.
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Packaging
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${SCIP_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${SCIP_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${SCIP_VERSION_PATCH}")
set(CPACK_PACKAGE_VENDOR "Zuse Institute Berlin")
set(CPACK_PACKAGE_CONTACT "http://scipopt.org")
include(CPack)


##### OPTIONS #####

# Whether to build shared libraries.
option(BUILD_SHARED_LIBS "Whether to build shared libraries" ON)
# For compatibility with the Makefile build system, the SHARED option exists.
set(SHARED AUTO CACHE STRING "Build shared libraries (AUTO: use BUILD_SHARED_LIBS)")
set_property(CACHE SHARED PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${SHARED}" SHARED)

option(AMPL "should AMPL interface be compiled" ON)

set(CHECKSTAGE AUTO CACHE STRING "enable stage checks; auto enables it in debug builds")
set_property(CACHE CHECKSTAGE PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${CHECKSTAGE}" CHECKSTAGE)

set(CHECKNAME AUTO CACHE STRING "enable name checks; auto enables it in debug builds")
set_property(CACHE CHECKNAME PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${CHECKNAME}" CHECKNAME)

set(CONOPT AUTO CACHE STRING "should SCIP link to CONOPT")
set_property(CACHE CONOPT PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${CONOPT}" CONOPT)

option(COVERAGE "enable coverage support" OFF)
set(COVERAGE_CTEST_ARGS "" CACHE STRING "additional ctest arguments for coverage")

option(DEBUGSOL "should the debug solution mechanism be enabled" OFF)

set(EXACTSOLVE AUTO CACHE STRING "compile with support for exact rational solve (requires GMP, MPFR and Boost)")
set_property(CACHE EXACTSOLVE PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${EXACTSOLVE}" EXACTSOLVE)

set(EXPRINT cppad CACHE STRING "options for expression interpreter")
set_property(CACHE EXPRINT PROPERTY STRINGS none cppad)
string(TOLOWER "${EXPRINT}" EXPRINT)

option(GITHASH_BRANCH "append the current branch to the githash" OFF)

set(GMP AUTO CACHE STRING "should SCIP link to GMP")
set_property(CACHE GMP PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${GMP}" GMP)

set(IPOPT AUTO CACHE STRING "should SCIP link to Ipopt")
set_property(CACHE IPOPT PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${IPOPT}" IPOPT)

set(LAPACK AUTO CACHE STRING "should SCIP link to LAPACK")
set_property(CACHE LAPACK PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${LAPACK}" LAPACK)

set(LPS AUTO CACHE STRING "options for floating-point LP solver")
set_property(CACHE LPS PROPERTY STRINGS spx clp qso glop highs cpx xprs msk grb none AUTO)
string(TOLOWER "${LPS}" LPS)

set(LPSEXACT AUTO CACHE STRING "options for exact LP solver")
set_property(CACHE LPSEXACT PROPERTY STRINGS AUTO spx qsoex none)
string(TOLOWER "${LPSEXACT}" LPSEXACT)
if(NOT (LPSEXACT MATCHES "auto|spx|qsoex|none"))
   message(FATAL_ERROR "Option LPSEXACT has invalid value \"${LPSEXACT}\"! Possible are `auto`, `spx`, `qsoex` and `none`.")
endif()

option(LPSCHECK "double check SoPlex results with CPLEX" OFF)

set(LTO OFF CACHE STRING "enable Link Time Optimization = Interprocedural Optimization")
set_property(CACHE LTO PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${LTO}" LTO)

set(MPFR AUTO CACHE STRING "should SCIP link to MPFR")
set_property(CACHE MPFR PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${MPFR}" MPFR)

option(NOBLKBUFMEM "should block and buffer memory be disabled" OFF)

option(NOBLKMEM "should block memory be disabled" OFF)

option(NOBUFMEM "should buffer memory be disabled" OFF)

set(PAPILO AUTO CACHE STRING "should SCIP link to PaPILO")
set_property(CACHE PAPILO PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${PAPILO}" PAPILO)

set(READLINE AUTO CACHE STRING "should SCIP link to readline")
set_property(CACHE READLINE PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${READLINE}" READLINE)

option(SANITIZE "should sanitizers be enabled if available" OFF)

set(SYM snauty CACHE STRING "options for symmetry computation")
set_property(CACHE SYM PROPERTY STRINGS nauty snauty bliss sbliss dejavu none)
string(TOLOWER "${SYM}" SYM)

option(THREADSAFE "should SCIP be compiled thread safe" ON)

set(TPI AUTO CACHE STRING "options for thread support library")
set_property(CACHE TPI PROPERTY STRINGS AUTO omp tny none)
string(TOLOWER "${TPI}" TPI)

set(WORHP AUTO CACHE STRING "should SCIP link to WORHP")
set_property(CACHE WORHP PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${WORHP}" WORHP)

set(ZIMPL AUTO CACHE STRING "should SCIP link to ZIMPL")
set_property(CACHE ZIMPL PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${ZIMPL}" ZIMPL)

set(ZLIB AUTO CACHE STRING "should SCIP link to ZLIB")
set_property(CACHE ZLIB PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${ZLIB}" ZLIB)

##### BASIC SETTINGS #####

# Build type can be set using SHARED or, if SHARED=auto, then via BUILD_SHARED_LIBS, the common CMake way.
if(NOT SHARED)
   set(BUILD_SHARED_LIBS 0)
elseif(NOT (SHARED STREQUAL "AUTO"))
   set(BUILD_SHARED_LIBS 1)
endif()
if(BUILD_SHARED_LIBS)
   message(STATUS "Building shared libraries.")
   set(SCIP_HAS_PIC 1)
else()
   message(STATUS "Building static libraries.")
   set(SCIP_HAS_PIC 0)
endif()

# Default build type is 'Release'.
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
elseif(NOT CMAKE_BUILD_TYPE MATCHES "Release|Debug|RelWithDebInfo|MinSizeRel")
   message(FATAL_ERROR "Option CMAKE_BUILD_TYPE has invalid value \"${CMAKE_BUILD_TYPE}\"! Possible are `Release`, `Debug`, `RelWithDebInfo` and `MinSizeRel`.")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Link-time optimization
if(LTO)
   include(CheckIPOSupported)
   check_ipo_supported(RESULT SCIP_HAS_LTO
      LANGUAGES C
      OUTPUT lto_error)
   if(SCIP_HAS_LTO)
      message(STATUS "LTO is enabled.")

      # with link-time optimization enabled, we get additional warnings at link time on
      # - some presumed malloc() calls with size 18446744073709551608 ((uint)-1)
      # - some presumed memmove() calls in boost::container::short_vector via papilo to a 0-length memory
      # we disable these warnings (need to pass them to the compiler when linking)
      if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "CLANG"))
         add_link_options(-Wno-alloc-size-larger-than -Wno-stringop-overflow)
      endif()
   elseif(LTO STREQUAL "AUTO")
      message(STATUS "LTO is not available: ${lto_error}")
   else()
      message(FATAL_ERROR "LTO is not available: ${lto_error}")
   endif()
else()
   message(STATUS "LTO is disabled.")
   set(SCIP_HAS_LTO 0)
endif()

# Query git hash and create a target that creates a preprocessor macro for it.
include(GitHash)
githash_query(SCIP_GITHASH "${PROJECT_SOURCE_DIR}/src/scip/githash.c" "SCIP_GITHASH" ${GITHASH_BRANCH})
githash_add_target(scip_update_githash "${PROJECT_SOURCE_DIR}/src/scip/githash.c" "SCIP_GITHASH" ${GITHASH_BRANCH})
if(SCIP_GITHASH)
   message(STATUS "Githash: ${SCIP_GITHASH}")
endif()


##### DEPENDENCIES #####

# Initialize SCIP_HAS_EXACTSOLVE. we disable it later if dependencies are missing.
if(EXACTSOLVE)
   set(SCIP_HAS_EXACTSOLVE 1)
else()
   set(SCIP_HAS_EXACTSOLVE 0)
endif()

# Search for GMP.
set(SCIP_HAS_GMP 0)
if(GMP)
   if(GMP STREQUAL "AUTO")
      find_package(GMP)
   else()
      find_package(GMP REQUIRED)
   endif()
   if(GMP_FOUND)
      set(SCIP_HAS_GMP 1)
   endif()
else()
   message(STATUS "GMP is disabled.")
endif()

# Support legacy FindGMP.cmake versions.
if(SCIP_HAS_GMP AND NOT TARGET GMP::GMP)
   if(TARGET GMP::gmp)
      message(WARNING "Found GMP version ${GMP_VERSION} has no GMP::GMP target; interfacing GMP::gmp target instead.")
      add_library(GMP::GMP INTERFACE IMPORTED)
      target_link_libraries(GMP::GMP INTERFACE GMP::gmp)
   elseif(TARGET GMP)
      message(WARNING "Found GMP version ${GMP_VERSION} has no GMP::GMP target; interfacing GMP target instead.")
      add_library(GMP::GMP INTERFACE IMPORTED)
      target_link_libraries(GMP::GMP INTERFACE GMP)
   elseif(TARGET gmp)
      message(WARNING "Found GMP version ${GMP_VERSION} has no GMP::GMP target; interfacing gmp target instead.")
      add_library(GMP::GMP INTERFACE IMPORTED)
      target_link_libraries(GMP::GMP INTERFACE gmp)
   elseif(GMP_INCLUDE_DIRS AND GMP_LIBRARIES)
      message(WARNING "Found GMP version ${GMP_VERSION} has no GMP::GMP target; using legacy variables GMP_INCLUDE_DIRS=\"${GMP_INCLUDE_DIRS}\" and GMP_LIBRARIES=\"${GMP_LIBRARIES}\" instead.")
      # GMP_LIBRARIES may be a list, so we need to create each library seperately and glue them via an INTERFACE library.
      list(LENGTH GMP_LIBRARIES num_gmp_lib)
      math(EXPR last_gmp_lib ${num_gmp_lib}-1)
      add_library(GMP_libgmp INTERFACE IMPORTED)
      foreach(id_gmp_lib RANGE ${last_gmp_lib})
         add_library(GMP::libgmp_${id_gmp_lib} UNKNOWN IMPORTED)
         list(GET GMP_LIBRARIES ${id_gmp_lib} the_gmp_lib)
         set_target_properties(GMP::libgmp_${id_gmp_lib} PROPERTIES
            IMPORTED_LOCATION "${the_gmp_lib}"
            INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIRS}")
         target_link_libraries(GMP_libgmp INTERFACE GMP::libgmp_${id_gmp_lib})
      endforeach()
      add_library(GMP::GMP INTERFACE IMPORTED)
      target_link_libraries(GMP::GMP INTERFACE GMP_libgmp)
   else()
      message(FATAL_ERROR "Found GMP version ${GMP_VERSION} has no GMP::GMP target; unable to recover!")
   endif()
endif()

if(NOT SCIP_HAS_GMP AND SCIP_HAS_EXACTSOLVE)
   if(EXACTSOLVE STREQUAL "AUTO")
      message(WARNING "EXACTSOLVE is disabled because it requires GMP.")
      set(SCIP_HAS_EXACTSOLVE 0)
   else()
      message(FATAL_ERROR "EXACTSOLVE requires GMP.")
   endif()
endif()

# Options for memory management.
if(NOBLKBUFMEM)
   set(NOBLKMEM 1 FORCE)
   set(NOBUFMEM 1 FORCE)
endif()

if(NOBLKMEM)
   set(BMS_HAS_NOBLOCKMEM 1)
   message(STATUS "Block memory is disabled.")
else()
   set(BMS_HAS_NOBLOCKMEM 0)
endif()

if(NOBUFMEM)
   set(SCIP_HAS_NOBUFFERMEM 1)
   message(STATUS "Buffer memory is disabled.")
else()
   set(SCIP_HAS_NOBUFFERMEM 0)
endif()

if(DEBUGSOL)
   set(SCIP_HAS_DEBUG_SOLUTION 1)
   message(STATUS "Debug solutions are enabled.")
else()
   set(SCIP_HAS_DEBUG_SOLUTION 0)
endif()

# Thread-safety.
if(THREADSAFE)
   set(SCIP_IS_THREADSAFE 1)
else()
   set(SCIP_IS_THREADSAFE 0)
endif()

# Thread processing interface.
set(SCIP_HAS_TPI 0)
set(SCIP_TPI_NONE 0)
set(SCIP_TPI_OMP 0)
set(SCIP_TPI_TNY 0)

# Search for Threads and use it via TinyCThread.
if((TPI STREQUAL "tny") OR (NOT SCIP_HAS_TPI AND (TPI STREQUAL "auto")))
   if(TPI STREQUAL "tny")
      find_package(Threads REQUIRED)
   else()
      find_package(Threads)
   endif()
   if(Threads_FOUND)
      set(SCIP_HAS_TPI 1)
      set(SCIP_TPI_TNY 1)
      set(SCIP_IS_THREADSAFE 1)
      message(STATUS "Task Processing Interface: TinyCThread")
   endif()
endif()

# Search for OpenMP.
if((TPI STREQUAL "omp") OR (NOT SCIP_HAS_TPI AND (TPI STREQUAL "auto")))
   set(OpenMP_RUNTIME_MSVC "llvm")
   if(TPI STREQUAL "omp")
      find_package(OpenMP REQUIRED COMPONENTS C CXX)
   else()
      find_package(OpenMP COMPONENTS C CXX)
   endif()
   if(OpenMP_C_FOUND)
      set(SCIP_HAS_TPI 1)
      set(SCIP_TPI_OMP 1)
      set(SCIP_IS_THREADSAFE 1)
      message(STATUS "Task Processing Interface: OpenMP")
   endif()
endif()

# No thread library.
if((TPI STREQUAL "none") OR (NOT SCIP_HAS_TPI AND (TPI STREQUAL "auto")))
   set(SCIP_TPI_NONE 1)
   message(STATUS "Task Processing Interface: none")
elseif(NOT SCIP_HAS_TPI)
   message(FATAL_ERROR "Option TPI has invalid value \"${TPI}\"! Possible are `auto`, `omp`, `tny` and `none`.")
endif()

# Search for ZLIB.
set(SCIP_HAS_ZLIB 0)
if(ZLIB)
   message(STATUS "Searching for ZLIB.")
   if(ZLIB STREQUAL "AUTO")
      find_package(ZLIB)
   else()
      find_package(ZLIB REQUIRED)
   endif()
   if(ZLIB_FOUND)
      set(SCIP_HAS_ZLIB 1)
      # Workaround for macOS SDK >= 26.4 not happy about an -isystem /path/to/macOS26.4.sdk/usr/include without
      # a prior -isystem for C++ system headers.
      # Removing the SYSTEM property for ZLIB (which is part of the macOS SDK) changes the -isystem to an -I
      if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
         set_property(TARGET ZLIB::ZLIB PROPERTY SYSTEM FALSE)
      endif()
   endif()
else()
   message(STATUS "ZLIB is disabled.")
endif()

# Search for Readline and Curses.
set(SCIP_HAS_READLINE 0)
set(SCIP_HAS_CURSES 0)
if(READLINE)
   if(READLINE STREQUAL "AUTO")
      find_package(Readline)
   else()
      find_package(Readline REQUIRED)
   endif()
   if(Readline_FOUND)
      set(SCIP_HAS_READLINE 1)
      message(STATUS "Searching for Curses.")
      find_package(Curses)
      if(Curses_FOUND)
         if(NOT TARGET Curses::Curses)
            add_library(Curses::Curses INTERFACE IMPORTED)
            set_target_properties(Curses::Curses
               PROPERTIES
                  INTERFACE_LINK_LIBRARIES "${CURSES_LIBRARIES}"
                  INTERFACE_INCLUDE_DIRECTORIES "${CURSES_INCLUDE_DIRS}"
               )
         endif()
         set(SCIP_HAS_CURSES 1)
      endif()
   endif()
else()
   message(STATUS "Readline is disabled.")
endif()

# Search for MPFR.
set(SCIP_HAS_MPFR 0)
if(MPFR)
   if(MPFR STREQUAL "AUTO")
      find_package(MPFR)
   else()
      find_package(MPFR REQUIRED)
   endif()
   if(MPFR_FOUND)
      set(SCIP_HAS_MPFR 1)
   endif()
else()
   message(STATUS "MPFR is disabled.")
endif()

# Support legacy FindMPFR.cmake versions.
if(SCIP_HAS_MPFR AND NOT TARGET MPFR::MPFR)
   if(TARGET MPFR::mpfr)
      message(WARNING "Found MPFR version ${MPFR_VERSION} has no MPFR::MPFR target; interfacing MPFR::mpfr instead.")
      add_library(MPFR::MPFR INTERFACE IMPORTED)
      target_link_libraries(MPFR::MPFR INTERFACE MPFR::mpfr)
   elseif(TARGET MPFR)
      message(WARNING "Found MPFR version ${MPFR_VERSION} has no MPFR::MPFR target; interfacing MPFR instead.")
      add_library(MPFR::MPFR INTERFACE IMPORTED)
      target_link_libraries(MPFR::MPFR INTERFACE MPFR)
   elseif(TARGET mpfr)
      message(WARNING "Found MPFR version ${MPFR_VERSION} has no MPFR::MPFR target; interfacing mpfr instead.")
      add_library(MPFR::MPFR INTERFACE IMPORTED)
      target_link_libraries(MPFR::MPFR INTERFACE mpfr)
   elseif(MPFR_INCLUDE_DIRS AND MPFR_LIBRARIES)
      message(WARNING "Found MPFR version ${MPFR_VERSION} has no MPFR::MPFR target; using legacy variables MPFR_INCLUDE_DIRS=\"${MPFR_INCLUDE_DIRS}\" and MPFR_LIBRARIES=\"${MPFR_LIBRARIES}\" instead.")
      # MPFR_LIBRARIES may be a list, so we need to create each library seperately and glue them via an INTERFACE library.
      list(LENGTH MPFR_LIBRARIES num_mpfr_lib)
      math(EXPR last_mpfr_lib ${num_mpfr_lib}-1)
      add_library(MPFR_libmpfr INTERFACE IMPORTED)
      foreach(id_mpfr_lib RANGE ${last_mpfr_lib})
         add_library(MPFR::libmpfr_${id_mpfr_lib} UNKNOWN IMPORTED)
         list(GET MPFR_LIBRARIES ${id_mpfr_lib} the_mpfr_lib)
         set_target_properties(MPFR::libmpfr_${id_mpfr_lib} PROPERTIES
            IMPORTED_LOCATION "${the_mpfr_lib}"
            INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE_DIRS}")
         target_link_libraries(MPFR_libmpfr INTERFACE MPFR::libmpfr_${id_mpfr_lib})
      endforeach()
      add_library(MPFR::MPFR INTERFACE IMPORTED)
      target_link_libraries(MPFR::MPFR INTERFACE MPFR_libmpfr)
   else()
      message(FATAL_ERROR "Found MPFR version ${MPFR_VERSION} has no MPFR::MPFR target; unable to recover!")
   endif()
endif()

# No MPFR but EXACTSOLVE.
if(NOT SCIP_HAS_MPFR AND SCIP_HAS_EXACTSOLVE)
   if(EXACTSOLVE STREQUAL "AUTO")
      message(WARNING "Disabling EXACTSOLVE since MPFR was not found.")
      set(SCIP_HAS_EXACTSOLVE 0)
   else()
      message(FATAL_ERROR "MPFR is required for EXACTSOLVE.")
   endif()
endif()

# Search for PaPILO.
set(SCIP_HAS_PAPILO 0)
if(PAPILO)
   if(PAPILO STREQUAL "AUTO")
      find_package(PAPILO ${PAPILO_MIN_VERSION} CONFIG HINTS "${PAPILO_DIR}")
   else()
      find_package(PAPILO ${PAPILO_MIN_VERSION} REQUIRED CONFIG HINTS "${PAPILO_DIR}")
   endif()
   if(PAPILO_FOUND)
      set(SCIP_HAS_PAPILO 1)
   endif()
else()
   message(STATUS "PaPILO is disabled.")
endif()

# Support for legacy PaPILO version.
if(SCIP_HAS_PAPILO AND NOT TARGET PAPILO::PAPILO)
   if(TARGET papilo)
      message(WARNING "Found PaPILO version ${PAPILO_VERSION} has no PAPILO::PAPILO target; interfacing papilo target instead.")
      add_library(PAPILO::PAPILO INTERFACE IMPORTED)
      target_link_libraries(PAPILO::PAPILO INTERFACE papilo)
   else()
      message(FATAL_ERROR "Found PaPILO version ${PAPILO_VERSION} has no PAPILO::PAPILO target; unable to recover!")
   endif()
endif()

# Check if PaPILO is compiled with position-independent code.
if(BUILD_SHARED_LIBS AND SCIP_HAS_PAPILO AND TARGET PAPILO::libpapilo-core)
   get_target_property(papilo_has_pic PAPILO::libpapilo-core POSITION_INDEPENDENT_CODE)
   if(NOT papilo_has_pic)
      message(FATAL_ERROR "Cannot build shared SCIP library with PaPILO's libpapilo-core built without position-independent code.")
   endif()
   if(PAPILO_HAS_LUSOL)
      get_target_property(papilo_has_pic PAPILO::libclusol POSITION_INDEPENDENT_CODE)
      if(NOT papilo_has_pic)
         message(FATAL_ERROR "Cannot build shared SCIP library with PaPILO's libclusol built without position-independent code.")
      endif()
   endif()
endif()

# Search the selected floating-point LP solver library.
if(LPS STREQUAL "auto")
   message(STATUS "Searching for any LP solver.")
elseif(NOT (LPS STREQUAL "none"))
   message(STATUS "Searching for LP solver \"${LPS}\".")
endif()
set(SCIP_HAS_LPS 0)
set(SCIP_LPS_NAME)
set(SCIP_LPS_SPX 0)
set(SCIP_HAS_SOPLEX 0)
set(SCIP_LPS_CLP 0)
set(SCIP_HAS_CLP 0)
set(SCIP_LPS_QSO 0)
set(SCIP_HAS_QSOPTEX 0)
set(SCIP_LPS_GLOP 0)
set(SCIP_HAS_ORTOOLS 0)
set(SCIP_LPS_HIGHS 0)
set(SCIP_HAS_HIGHS 0)
set(SCIP_LPS_CPX 0)
set(SCIP_HAS_CPLEX 0)
set(SCIP_LPS_XPRS 0)
set(SCIP_HAS_XPRESS 0)
set(SCIP_LPS_MSK 0)
set(SCIP_HAS_MOSEK 0)
set(SCIP_LPS_GRB 0)
set(SCIP_HAS_GUROBI 0)
set(SCIP_LPS_NONE 0)
set(SCIP_HAS_LPSCHECK 0)

# Search for SoPlex.
if((LPS STREQUAL "spx") OR (NOT SCIP_HAS_LPS AND (LPS STREQUAL "auto")))
   if(LPS STREQUAL "spx")
      find_package(SOPLEX REQUIRED CONFIG HINTS "${SOPLEX_DIR}")
   else()
      find_package(SOPLEX CONFIG HINTS "${SOPLEX_DIR}")
   endif()
   if(SOPLEX_FOUND)
      set(SCIP_HAS_SOPLEX 1)
      set(SCIP_HAS_LPS 1)
      set(SCIP_LPS_SPX 1)
      set(SCIP_LPS_NAME "SoPlex")
   endif()
   if(LPSCHECK)
      find_package(CPLEX REQUIRED)
      if(NOT CPLEX_FOUND)
         message(FATAL_ERROR "Requested LP check solver CPLEX not found.")
      endif()
      set(SCIP_HAS_LPSCHECK 1)
   endif()
endif()

# Search for CLP.
if((LPS STREQUAL "clp") OR (NOT SCIP_HAS_LPS AND (LPS STREQUAL "auto")))
   if(LPS STREQUAL "clp")
      find_package(CLP REQUIRED)
   else()
      find_package(CLP)
   endif()
   if(CLP_FOUND)
      set(SCIP_HAS_CLP 1)
      set(SCIP_HAS_LPS 1)
      set(SCIP_LPS_CLP 1)
      set(SCIP_LPS_NAME "CLP")
   endif()
endif()

# Search for QSopt_ex.
if((LPS STREQUAL "qso") OR (NOT SCIP_HAS_LPS AND (LPS STREQUAL "auto")))
   if(LPS STREQUAL "qso")
      find_package(QSoptex REQUIRED)
   else()
      find_package(QSoptex)
   endif()
   if(QSoptex_FOUND)
      set(SCIP_HAS_QSOPTEX 1)
      set(SCIP_HAS_LPS 1)
      set(SCIP_LPS_QSO 1)
      set(SCIP_LPS_NAME "QSopt_ex")
   endif()
endif()

# Search for OR-Tools.
if((LPS STREQUAL "glop") OR (NOT SCIP_HAS_LPS AND (LPS STREQUAL "auto") AND SCIP_TPI_TNY))
   if(NOT SCIP_TPI_TNY)
      message(FATAL_ERROR "When using GLOP as an LP solver, TPI needs to be set to 'tny'.")
   endif()
   if(ortools_FOUND)
      message(STATUS "Already found OR-Tools version ${ortools_VERSION}.")
   else()
      message(STATUS "Searching for OR-Tools.")
      # Indent messages for dependency search.
      list(APPEND CMAKE_MESSAGE_INDENT "  ")
      set(OLD_CMAKE_FIND_PACKAGE_PREFER_CONFIG ${CMAKE_FIND_PACKAGE_PREFER_CONFIG})
      set(CMAKE_FIND_PACKAGE_PREFER_CONFIG 1)
      if(LPS STREQUAL "glop")
         find_package(ortools REQUIRED CONFIG HINTS "${ORTOOLS_DIR}")
      else()
         find_package(ortools CONFIG HINTS "${ORTOOLS_DIR}")
      endif()
      set(CMAKE_FIND_PACKAGE_PREFER_CONFIG ${OLD_CMAKE_FIND_PACKAGE_PREFER_CONFIG})

      # Unindent messages again; equivalent to list(POP_BACK, CMAKE_MESSAGE_INDENT) which requires CMake 3.15.
      list(REVERSE CMAKE_MESSAGE_INDENT)
      list(REMOVE_AT CMAKE_MESSAGE_INDENT 0)
      list(REVERSE CMAKE_MESSAGE_INDENT)

      if(ortools_FOUND)
         message(STATUS "Found OR-Tools version ${ortools_VERSION}: ${ortools_CONFIG}")
      endif()
   endif()
   if(ortools_FOUND)
      get_target_property(ortools_pic ortools::ortools INTERFACE_POSITION_INDEPENDENT_CODE)
      if(NOT (ortools_pic STREQUAL "ortools_pic-NOTFOUND"))
         if(ortools_pic AND NOT BUILD_SHARED_LIBS)
            message(WARNING "Disabling OR-Tools since it requires position-independent code and we build a static library.")
         elseif(NOT ortools_pic AND BUILD_SHARED_LIBS)
            message(WARNING "Disabling OR-Tools since it has position-dependent code and we build a shared library.")
         else()
            set(ortools_pic "ortools_pic-NOTFOUND")
         endif()
      endif()
      if(ortools_pic STREQUAL "ortools_pic-NOTFOUND")
         set(SCIP_HAS_ORTOOLS 1)
         set(SCIP_HAS_LPS 1)
         set(SCIP_LPS_GLOP 1)
         set(SCIP_LPS_NAME "GLOP")
      elseif(LPS STREQUAL "glop")
         message(FATAL_ERROR "Choose a different OR-Tools build or switch static/shared library!")
      endif()
   else()
      message(STATUS "Could not find OR-Tools.")
   endif()
endif()

# Search for Highs.
if((LPS STREQUAL "highs") OR (NOT SCIP_HAS_LPS AND (LPS STREQUAL "auto")))
   if(HIGHS_FOUND)
      message(STATUS "Already found HiGHS version ${HIGHS_VERSION}.")
   else()
      message(STATUS "Searching for HiGHS.")
      if(LPS STREQUAL "highs")
         find_package(HIGHS REQUIRED CONFIG HINTS "${HIGHS_DIR}")
      else()
         find_package(HIGHS CONFIG HINTS "${HIGHS_DIR}")
      endif()
      if(HIGHS_FOUND)
         message(STATUS "Found HiGHS version ${HIGHS_VERSION}: ${HIGHS_CONFIG}")
      endif()
   endif()
   if(HIGHS_FOUND)
      get_target_property(highs_pic highs::highs INTERFACE_POSITION_INDEPENDENT_CODE)
      if(NOT (highs_pic STREQUAL "highs_pic-NOTFOUND"))
         if(highs_pic AND NOT BUILD_SHARED_LIBS)
            message(WARNING "Disabling HiGHS since it requires position-independent code and we build a static library.")
         elseif(NOT highs_pic AND BUILD_SHARED_LIBS)
            message(WARNING "Disabling HiGHS since it has position-dependent code and we build a shared library.")
         else()
            set(highs_pic "highs_pic-NOTFOUND")
         endif()
      endif()
      if(highs_pic STREQUAL "highs_pic-NOTFOUND")
         set(SCIP_HAS_HIGHS 1)
         set(SCIP_HAS_LPS 1)
         set(SCIP_LPS_HIGHS 1)
         set(SCIP_LPS_NAME "HiGHS")
         message(WARNING "EXPERIMENTAL: You are using the experimental HiGHS interface. Please report any issues on https://github.com/scipopt/scip.")
      elseif(LPS STREQUAL "highs")
         message(FATAL_ERROR "Choose a different HiGHS build or switch static/shared library!")
      endif()
   endif()
endif()

# Search for CPLEX.
if((LPS STREQUAL "cpx") OR (NOT SCIP_HAS_LPS AND (LPS STREQUAL "auto")))
   if(LPS STREQUAL "cpx")
      find_package(CPLEX REQUIRED)
   else()
      find_package(CPLEX)
   endif()
   if(CPLEX_FOUND)
      set(SCIP_HAS_CPLEX 1)
      set(SCIP_HAS_LPS 1)
      set(SCIP_LPS_CPX 1)
      set(SCIP_LPS_NAME "CPLEX")
   endif()
endif()

# Search for Xpress.
if((LPS STREQUAL "xprs") OR (NOT SCIP_HAS_LPS AND (LPS STREQUAL "auto")))
   if(LPS STREQUAL "xprs")
      find_package(Xpress REQUIRED)
   else()
      find_package(Xpress)
   endif()
   if(Xpress_FOUND)
      set(SCIP_HAS_XPRESS 1)
      set(SCIP_HAS_LPS 1)
      set(SCIP_LPS_XPRS 1)
      set(SCIP_LPS_NAME "Xpress")
   endif()
endif()

# Search for MOSEK.
if((LPS STREQUAL "msk") OR (NOT SCIP_HAS_LPS AND (LPS STREQUAL "auto")))
   if(LPS STREQUAL "msk")
      find_package(MOSEK REQUIRED)
   else()
      find_package(MOSEK)
   endif()
   if(MOSEK_FOUND)
      set(SCIP_HAS_MOSEK 1)
      set(SCIP_HAS_LPS 1)
      set(SCIP_LPS_MSK 1)
      set(SCIP_LPS_NAME "MOSEK")
   endif()
endif()

# Search for Gurobi.
if((LPS STREQUAL "grb") OR (NOT SCIP_HAS_LPS AND (LPS STREQUAL "auto")))
   if(LPS STREQUAL "grb")
      find_package(Gurobi REQUIRED)
   else()
      find_package(Gurobi)
   endif()
   if(Gurobi_FOUND)
      set(SCIP_HAS_GUROBI 1)
      set(SCIP_HAS_LPS 1)
      set(SCIP_LPS_GRB 1)
      set(SCIP_LPS_NAME "Gurobi")
   endif()
endif()

# No LP solver found.
if(NOT SCIP_HAS_LPS)
   if(LPS STREQUAL "auto")
      message(WARNING "No LP solver found!")
   elseif(NOT (LPS STREQUAL "none"))
      message(FATAL_ERROR "Option LPS has invalid value \"${LPS}\"! Possible are `auto`, `spx`, `clp`, `qso`, `glop`, `highs`, `cpx`, `xprs`, `msk`, `grb`, `none`.")
   endif()
   set(SCIP_LPS_NONE 1)
endif()

if(SCIP_LPS_NONE)
   message(STATUS "Building SCIP without floating-point LP solver.")
else()
   message(STATUS "Building SCIP with floating-point LP solver ${SCIP_LPS_NAME}.")
endif()

# Boost
set(SCIP_HAS_BOOST 0)
if(SCIP_HAS_EXACTSOLVE)
   if(Boost_FOUND AND ("${Boost_VERSION}" VERSION_GREATER_EQUAL ${BOOST_MIN_VERSION}))
      message(STATUS "Already found Boost version ${Boost_VERSION_STRING}, minimum required is ${BOOST_MIN_VERSION}.")
   elseif(EXACTSOLVE STREQUAL "AUTO")
      message(STATUS "Searching for Boost.")
      find_package(Boost ${BOOST_MIN_VERSION})
   else()
      message(STATUS "Searching for Boost.")
      find_package(Boost ${BOOST_MIN_VERSION} REQUIRED)
   endif()
   if(Boost_FOUND)
      if(NOT DEFINED Boost_VERSION_STRING)
         math(EXPR _boost_major "${Boost_VERSION} / 100000")
         math(EXPR _boost_minor "${Boost_VERSION} / 100 % 1000")
         math(EXPR _boost_patch "${Boost_VERSION} % 100")
         set(Boost_VERSION_STRING "${_boost_major}.${_boost_minor}.${_boost_patch}")
      endif()
      message(STATUS "Found Boost version ${Boost_VERSION_STRING}.")
   endif()
   if(Boost_FOUND AND (Boost_VERSION_STRING VERSION_GREATER_EQUAL ${BOOST_MIN_VERSION}))
      set(SCIP_HAS_BOOST 1)
   else()
      set(SCIP_HAS_EXACTSOLVE 0)
      if(EXACTSOLVE AND NOT (EXACTSOLVE STREQUAL "AUTO")) # EXACTSOLVE was explicitly requested.
         message(FATAL_ERROR "Boost version ${BOOST_MIN_VERSION} or higher not found, which is required for EXACTSOLVE; try specifying Boost_DIR.")
      endif()
   endif()
endif()

# Find exact LP solvers.
set(SCIP_HAS_LPSEXACT 0)
set(SCIP_LPSEXACT_NAME)
set(SCIP_LPSEXACT_SPX 0)
set(SCIP_LPSEXACT_QSOEX 0)
set(SCIP_LPSEXACT_NONE 0)

# Search for exact SoPlex.
if(SCIP_HAS_EXACTSOLVE AND ((LPSEXACT STREQUAL "spx") OR (NOT SCIP_HAS_LPSEXACT AND (LPSEXACT STREQUAL "auto"))))
   set(bad_soplex_reason 0)
   # TODO: Remove SOPLEX_WITH_MPFR once SoPlex versions <= 8 are not compatible anymore.
   if(SOPLEX_FOUND AND (SOPLEX_HAS_MPFR OR SOPLEX_WITH_MPFR))
      message(STATUS "Already found LP solver SoPlex ${SOPLEX_VERSION}.")
   elseif(SOPLEX_FOUND)
      set(bad_soplex_reason "does not support MPFR")
   else()
      message(STATUS "Searching exact LP solver SoPlex.")
      if(EXACTSOLVE AND (NOT EXACTSOLVE STREQUAL "AUTO") AND (LPSEXACT STREQUAL "spx"))
         find_package(SOPLEX CONFIG REQUIRED HINTS "${SOPLEX_DIR}")
      else()
         find_package(SOPLEX CONFIG HINTS "${SOPLEX_DIR}")
      endif()
   endif()
   if(NOT SOPLEX_FOUND)
      set(bad_soplex_reason "not found")
   elseif(NOT SOPLEX_HAS_MPFR AND NOT SOPLEX_WITH_MPFR)
      set(bad_soplex_reason "does not support MPFR")
   endif()
   if(SOPLEX_FOUND AND (SOPLEX_HAS_MPFR OR SOPLEX_WITH_MPFR))
      set(SCIP_HAS_LPSEXACT 1)
      set(SCIP_LPSEXACT_SPX 1)
      set(SCIP_LPSEXACT_NAME "SoPlex")
      set(SCIP_HAS_SOPLEX 1)
   elseif(LPSEXACT STREQUAL "spx" AND NOT (EXACTSOLVE STREQUAL "AUTO")) # SoPlex and EXACTSOLVE were explicitly requested.
      message(FATAL_ERROR "LP solver SoPlex ${bad_soplex_reason}, but is required for EXACTSOLVE; try specifying suitable SOPLEX_DIR.")
   elseif(LPSEXACT STREQUAL "spx") # SoPlex was explicitly requested, but EXACTSOLVE is AUTO.
      message(WARNING "LP solver SoPlex ${bad_soplex_reason}, but is required for EXACTSOLVE; try specifying suitable SOPLEX_DIR.")
      set(SCIP_HAS_EXACTSOLVE 0)
   else()
      message(STATUS "LP solver SoPlex ${bad_soplex_reason}.")
   endif()
endif()

# Search for QSopt_ex.
if(SCIP_HAS_EXACTSOLVE AND ((LPSEXACT STREQUAL "qsoex") OR (NOT SCIP_HAS_LPSEXACT AND (LPSEXACT STREQUAL "auto"))))
   message(STATUS "Searching exact LP solver QSopt_ex.")
   if(QSoptex_FOUND)
      message(STATUS "Already found LP solver QSopt_ex ${QSoptex_VERSION}.")
   elseif(EXACTSOLVE AND (NOT EXACTSOLVE STREQUAL "AUTO") AND (LPSEXACT STREQUAL "qsoex"))
      find_package(QSoptex REQUIRED)
   else()
      find_package(QSoptex)
   endif()
   if(QSoptex_FOUND)
      set(SCIP_HAS_LPSEXACT 1)
      set(SCIP_LPSEXACT_QSOEX 1)
      set(SCIP_LPSEXACT_NAME "QSopt_ex")
      set(SCIP_HAS_QSOPTEX 1)
   elseif(LPSEXACT STREQUAL "qsoex" AND NOT (EXACTSOLVE STREQUAL "AUTO")) # QSopt_ex and EXACTSOLVE were explicitly requested.
      message(FATAL_ERROR "Exact LP solver QSopt_ex not found, which is required for EXACTSOLVE; try specifying QSOPTEX_DIR.")
   elseif(LPSEXACT STREQUAL "qsoex") # QSopt_ex was explicitly requested, but EXACTSOLVE is AUTO.
      message(WARNING "Exact LP solver QSopt_ex not found, which is required for EXACTSOLVE; try specifying QSOPTEX_DIR.")
      set(SCIP_HAS_EXACTSOLVE 0)
   endif()
endif()

# No exact LP solver found.
if(NOT SCIP_HAS_LPSEXACT)
   # We tried to find one though.
   if(SCIP_HAS_EXACTSOLVE)
      if(LPSEXACT STREQUAL "auto" AND EXACTSOLVE STREQUAL "AUTO")
         message(STATUS "No exact LP solver found, which is required for EXACTSOLVE.")
      elseif(LPSEXACT STREQUAL "auto")
         message(FATAL_ERROR "No exact LP solver found, which is required for EXACTSOLVE.")
      endif()
   endif()
   set(SCIP_LPSEXACT_NONE 1)
   set(SCIP_HAS_EXACTSOLVE 0)
endif()



# No exact LP solver found.
if(NOT SCIP_HAS_EXACTSOLVE AND (EXACTSOLVE STREQUAL "AUTO"))
   # Disable dependencies only required for exact solving.
   if(SCIP_HAS_BOOST)
      message(STATUS "Disabling Boost since EXACTSOLVE is disabled.")
      set(SCIP_HAS_BOOST 0)
   endif()
   if(SCIP_HAS_MPFR)
      message(STATUS "Disabling MPFR since EXACTSOLVE is disabled.")
      set(SCIP_HAS_MPFR 0)
   endif()
endif()

# Now we are sure about EXACTSOLVE.
if(SCIP_HAS_EXACTSOLVE)
   message(STATUS "Building SCIP with support for exact solving mode using exact LP solver ${SCIP_LPSEXACT_NAME}.")

   message(STATUS "Searching for VIPR checker")
   find_program(VIPRCOMP "viprcomp" HINTS ${VIPR_DIR} $ENV{VIPR_DIR})
   find_program(VIPRCHK "viprchk" HINTS ${VIPR_DIR} $ENV{VIPR_DIR})
   if(VIPRCOMP AND VIPRCHK)
      message(STATUS "Found VIPR checker.")
   else()
      message(STATUS "Could NOT find VIPR checker: verification is not enabled for exact ctest. Make sure to install VIPR (https://github.com/scipopt/vipr) to enable certificate checking.")
   endif()

else()
    message(STATUS "Building SCIP without support for exact solving mode.")
    set(lpiexact "lpiexact/lpiexact_none.c")
endif()

# Support for legacy SoPlex version.
if(SCIP_HAS_SOPLEX AND NOT TARGET SOPLEX::SOPLEX)
   if(TARGET libsoplex-pic)
      message(WARNING "Found SoPlex version ${SOPLEX_VERSION} has no SOPLEX::SOPLEX target; interfacing libsoplex-pic target instead.")
      add_library(SOPLEX::SOPLEX INTERFACE IMPORTED)
      target_link_libraries(SOPLEX::SOPLEX INTERFACE libsoplex-pic)
   elseif(TARGET libsoplexshared)
      message(WARNING "Found SoPlex version ${SOPLEX_VERSION} has no SOPLEX::SOPLEX target; interfacing libsoplexshared target instead.")
      add_library(SOPLEX::SOPLEX INTERFACE IMPORTED)
      target_link_libraries(SOPLEX::SOPLEX INTERFACE libsoplexshared)
   else()
      message(FATAL_ERROR "Found SoPlex version ${SOPLEX_VERSION} has no SOPLEX::SOPLEX target; unable to recover!")
   endif()
endif()

# Search for ZIMPL.
set(SCIP_HAS_ZIMPL 0)
if(ZIMPL)
   if(ZIMPL STREQUAL "AUTO")
      find_package(ZIMPL ${ZIMPL_MIN_VERSION} CONFIG HINTS "${ZIMPL_DIR}")
   else()
      find_package(ZIMPL ${ZIMPL_MIN_VERSION} REQUIRED CONFIG HINTS "${ZIMPL_DIR}")
   endif()
   if(ZIMPL_FOUND)
      set(SCIP_HAS_ZIMPL 1)
   endif()
else()
   message(STATUS "ZIMPL is disabled.")
endif()

# We only support ZIMPL with GMP.
if(SCIP_HAS_ZIMPL AND DEFINED ZIMPL_HAS_GMP AND NOT ZIMPL_HAS_GMP)
   set(SCIP_HAS_ZIMPL 0)
   if(ZIMPL STREQUAL "AUTO")
      message(STATUS "Disabling ZIMPL since it is compiled without GMP.")
   else()
      message(FATAL_ERROR "ZIMPL is compiled without GMP!")
   endif()
endif()

# Support for legacy ZIMPL version.
if(SCIP_HAS_ZIMPL AND NOT TARGET ZIMPL::ZIMPL)
   if(TARGET libzimpl-pic)
      message(WARNING "Found ZIMPL version ${ZIMPL_VERSION} has no ZIMPL::ZIMPL target; interfacing libzimpl-pic target instead.")
      add_library(ZIMPL::ZIMPL INTERFACE IMPORTED)
      target_link_libraries(ZIMPL::ZIMPL INTERFACE libzimpl-pic)
   else()
      message(FATAL_ERROR "Found ZIMPL version ${ZIMPL_VERSION} has no ZIMPL::ZIMPL target; unable to recover!")
   endif()
endif()

# AMPL interface.
if(AMPL)
  set(SCIP_HAS_AMPL 1)
  message(STATUS "Building AMPL interface.")
else()
  set(SCIP_HAS_AMPL 0)
  message(STATUS "Disabling AMPL interface.")
endif()

# Expression interpreter.
set(SCIP_HAS_EXPRINT 1)
set(SCIP_HAS_CPPAD 0)
set(SCIP_EXPRINT_NAME NONE)
set(SCIP_EXPRINT_VERSION NONE)
if(EXPRINT STREQUAL "cppad")
   set(SCIP_HAS_CPPAD 1)
   set(SCIP_EXPRINT_NAME "CppAD")
   set(SCIP_EXPRINT_VERSION ${CPPAD_SHIPPED_VERSION})
elseif(EXPRINT STREQUAL "none")
   set(SCIP_HAS_EXPRINT 0)
else()
   message(FATAL_ERROR "EXPRINT option has wrong value \"${EXPRINT}\"! Possible are `cppad` and `none`.")
endif()

if(SCIP_HAS_EXPRINT)
   message(STATUS "Expression interpreter: ${SCIP_EXPRINT_NAME} version ${SCIP_EXPRINT_VERSION}.")
else()
   message(STATUS "Expression interpreter disabled.")
endif()

# IPOPT and LAPACK are incompatible.
if(IPOPT AND NOT (IPOPT STREQUAL "AUTO") AND LAPACK AND NOT (LAPACK STREQUAL "AUTO"))
   message(FATAL_ERROR "Ipopt and LAPACK are incompatible!")
endif()

# Search for Ipopt.
set(SCIP_HAS_IPOPT 0)
if(IPOPT AND NOT (LAPACK AND NOT (LAPACK STREQUAL "AUTO")))
   if(IPOPT STREQUAL "AUTO")
      find_package(Ipopt ${IPOPT_MIN_VERSION})
   else()
      find_package(Ipopt ${IPOPT_MIN_VERSION} REQUIRED)
   endif()
   if(Ipopt_FOUND)
      set(SCIP_HAS_IPOPT 1)
   endif()
elseif(IPOPT)
   message(STATUS "Ipopt is disabled due to LAPACK.")
else()
   message(STATUS "Ipopt is disabled.")
endif()

# Search for LAPACK.
set(SCIP_HAS_LAPACK 0)
if(LAPACK AND NOT SCIP_HAS_IPOPT)
   message(STATUS "Searching for LAPACK.")
   list(APPEND CMAKE_MESSAGE_INDENT "  ") # Indent messages for dependency search.
   if(LAPACK STREQUAL "AUTO")
      find_package(LAPACK)
   else()
      find_package(LAPACK REQUIRED)
   endif()
   list(POP_BACK CMAKE_MESSAGE_INDENT) # Unindent again.
   if(LAPACK_FOUND)
      set(SCIP_HAS_LAPACK 1)
   endif()
elseif(LAPACK)
   message(STATUS "LAPACK is disabled due to Ipopt.")
else()
   message(STATUS "LAPACK is disabled.")
endif()

# Search for WORHP.
set(SCIP_HAS_WORHP 0)
if(WORHP)
   if(WORHP STREQUAL "AUTO")
      find_package(WORHP)
   else()
      find_package(WORHP REQUIRED)
   endif()
   if(WORHP_FOUND)
      set(SCIP_HAS_WORHP 1)
   endif()
else()
   message(STATUS "WORHP is disabled.")
endif()

# Search for CONOPT.
set(SCIP_HAS_CONOPT 0)
if(CONOPT)
   if(CONOPT STREQUAL "AUTO")
      find_package(CONOPT)
   else()
      find_package(CONOPT REQUIRED)
   endif()
   if(CONOPT_FOUND)
      set(SCIP_HAS_CONOPT 1)
   endif()
else()
   message(STATUS "CONOPT is disabled.")
endif()

# Symmetry computation tools.
set(SCIP_HAS_BLISS 0)
set(SCIP_HAS_NAUTY 0)
set(SCIP_HAS_EXTERNAL_NAUTY 0)
set(SCIP_HAS_SASSY 0)
set(SCIP_HAS_DEJAVU 0)
set(SCIP_HAS_SYM 1)
set(SCIP_SYM_NAME NONE)
set(SCIP_SYM_VERSION NONE)
if(SYM MATCHES "s?nauty")
   if(NAUTY_DIR)
      find_package(nauty)
   endif()
   if(nauty_FOUND)
      set(SCIP_HAS_EXTERNAL_NAUTY 1)
      set(SCIP_SYM_VERSION "${nauty_VERSION}")
   else()
      message(STATUS "Using nauty version shipped with SCIP; if you prefer a different version, please specify NAUTY_DIR.")
      set(SCIP_SYM_VERSION ${NAUTY_SHIPPED_VERSION})
   endif()
   set(SCIP_HAS_NAUTY 1)
   if(SYM STREQUAL "snauty")
      set(SCIP_HAS_SASSY 1)
      set(SCIP_SYM_NAME "sassy+nauty")
      set(SCIP_SYM_VERSION "${DEJAVU_SHIPPED_VERSION}+${SCIP_SYM_VERSION}")
   else()
      set(SCIP_SYM_NAME "nauty")
   endif()
elseif(SYM MATCHES "s?bliss")
   find_package(Bliss REQUIRED CONFIG HINTS "${BLISS_DIR}")
   set(SCIP_HAS_BLISS 1)
   if(SYM STREQUAL "sbliss")
      set(SCIP_HAS_SASSY 1)
      set(SCIP_SYM_NAME "sassy+bliss")
      set(SCIP_SYM_VERSION "${DEJAVU_SHIPPED_VERSION}+${Bliss_VERSION}")
   else()
      set(SCIP_SYM_NAME "bliss")
      set(SCIP_SYM_VERSION "${Bliss_VERSION}")
   endif()
elseif(SYM STREQUAL "dejavu")
   set(SCIP_HAS_DEJAVU 1)
   set(SCIP_HAS_SASSY 1)
   set(SCIP_SYM_NAME "dejavu")
   set(SCIP_SYM_VERSION ${DEJAVU_SHIPPED_VERSION})
elseif(SYM STREQUAL "none")
   set(SCIP_HAS_SYM 0)
   set(SCIP_SYM_NAME "none")
else()
   message(FATAL_ERROR "Option SYM has invalid value \"${SYM}\"! Possible are `nauty`, `snauty`, `bliss`, `sbliss`, `dejavy` and `none`.")
endif()

if(SCIP_HAS_SYM)
   message(STATUS "Symmetry computation via ${SCIP_SYM_NAME} version ${SCIP_SYM_VERSION}.")
else()
   message(STATUS "Symmetry computation disabled.")
endif()

# Check if dependencies are compiled with position-independent code.
if(BUILD_SHARED_LIBS)
   if(SCIP_HAS_PAPILO AND (DEFINED PAPILO_HAS_PIC) AND (NOT PAPILO_HAS_PIC))
      message(FATAL_ERROR "Cannot link shared SCIP library to a PaPILO library without position-independent code!")
   endif()
   if(SCIP_HAS_SOPLEX AND (DEFINED SOPLEX_HAS_PIC) AND (NOT SOPLEX_HAS_PIC))
      message(FATAL_ERROR "Cannot link shared SCIP library to a SoPlex library without position-independent code!")
   endif()
   if(SCIP_HAS_ZIMPL AND (DEFINED ZIMPL_HAS_PIC) AND (NOT ZIMPL_HAS_PIC))
      message(FATAL_ERROR "Cannot link shared SCIP library to a ZIMPL library without position-independent code!")
   endif()
endif()

# Warn if we build static libraries and activate testing.
if(NOT BUILD_SHARED_LIBS AND BUILD_TESTING)
   message(WARNING "Statically linking libscip produces many large unittest executables!")
endif()

##### TOOLS #####

# Enable coverage support.
#
# It is very important to execute this prior to adding subdirectories because of
# compiler flag changes.
#
if(COVERAGE)
   message(STATUS "CodeCoverage is enabled.")
   include(CodeCoverage)
   append_coverage_compiler_flags()
   set(COVERAGE_EXCLUDES "/usr*")

   # Create a CMake script file that is executed to run the coverage test.
   # With a script file, the return code of the tests is simply ignored,
   # and a coverage report is generated even if some tests fail currently.
   file(WRITE "${PROJECT_BINARY_DIR}/RunCoverage.cmake" "execute_process(COMMAND ctest ${CMAKE_CTEST_COMMAND} ${COVERAGE_CTEST_ARGS})")

   # Setup the coverage target to execute the RunCoverage script.
   setup_target_for_coverage(NAME coverage
      EXECUTABLE ${CMAKE_COMMAND} -P "${PROJECT_BINARY_DIR}/RunCoverage.cmake"
      DEPENDENCIES all_executables
      )
else()
   message(STATUS "CodeCoverage is disabled.")
endif()

# Determine C/C++ standard.
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED TRUE)
set(CMAKE_C_EXTENSIONS OFF)
if(SCIP_HAS_ORTOOLS)
   set(CMAKE_CXX_STANDARD 17)
else()
   set(CMAKE_CXX_STANDARD 14)
endif()

# Set C/C++ flags for Release with debug info configuration.
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}")

# Set defines for Windows
if(WIN32)
    set(SCIP_NO_SIGACTION on)
    set(SCIP_NO_STRTOK_R on)
endif()
set(WITH_SCIPDEF on)

# if changing these flags, also update GCCWARN/GXXWARN in make/make.project
set(ADD_C_AND_CXX_FLAGS -pedantic -Wall -W -Wpointer-arith -Wcast-align -Wwrite-strings -Wshadow  -Wredundant-decls -Wdisabled-optimization -Wno-long-long -Wno-unknown-pragmas -Wno-unused-parameter)
set(ADD_C_FLAGS -Wsign-compare -Wstrict-prototypes -Wmissing-declarations -Wmissing-prototypes -Wno-override-init -Wno-uninitialized -Wno-maybe-uninitialized)
set(ADD_CXX_FLAGS -Wnon-virtual-dtor -Wreorder -Woverloaded-virtual -Wsign-promo -Wsynth -Wcast-qual -Wno-strict-overflow -Wno-psabi)

if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
   # disable fused floating point contraction to enhance reproducibility across compilers and architectures
   add_compile_options(/fp:precise)
   if((LPS STREQUAL "spx") AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
      # add /bigobj to avoid error C1128: number of sections exceeded object file format limit
      # (see also https://github.com/scipopt/scip/pull/200)
      add_compile_options(/bigobj)
   endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
   # disable fused floating point contraction to enhance reproducibility across compilers and architectures
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

if(MSVC)
#   add_definitions(/W4)
   add_definitions(/wd4100)
   add_definitions(/wd4244)
   add_definitions(-D_CRT_SECURE_NO_WARNINGS)
   add_definitions(/Zc:__cplusplus)
   add_definitions(/EHsc)
endif()

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

if(CHECKSTAGE AND
   ((NOT CHECKSTAGE STREQUAL "AUTO") OR (CMAKE_BUILD_TYPE STREQUAL "Debug")))
   if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      message(STATUS "Enabling stage checks.")
   else()
      message(WARNING "Enabling stage checks decreases performance.")
   endif()
   set(SCIP_CHECK_STAGE 1)
else()
   set(SCIP_CHECK_STAGE 0)
endif()

if(CHECKNAME AND
   ((NOT CHECKNAME STREQUAL "AUTO") OR (CMAKE_BUILD_TYPE STREQUAL "Debug")))
   if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      message(STATUS "Enabling name checks.")
   else()
      message(WARNING "Enabling name checks decreases performance.")
   endif()
   set(SCIP_CHECK_NAME 1)
else()
   set(SCIP_CHECK_NAME 0)
endif()

##### SUBDIRECTORIES #####

# Main code
add_subdirectory(src)

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

    #
    # add examples
    #
    # use sub directory bin/examples for executables of examples
    #
    add_subdirectory(examples)

    #
    # add applications
    #
    add_subdirectory(applications)

    add_custom_target(all_executables DEPENDS scip unittests examples applications)
endif()

# Documentation
add_subdirectory(doc EXCLUDE_FROM_ALL)

# Feature summary
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    include(FeatureSummary)
    feature_summary(WHAT ALL)
endif()

##### INSTALLATION #####

# Install license files.
install(FILES "${PROJECT_SOURCE_DIR}/LICENSE"
   DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/scip")
install(FILES "${PROJECT_SOURCE_DIR}/src/tclique/LICENSE"
   DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/scip/tclique")

if(SCIP_HAS_AMPL)
   install(FILES "${PROJECT_SOURCE_DIR}/src/amplmp/LICENSE.rst"
      DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/scip/amplmp")
endif()

if(SCIP_TPI_TNY)
   install(FILES "${PROJECT_SOURCE_DIR}/src/tinycthread/COPYRIGHT"
      DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/scip/tinycthread")
endif()

if(SCIP_HAS_CPPAD)
   install(FILES "${PROJECT_SOURCE_DIR}/src/cppad/COPYING"
      DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/scip/cppad")
endif()

# Create scip-config-version.cmake.
include(CMakePackageConfigHelpers)
write_basic_package_version_file("${CMAKE_BINARY_DIR}/scip-config-version.cmake"
   VERSION ${SCIP_VERSION_MAJOR}.${SCIP_VERSION_MINOR}.${SCIP_VERSION_PATCH}
   COMPATIBILITY SameMajorVersion
   )

# Unless we build a package, we populate search directories for dependencies.
include(PackageConfigDependencies)

if(NOT BUILD_PACKAGE AND SCIP_HAS_BOOST AND Boost_DIR)
   package_config_dependency_dir(SCIP_CONFIG_BOOST_DIR Boost_FOUND Boost_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_CLP AND CLP_DIR)
   package_config_dependency_dir(SCIP_CONFIG_CLP_DIR CLP_FOUND CLP_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_CLP AND COINUTILS_DIR)
   package_config_dependency_dir(SCIP_CONFIG_COINUTILS_DIR COINUTILS_FOUND COINUTILS_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_CONOPT AND CONOPT_DIR)
   package_config_dependency_dir(SCIP_CONFIG_CONOPT_DIR CONOPT_FOUND CONOPT_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_CPLEX AND CPLEX_DIR)
   package_config_dependency_dir(SCIP_CONFIG_CPLEX_DIR CPLEX_FOUND CPLEX_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_GMP AND GMP_DIR)
   package_config_dependency_dir(SCIP_CONFIG_GMP_DIR GMP_FOUND GMP_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_GUROBI AND GUROBI_DIR)
   package_config_dependency_dir(SCIP_CONFIG_GUROBI_DIR Gurobi_FOUND GUROBI_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_HIGHS AND HIGHS_DIR)
   package_config_dependency_dir(SCIP_CONFIG_HIGHS_DIR HIGHS_FOUND HIGHS_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_IPOPT AND IPOPT_DIR)
   package_config_dependency_dir(SCIP_CONFIG_IPOPT_DIR Ipopt_FOUND IPOPT_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_MOSEK AND MOSEK_DIR)
   package_config_dependency_dir(SCIP_CONFIG_MOSEK_DIR MOSEK_FOUND MOSEK_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_MPFR AND MPFR_DIR)
   package_config_dependency_dir(SCIP_CONFIG_MPFR_DIR MPFR_FOUND MPFR_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_ORTOOLS AND ORTOOLS_DIR)
   package_config_dependency_dir(SCIP_CONFIG_ORTOOLS_DIR ortools_FOUND ORTOOLS_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_PAPILO AND PAPILO_DIR)
   package_config_dependency_dir(SCIP_CONFIG_PAPILO_DIR PAPILO_FOUND PAPILO_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_QSOPTEX AND QSOPTEX_DIR)
   package_config_dependency_dir(SCIP_CONFIG_QSOPTEX_DIR QSoptex_FOUND QSOPTEX_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_SOPLEX AND SOPLEX_DIR)
   package_config_dependency_dir(SCIP_CONFIG_SOPLEX_DIR SOPLEX_FOUND SOPLEX_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_XPRESS AND XPRESS_DIR)
   package_config_dependency_dir(SCIP_CONFIG_XPRESS_DIR Xpress_FOUND XPRESS_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_WORHP AND WORHP_DIR)
   package_config_dependency_dir(SCIP_CONFIG_WORHP_DIR WORHP_FOUND WORHP_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_ZIMPL AND ZIMPL_DIR)
   package_config_dependency_dir(SCIP_CONFIG_ZIMPL_DIR ZIMPL_FOUND ZIMPL_DIR)
endif()
if(NOT BUILD_PACKAGE AND SCIP_HAS_ZLIB AND ZLIB_DIR)
   package_config_dependency_dir(SCIP_CONFIG_ZLIB_DIR ZLIB_FOUND ZLIB_DIR)
endif()

# Configure scip-config.cmake for the build-tree:
# We add the CMake module path from the sources.
set(EXTRA_CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
configure_file("${PROJECT_SOURCE_DIR}/scip-config.cmake.in"
   "${CMAKE_BINARY_DIR}/scip-config.cmake"
   @ONLY
   )

# Configure scip-config.cmake for the installation.
# We add the current directory of the installed scip-config.cmake.
set(EXTRA_CMAKE_MODULE_PATH "\${CMAKE_CURRENT_LIST_DIR}")
configure_file("${PROJECT_SOURCE_DIR}/scip-config.cmake.in"
   "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/scip-config.cmake"
   @ONLY
   )

# Install CMake config and modules.
install(FILES
   "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/scip-config.cmake"
   "${CMAKE_BINARY_DIR}/scip-config-version.cmake"
   "$<$<BOOL:${SCIP_HAS_CONOPT}>:${PROJECT_SOURCE_DIR}/cmake/FindCONOPT.cmake>"
   "$<$<BOOL:${SCIP_HAS_WORHP}>:${PROJECT_SOURCE_DIR}/cmake/FindWORHP.cmake>"
   "$<$<BOOL:${SCIP_HAS_GMP}>:${PROJECT_SOURCE_DIR}/cmake/FindGMP.cmake>"
   "$<$<BOOL:${SCIP_HAS_MPFR}>:${PROJECT_SOURCE_DIR}/cmake/FindMPFR.cmake>"
   "$<$<BOOL:${SCIP_HAS_EXTERNAL_NAUTY}>:${PROJECT_SOURCE_DIR}/cmake/Findnauty.cmake>"
   "$<$<BOOL:${SCIP_HAS_IPOPT}>:${PROJECT_SOURCE_DIR}/cmake/FindIpopt.cmake>"
   "$<$<BOOL:${SCIP_HAS_CLP}>:${PROJECT_SOURCE_DIR}/cmake/FindCLP.cmake>"
   "$<$<BOOL:${SCIP_HAS_CLP}>:${PROJECT_SOURCE_DIR}/cmake/FindCoinUtils.cmake>"
   "$<$<BOOL:${SCIP_HAS_QSOPTEX}>:${PROJECT_SOURCE_DIR}/cmake/FindQSoptex.cmake>"
   "$<$<BOOL:${SCIP_HAS_READLINE}>:${PROJECT_SOURCE_DIR}/cmake/FindReadline.cmake>"
   "$<$<BOOL:${SCIP_HAS_CPLEX}>:${PROJECT_SOURCE_DIR}/cmake/FindCPLEX.cmake>"
   "$<$<BOOL:${SCIP_HAS_XPRESS}>:${PROJECT_SOURCE_DIR}/cmake/FindXpress.cmake>"
   "$<$<BOOL:${SCIP_HAS_MOSEK}>:${PROJECT_SOURCE_DIR}/cmake/FindMOSEK.cmake>"
   "$<$<BOOL:${SCIP_HAS_GUROBI}>:${PROJECT_SOURCE_DIR}/cmake/FindGurobi.cmake>"
   "${PROJECT_SOURCE_DIR}/cmake/FindPackageMessage.cmake"
   DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/scip"
   )

# Install the scip-targets export set.
install(EXPORT scip-targets
   NAMESPACE SCIP::
   DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/scip"
   )

# Export targets for the build-tree.
export(EXPORT scip-targets
   NAMESPACE SCIP::
   FILE "${CMAKE_BINARY_DIR}/scip-targets.cmake"
   )

##### LOGGING #####

# Write log file of non-default parameter settings
file(WRITE "${CMAKE_BINARY_DIR}/cmake.log" "")

function(log_option option default_value optional_extra)
   string(TOUPPER ${${option}} current_value)
   set(true_values "TRUE;ON;YES;Y;1")
   set(false_values "FALSE;OFF;NO;N;0")
   if(current_value STREQUAL default_value)
      return()
   elseif((default_value IN_LIST true_values) AND (current_value IN_LIST true_values))
      return()
   elseif((default_value IN_LIST false_values) AND (current_value IN_LIST false_values))
      return()
   else()
      file(APPEND "${CMAKE_BINARY_DIR}/cmake.log" " -D${option}=${${option}}")
   endif()
   foreach(extra ${optional_extra})
      if(${extra})
         file(APPEND "${CMAKE_BINARY_DIR}/cmake.log" " -D${extra}=${${extra}}")
      endif()
   endforeach()
endfunction()

message(STATUS "Writing file with modified options to ${CMAKE_BINARY_DIR}/cmake.log.")
log_option(AMPL ON "")
log_option(BUILD_SHARED_LIBS ON "")
log_option(CHECKNAME AUTO "")
log_option(CHECKSTAGE AUTO "")
log_option(CONOPT AUTO "CONOPT_DIR;CONOPT_LICENSE_INT_1;CONOPT_LICENSE_INT_2;CONOPT_LICENSE_INT_3;CONOPT_LICENSE_TEXT")
log_option(COVERAGE OFF "")
log_option(DEBUGSOL OFF "")
log_option(EXACTSOLVE AUTO "")
log_option(EXPRINT CPPAD "")
log_option(GMP AUTO "GMP_DIR")
log_option(IPOPT AUTO "IPOPT_DIR")
log_option(LAPACK AUTO "")
log_option(LPS SPX "SOPLEX_DIR;CLP_DIR;COINUTILS_DIR;CPLEX_DIR;ORTOOLS_DIR;QSOPTEX_DIR;HIGHS_DIR;XPRESS_DIR;MOSEK_DIR;GUROBI_DIR")
log_option(LPSEXACT AUTO "SOPLEX_DIR;QSOPTEX_DIR")
log_option(LPSCHECK OFF "")
log_option(LTO OFF "")
log_option(MPFR AUTO "MPFR_DIR")
log_option(NOBLKBUFMEM OFF "")
log_option(NOBLKMEM OFF "")
log_option(NOBUFMEM OFF "")
log_option(PAPILO AUTO "PAPILO_DIR")
log_option(READLINE AUTO "READLINE_DIR")
log_option(SANITIZE OFF "")
log_option(SHARED AUTO "")
log_option(SYM "SNAUTY" "BLISS_DIR;NAUTY_DIR")
log_option(THREADSAFE ON "")
log_option(TPI AUTO "")
log_option(WORHP AUTO "WORHP_DIR")
log_option(ZIMPL AUTO "ZIMPL_DIR")
log_option(ZLIB AUTO "ZLIB_DIR")
file(APPEND ${CMAKE_BINARY_DIR}/cmake.log "\n")

# If we include this directory then we mark SCIP as found.
set(SCIP_FOUND 1)

