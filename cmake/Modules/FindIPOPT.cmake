# SPDX-FileCopyrightText: 2012-2021 Istituto Italiano di Tecnologia (IIT)
# SPDX-FileCopyrightText: 2008-2010 RobotCub Consortium
# SPDX-License-Identifier: BSD-3-Clause

#
# Try to locate the IPOPT library
#
# On non Windows systems, use pkg-config to try to locate the library,
# if this fails then try to locate the library in the directory pointed by
# the IPOPT_DIR enviromental variable.
#
# On Windows systems,  just try to find the library using the IPOPT_DIR
# enviromental variable.
#
# Create the following variables::
#
#  IPOPT_INCLUDE_DIRS - Directories to include to use IPOPT
#  IPOPT_LIBRARIES    - Default library to link against to use IPOPT
#  IPOPT_FOUND        - If false, don't try to use IPOPT

#=============================================================================
# Copyright (C) 2008-2010 RobotCub Consortium
# Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
#   Authors: Ugo Pattacini <ugo.pattacini@iit.it>
#   Authors: Daniele E. Domenichelli <daniele.domenichelli@iit.it>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of YCM, substitute the full
#  License text for the above reference.)

if(NOT WIN32)
  # On non Windows systems we use PkgConfig to find IPOPT
  find_package(PkgConfig QUIET)

  if(DEFINED IPOPT_DIR)
      message("Searching for IPOPT in ${IPOPT_DIR}")
  else()
    set(IPOPT_DIR_TEST $ENV{IPOPT_DIR})
    if(DEFINED ENV{IPOPT_DIR})
      set(IPOPT_DIR $ENV{IPOPT_DIR} CACHE PATH "Path to IPOPT build directory")
    endif()
  endif()

  if(PKG_CONFIG_FOUND)
    if(IPOPT_DIR)
      set(ENV{PKG_CONFIG_PATH} "${IPOPT_DIR}/lib/pkgconfig/:$ENV{PKG_CONFIG_PATH}")
    endif()

    if(IPOPT_FIND_VERSION)
      if(IPOPT_FIND_VERSION_EXACT)
        pkg_check_modules(_PC_IPOPT QUIET IMPORTED_TARGET ipopt=${IPOPT_FIND_VERSION})
      else()
        pkg_check_modules(_PC_IPOPT QUIET IMPORTED_TARGET ipopt>=${IPOPT_FIND_VERSION})
      endif()
    else()
      pkg_check_modules(_PC_IPOPT QUIET IMPORTED_TARGET ipopt)
    endif()
  endif()

  # package config installations seem to work only with ipopt 3.13 or higher
  if(_PC_IPOPT_FOUND AND _PC_IPOPT_VERSION VERSION_GREATER_EQUAL 3.13)
    set(IPOPT_INCLUDE_DIRS ${_PC_IPOPT_INCLUDE_DIRS} CACHE PATH "IPOPT include directory")
    set(IPOPT_LIBRARIES PkgConfig::_PC_IPOPT CACHE STRING "IPOPT libraries" FORCE)
  else()
  # if no package config or it found a version too old (<= 3.12)
  # then try to find the package manually using IPOPT_DIR

    find_path(IPOPT_INCLUDE_DIRS
        NAMES IpoptConfig.h IpIpoptApplication.hpp
        HINTS ${IPOPT_DIR} /usr /opt/homebrew
        PATH_SUFFIXES include/coin include/coin-or
        )

    find_library(IPOPT_LIBRARIES
        NAMES ipopt
        HINTS ${IPOPT_DIR} /usr /opt/homebrew
        PATH_SUFFIXES lib
        )

    if(IPOPT_LIBRARIES)

      find_file(IPOPT_DEP_FILE
          NAMES ipopt_addlibs_cpp.txt
          HINTS ${IPOPT_DIR} /usr /opt/homebrew
          PATH_SUFFIXES share/coin/doc/Ipopt
          )
      mark_as_advanced(IPOPT_DEP_FILE)

      if(IPOPT_DEP_FILE)
        # add libraries from ipopt_addlibs_cpp.txt
        file(READ ${IPOPT_DEP_FILE} IPOPT_DEP)
        string(REGEX REPLACE "-[^l][^ ]* " "" IPOPT_DEP ${IPOPT_DEP})
        string(REPLACE "-l"                "" IPOPT_DEP ${IPOPT_DEP})
        string(REPLACE "\n"                "" IPOPT_DEP ${IPOPT_DEP})
        string(REPLACE "ipopt"             "" IPOPT_DEP ${IPOPT_DEP})       # remove any possible auto-dependency
        separate_arguments(IPOPT_DEP)

        # use the find_library command in order to prepare rpath correctly
        foreach(LIB ${IPOPT_DEP})
          find_library(IPOPT_SEARCH_FOR_${LIB} ${LIB} ${IPOPT_DIR}/lib)
          if(IPOPT_SEARCH_FOR_${LIB})
            # handle non-system libraries (e.g. coinblas)
            set(IPOPT_LIBRARIES ${IPOPT_LIBRARIES} ${IPOPT_SEARCH_FOR_${LIB}})
          else()
            # handle system libraries (e.g. gfortran)
            set(IPOPT_LIBRARIES ${IPOPT_LIBRARIES} ${LIB})
          endif()
          mark_as_advanced(IPOPT_SEARCH_FOR_${LIB})
        endforeach()
      endif()
    endif()

    set(IPOPT_DEFINITIONS "")
    set(IPOPT_LINK_FLAGS "")
  endif()

# Windows platforms
else()
  include(SelectLibraryConfigurations)

  set(IPOPT_DIR $ENV{IPOPT_DIR} CACHE PATH "Path to IPOPT build directory")

  find_path(IPOPT_INCLUDE_DIRS
      NAMES IpIpoptApplication.hpp
      PATH_SUFFIXES coin coin-or
      PATHS ${IPOPT_DIR}/include/coin)
  if(NOT IPOPT_INCLUDE_DIRS)
      # version ipopt >= 3.13
      set(IPOPT_INCLUDE_DIRS ${IPOPT_DIR}/include/coin-or)
      if(NOT EXISTS "${IPOPT_INCLUDE_DIRS}")
          # version ipopt <= 3.12
          set(IPOPT_INCLUDE_DIRS ${IPOPT_DIR}/include/coin)
      endif()
  endif()


  find_library(IPOPT_IPOPT_LIBRARY_RELEASE
      NAMES libipopt ipopt ipopt-3 ipopt-0 libipopt-3 libipopt-0
      HINTS ${IPOPT_DIR}/lib)

  select_library_configurations(IPOPT_IPOPT)
  set(IPOPT_LIBRARIES ${IPOPT_IPOPT_LIBRARY})

  find_file(IPOPT_DEP_FILE ipopt_addlibs_cpp.txt ${IPOPT_DIR}/share/doc/coin-or/Ipopt
                                                 ${IPOPT_DIR}/share/coin-or/doc/Ipopt
                                                 ${IPOPT_DIR}/share/doc/coin/Ipopt
                                                 ${IPOPT_DIR}/share/coin/doc/Ipopt
                                                 NO_DEFAULT_PATH)
  mark_as_advanced(IPOPT_DEP_FILE)

  if(IPOPT_DEP_FILE)
    # parse the file and acquire the dependencies
    file(READ ${IPOPT_DEP_FILE} IPOPT_DEP)

    string(REGEX REPLACE "-[^l][^ ]* " "" IPOPT_DEP ${IPOPT_DEP})
    string(REPLACE "\n"                "" IPOPT_DEP ${IPOPT_DEP})
    string(REPLACE "\n"                "" IPOPT_DEP ${IPOPT_DEP})
    string(REPLACE "ipopt"             "" IPOPT_DEP ${IPOPT_DEP})       # remove any possible auto-dependency
    separate_arguments(IPOPT_DEP)

    # use the find_library command in order to prepare rpath correctly
    foreach(LIB ${IPOPT_DEP})

      # skip LD library flags (this can be either -libflags or -l)
      if(${LIB} MATCHES "-l*")
        continue()
      endif()

      # check whether we compile for x86 or x64
      if(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
        set(MKL_ARCH_DIR "intel64")
      else()
        set(MKL_ARCH_DIR "ia32")
      endif()

      find_library(IPOPT_SEARCH_FOR_${LIB} ${LIB} $ENV{MKLROOT}/lib/${MKL_ARCH_DIR}
                                                  ${IPOPT_DIR}/lib
                                                  ${IPOPT_DIR}/lib/coin-or
                                                  ${IPOPT_DIR}/lib/coin-or/ThirdParty
                                                  ${IPOPT_DIR}/lib/coin
                                                  ${IPOPT_DIR}/lib/coin/ThirdParty
                                                  NO_DEFAULT_PATH)

      if(IPOPT_SEARCH_FOR_${LIB})
        set(IPOPT_LIBRARIES ${IPOPT_LIBRARIES} ${IPOPT_SEARCH_FOR_${LIB}})
      endif()
      mark_as_advanced(IPOPT_SEARCH_FOR_${LIB})
    endforeach()

  else()
    # if no dep file was found, i.e. newer versions >= 3.13 of ipopt

    # Some very old version of binary releases of IPOPT have Intel fortran
    # libraries embedded in the library, some less older releases require them to
    # be explicitly linked. Current releases include the Intel fortran libraries
    # as dependencies of the DLL.
    if(IPOPT_IPOPT_LIBRARY)
      get_filename_component(_MSVC_DIR "${CMAKE_LINKER}" DIRECTORY)

      # Find the lib.exe executable
      find_program(LIB_EXECUTABLE
                   NAMES lib.exe
                   HINTS "${_MSVC_BINDIR}"
                         "C:/Program Files/Microsoft Visual Studio 10.0/VC/bin"
                         "C:/Program Files (x86)/Microsoft Visual Studio 10.0/VC/bin"
                         "C:/Program Files/Microsoft Visual Studio 11.0/VC/bin"
                         "C:/Program Files (x86)/Microsoft Visual Studio 11.0/VC/bin"
                         "C:/Program Files/Microsoft Visual Studio 12.0/VC/bin"
                         "C:/Program Files (x86)/Microsoft Visual Studio 12.0/VC/bin"
                         "C:/Program Files/Microsoft Visual Studio 14.0/VC/bin"
                         "C:/Program Files (x86)/Microsoft Visual Studio 14.0/VC/bin"
                   DOC "Path to the lib.exe executable")
      mark_as_advanced(LIB_EXECUTABLE)

      # backup PATH environment variable
      set(_path $ENV{PATH})

      # Add th MSVC "Common7/IDE" dir containing the dlls in the PATH when needed.
      get_filename_component(_MSVC_LIBDIR "${_MSVC_BINDIR}/../../Common7/IDE" ABSOLUTE)
      if(NOT EXISTS "${_MSVC_LIBDIR}")
        get_filename_component(_MSVC_LIBDIR "${_MSVC_BINDIR}/../../../Common7/IDE" ABSOLUTE)
      endif()

      if(EXISTS "${_MSVC_LIBDIR}")
        set(_MSVC_LIBDIR_FOUND 0)
        file(TO_CMAKE_PATH "$ENV{PATH}" _env_path)
        foreach(_dir ${_env_path})
          if("${_dir}" STREQUAL ${_MSVC_LIBDIR})
            set(_MSVC_LIBDIR_FOUND 1)
          endif()
        endforeach()
        if(NOT _MSVC_LIBDIR_FOUND)
          file(TO_NATIVE_PATH "${_MSVC_LIBDIR}" _MSVC_LIBDIR)
          set(ENV{PATH} "$ENV{PATH};${_MSVC_LIBDIR}")
        endif()
      endif()

      set(_IPOPT_LIB ${IPOPT_IPOPT_LIBRARY_RELEASE})

      execute_process(COMMAND ${LIB_EXECUTABLE} /list "${_IPOPT_LIB}"
                      OUTPUT_VARIABLE _lib_output)

      set(ENV{PATH} "${_path}")
      unset(_path)

      if(NOT "${_lib_output}" MATCHES "libifcoremd.dll")
        get_filename_component(_IPOPT_IPOPT_LIBRARY_DIR "${_IPOPT_LIB}" DIRECTORY)

        foreach(_lib ifconsol
                     libifcoremd
                     libifportmd
                     libmmd
                     libirc
                     svml_dispmd)
          string(TOUPPER "${_lib}" _LIB)
          find_library(IPOPT_${_LIB}_LIBRARY_RELEASE ${_lib} ${_IPOPT_IPOPT_LIBRARY_DIR})
          select_library_configurations(IPOPT_${_LIB})
          if(NOT "${IPOPT_${_LIB}_LIBRARY}" MATCHES "NOTFOUND$")
            list(APPEND IPOPT_LIBRARIES ${IPOPT_${_LIB}_LIBRARY})
          endif()
        endforeach()
      endif()
    endif()

    set(IPOPT_DEFINITIONS "")
    if(MSVC)
      set(IPOPT_LINK_FLAGS "/NODEFAULTLIB:libcmt.lib;libcmtd.lib")
    else()
      set(IPOPT_LINK_FLAGS "")
    endif()
  endif()

endif()

# parse the version number
foreach( INCLUDE_DIR ${IPOPT_INCLUDE_DIRS} )
  if( EXISTS ${INCLUDE_DIR}/IpoptConfig.h )
  file(STRINGS ${INCLUDE_DIR}/IpoptConfig.h CONFIGFILE)
    foreach(STR ${CONFIGFILE})
      if("${STR}" MATCHES "^#define IPOPT_VERSION ")
        string(REGEX REPLACE "#define IPOPT_VERSION " "" IPOPT_VERSION ${STR})
        string(REGEX REPLACE "\"" "" IPOPT_VERSION ${IPOPT_VERSION})
      endif()
    endforeach()
  endif()
endforeach()

mark_as_advanced(IPOPT_INCLUDE_DIRS
                 IPOPT_LIBRARIES
                 IPOPT_LINK_FLAGS)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IPOPT
  FOUND_VAR IPOPT_FOUND
  REQUIRED_VARS IPOPT_LIBRARIES
  VERSION_VAR IPOPT_VERSION)
