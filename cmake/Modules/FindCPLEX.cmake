# Finds the CPLEX C library.
#
# == Imported targets ==
#
# This module defines the following IMPORTED targets:
# - CPLEX::CPLEX (interfaces CPLEX::libcplex)
# - CPLEX::libcplex (C library)
#
# == Result variables ==
#
# This module will set the following variables in your project:
# - CPLEX_FOUND (bool indicating presence of CPLEX::libcplex)
# - CPLEX_VERSION (version string)
#
# == Hints ==
#
# CPLEX_DIR (CMake or environment variable) is inspected when searching for the library.

include(FindPackageHandleStandardArgs)
if(NOT COMMAND find_package_check_version)
   function(find_package_check_version version output_var HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE result_var)
      if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION)
         message(WARNING "Accepting version ${version} as compatible with ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION}; CMake version >= 3.19 needed to actually check versions.")
      endif()
      set(${output_var} 1 PARENT_SCOPE)
   endfunction()
endif()

# Parse CPLEX version from an include directory.
function(cplex_parse_version result include_dir)
   file(READ "${include_dir}/cpxconst.h" header_file)
   string(REGEX MATCH "#define CPX_VERSION *([0-9]*)" _ ${header_file})
   set(version_int ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+CPX_VERSION_VERSION[ \t]+([0-9]+)" _ ${header_file})
   set(version_major ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+CPX_VERSION_RELEASE[ \t]+([0-9]+)" _ ${header_file})
   set(version_minor ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+CPX_VERSION_MODIFICATION[ \t]+([0-9]+)" _ ${header_file})
   set(version_patch ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+CPX_VERSION_FIX[ \t]+([0-9]+)" _ ${header_file})
   set(version_tweak ${CMAKE_MATCH_1})
   if(version_major)
      set(${result} "${version_int};${version_major}.${version_minor}.${version_patch}.${version_tweak}" PARENT_SCOPE)
      return()
   else()
      set(${result} NONE PARENT_SCOPE)
      return()
   endif()
endfunction()

# Search for the C library.
function(cplex_find_main_lib result search_dir version_list)
   message(DEBUG "    cplex_find_main_lib(\"${search_dir}\", \"${version_list}\") with cplex_no_default_path=<${cplex_no_default_path}>")

   if(MSVC)
      list(GET version_list 1 version_int)
      string(SUBSTRING "${version_int}" 0 4 version_int)
      set(lib_name "cplex${version_int}")

      string(REGEX REPLACE "/VC/bin/.*" "" VISUAL_STUDIO_PATH ${CMAKE_CXX_COMPILER})
      string(REGEX MATCH "Studio/[0-9]+/" CPLEX_WIN_VS_VERSION ${VISUAL_STUDIO_PATH})
      string(REGEX REPLACE "Studio/" "" CPLEX_WIN_VS_VERSION ${CPLEX_WIN_VS_VERSION})
      string(REGEX REPLACE "/" "" CPLEX_WIN_VS_VERSION ${CPLEX_WIN_VS_VERSION})
      set(lib_dir "${search_dir}/../../lib/x64_windows_vs${CPLEX_WIN_VS_VERSION}/stat_m")
      if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
         set(lib_dir "${lib_dir}/stat_mt") # stat_mt{d,a}
      else()
         set(lib_dir "${lib_dir}/stat_md") # stat_md{d,a}
      endif()
      if(MT)
         set(lib_dir "${lib_dir}d") # stat_m{t,d}d
      else()
         set(lib_dir "${lib_dir}a") # stat_m{t,d}a
      endif()
   elseif(APPLE)
      set(lib_name "cplex")
      string(REGEX REPLACE "x86_64" "x86-64" lib_dir "${search_dir}/../../lib/${CMAKE_OSX_ARCHITECTURES}_osx/static_pic")
   else()
      set(lib_name "cplex")
      string(REGEX REPLACE "x86_64" "x86-64" lib_dir "${search_dir}/../../lib/${CMAKE_SYSTEM_PROCESSOR}_linux/static_pic")
   endif()

   find_library(cplex_main_lib
      NAMES
         "${lib_name}"
      PATHS
         "${lib_dir}"
      ${cplex_no_default_path}
      )
   message(DEBUG "    -> ${cplex_main_lib}, stored in ${result}")
   set(${result} "${cplex_main_lib}" PARENT_SCOPE)
endfunction()

# Inspection of a candidate include directory.
function(cplex_validate result include_dir)
   message(DEBUG "  cplex_validate(\"${include_dir}\") with cplex_no_default_path=<${cplex_no_default_path}>")

   # Parse version.
   cplex_parse_version(version_list "${include_dir}")
   if(version_list)
      list(GET version_list 1 version)
      find_package_check_version("${version}" version_okay HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE reason)
   else()
      set(version_okay FALSE)
      set(reason "unable to determine version.")
   endif()
   if(NOT version_okay)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} ignored: ${reason}")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   # Check for library.
   cplex_find_main_lib(cplex_main_lib "${include_dir}" "${version_list}")
   message(DEBUG "  -> cplex_main_lib=${cplex_main_lib}")
   if(NOT cplex_main_lib)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without a library.")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   set(${result} TRUE PARENT_SCOPE)
   return()
endfunction()

# If CPLEX was already found before, check version.
if(${CMAKE_FIND_PACKAGE_NAME}_FOUND
      AND ${CMAKE_FIND_PACKAGE_NAME}_VERSION
      AND TARGET CPLEX::CPLEX
      AND TARGET CPLEX::libcplex
   )
   find_package_check_version(${${CMAKE_FIND_PACKAGE_NAME}_VERSION} cplex_version_result HANDLE_VERSION_RANGE)
   if(cplex_version_result)
      if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
         message(STATUS "Already found suitable ${CMAKE_FIND_PACKAGE_NAME} version ${${CMAKE_FIND_PACKAGE_NAME}_VERSION}.")
      endif()
      unset(cplex_version_result)
      return()
   endif()
   unset(cplex_version_result)
endif()

# Main search function whose variables are local.
function(cplex_find)
   if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Searching for ${CMAKE_FIND_PACKAGE_NAME}.")
   endif()

   foreach(find_local TRUE;FALSE)
      # Whether to find libraries only at the location derived from the header.
      if(find_local)
         set(cplex_no_default_path NO_DEFAULT_PATH)
      else()
         set(cplex_no_default_path)
      endif()

      find_path(cplex_main_inc
         NAMES
            "cpxconst.h"
         HINTS
            "${CPLEX_DIR}"
            "$ENV{CPLEX_DIR}"
         PATH_SUFFIXES
            "include"
            "include/ilcplex"
            "cplex/include/ilcplex"
         VALIDATOR
            cplex_validate
         )

      message(DEBUG "With cplex_no_default_path=<${cplex_no_default_path}>, find_path returned cplex_main_inc <${cplex_main_inc}>.")
      if(cplex_main_inc)
         # The VALIDATOR feature does not exist before 3.25, so we just manually call the validate function.
         if(CMAKE_VERSION VERSION_LESS 3.25)
            message(DEBUG "A-posteriori validating candidate include directory <${cplex_main_inc}>.")
            cplex_validate(cplex_validation_result "${cplex_main_inc}")
            if(cplex_validation_result)
               break()
            else()
               unset(cplex_main_inc CACHE)
            endif()
         else()
            break()
         endif()
      endif()
   endforeach()

   if(cplex_main_inc)
      cplex_parse_version(cplex_version_list "${cplex_main_inc}")
      list(GET cplex_version_list 1 cplex_version)
      set(${CMAKE_FIND_PACKAGE_NAME}_VERSION "${cplex_version}" PARENT_SCOPE)
      cplex_find_main_lib(cplex_main_lib "${cplex_main_inc}" "${cplex_version_list}")
      set(cplex_list_libs "${cplex_main_lib}")
      if(cplex_main_lib)

         # C library target
         if(NOT TARGET CPLEX::libcplex)
            add_library(CPLEX::libcplex UNKNOWN IMPORTED)
            set_target_properties(CPLEX::libcplex PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${cplex_main_inc}"
               IMPORTED_LOCATION "${cplex_main_lib}"
               )
         endif()

         # Main interface target
         if(NOT TARGET CPLEX::CPLEX)
            # Note that we do not create an ALIAS to support CMake version < 3.15;
            # see https://cmake.org/cmake/help/latest/release/3.15.html#commands
            # Moreover, this allows a consuming package to ALIAS it.
            add_library(CPLEX::CPLEX INTERFACE IMPORTED)
            target_link_libraries(CPLEX::CPLEX INTERFACE CPLEX::libcplex)
         endif()

         # Deprecated variables refer to a deprecated interface target.
         if(NOT TARGET CPLEX_LIBRARIES)
            add_library(CPLEX_LIBRARIES INTERFACE IMPORTED)
            if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
               set_target_properties(CPLEX_LIBRARIES PROPERTIES
                  DEPRECATION "Do not use ${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS and ${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES anymore, but the CPLEX::CPLEX target instead.")
            else()
               message(DEBUG "CMake version ${CMAKE_VERSION} is too low to mark CPLEX_LIBRARIES target as deprecated!")
            endif()
            target_link_libraries(CPLEX_LIBRARIES INTERFACE CPLEX::libcplex)
            set(${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES CPLEX_LIBRARIES PARENT_SCOPE)
            get_target_property(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS CPLEX::libcplex INTERFACE_INCLUDE_DIRECTORIES)
            set(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS "${${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS}" PARENT_SCOPE)
         endif()

         if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME} version ${cplex_version}: ${cplex_main_inc} and ${cplex_list_libs}")
         endif()
         set(${CMAKE_FIND_PACKAGE_NAME}_FOUND TRUE PARENT_SCOPE)
         return()
      else()
         set(cplex_fail_reason ": unexpected failure (cplex_main_lib=${cplex_main_lib})")
      endif()
   else()
      set(cplex_fail_reason "")
   endif()

   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE)
      set(cplex_fail_version " version ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE}")
   else()
      set(cplex_fail_version "")
   endif()
   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find required dependency CPLEX${cplex_fail_version}${cplex_fail_reason}. Try specifying suitable CPLEX_DIR!")
   elseif(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(DEBUG "Could NOT find CPLEX${cplex_fail_version}${cplex_fail_reason}.")
   else()
      message(STATUS "Could NOT find optional dependency CPLEX${cplex_fail_version}${cplex_fail_reason}. Try specifying suitable CPLEX_DIR.")
   endif()
   set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE PARENT_SCOPE)

endfunction()

# Invoke main search.
cplex_find()

