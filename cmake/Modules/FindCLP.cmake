# Finds the CLP C++ library.
#
# == Imported targets ==
#
# This module defines the following IMPORTED targets:
# - CLP::CLP (interfaces CLP::libclp)
# - CLP::libclp (C++ library)
#
# == Result variables ==
#
# This module will set the following variables in your project:
# - CLP_FOUND (bool indicating presence of CLP::libclp)
# - CLP_VERSION (version string)
#
# == Hints ==
#
# CLP_DIR (CMake or environment variable) is inspected when searching for the library.

# For `find_dependency` function.
include(CMakeFindDependencyMacro)

include(FindPackageHandleStandardArgs)
if(NOT COMMAND find_package_check_version)
   function(find_package_check_version version output_var HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE result_var)
      if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION)
         message(WARNING "Accepting version ${version} as compatible with ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION}; CMake version >= 3.19 needed to actually check versions.")
      endif()
      set(${output_var} 1 PARENT_SCOPE)
   endfunction()
endif()

# Parse CLP version from an include directory.
function(clp_parse_version result include_dir)
   file(READ "${include_dir}/ClpConfig.h" header_file)
   string(REGEX MATCH "#define[ \t]+CLP_VERSION_MAJOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_major ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+CLP_VERSION_MINOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_minor ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+CLP_VERSION_RELEASE[ \t]+([0-9]+)" _ ${header_file})
   set(version_patch ${CMAKE_MATCH_1})
   if(version_major)
      set(${result} "${version_major}.${version_minor}.${version_patch}" PARENT_SCOPE)
      return()
   else()
      set(${result} NONE PARENT_SCOPE)
      return()
   endif()
endfunction()

# Search for the C++ library.
function(clp_find_main_lib result search_dir)
   message(DEBUG "    clp_find_main_lib(\"${search_dir}\") with clp_no_default_path=<${clp_no_default_path}>")
   find_library(clp_main_lib
      NAMES
         "Clp"
      PATHS
         "${search_dir}/../../lib/"
      ${clp_no_default_path}
      )
   message(DEBUG "    -> ${clp_main_lib}, stored in ${result}")
   set(${result} "${clp_main_lib}" PARENT_SCOPE)
endfunction()

# Inspection of a candidate include directory.
function(clp_validate result include_dir)
   message(DEBUG "  clp_validate(\"${include_dir}\") with clp_no_default_path=<${clp_no_default_path}>")

   # Parse version.
   clp_parse_version(version "${include_dir}")
   if(version)
      find_package_check_version(${version} version_okay HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE reason)
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
   clp_find_main_lib(clp_main_lib "${include_dir}")
   message(DEBUG "  -> clp_main_lib=${clp_main_lib}")
   if(NOT clp_main_lib)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without a library.")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   message(DEBUG "    Accepting ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} with library ${clp_main_lib} and version ${version}.")
   set(${result} TRUE PARENT_SCOPE)
   return()
endfunction()

# If CLP was already found before, check version.
if(${CMAKE_FIND_PACKAGE_NAME}_FOUND
      AND ${CMAKE_FIND_PACKAGE_NAME}_VERSION
      AND TARGET CLP::CLP
      AND TARGET CLP::libclp
   )
   find_package_check_version(${${CMAKE_FIND_PACKAGE_NAME}_VERSION} clp_version_result HANDLE_VERSION_RANGE)
   if(clp_version_result)
      if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
         message(STATUS "Already found suitable ${CMAKE_FIND_PACKAGE_NAME} version ${${CMAKE_FIND_PACKAGE_NAME}_VERSION}.")
      endif()
      unset(clp_version_result)
      return()
   endif()
   unset(clp_version_result)
endif()

# Main search function whose variables are local.
function(clp_find)
   if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Searching for ${CMAKE_FIND_PACKAGE_NAME}.")
   endif()

   # Indent messages for dependency search.
   list(APPEND CMAKE_MESSAGE_INDENT "  ")

   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      find_package(CoinUtils QUIET)
   else()
      find_package(CoinUtils)
   endif()

   # Unindent messages again; equivalent to list(POP_BACK, CMAKE_MESSAGE_INDENT) which requires CMake 3.15.
   list(REVERSE CMAKE_MESSAGE_INDENT)
   list(REMOVE_AT CMAKE_MESSAGE_INDENT 0)
   list(REVERSE CMAKE_MESSAGE_INDENT)
   
   # Stop if some dependency is missing.
   if(NOT CoinUtils_FOUND)
      if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
         message(FATAL_ERROR "Could NOT find ${CMAKE_FIND_PACKAGE_NAME} due to missing dependency.")
      elseif(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
         message(DEBUG "Could NOT find ${CMAKE_FIND_PACKAGE_NAME} due to missing dependency.")
      else()
         message(STATUS "Could NOT find ${CMAKE_FIND_PACKAGE_NAME} due to missing dependency.")
      endif()
      set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE PARENT_SCOPE)
      return()
   endif()

   foreach(find_local TRUE;FALSE)
      # Whether to find libraries only at the location derived from the header.
      if(find_local)
         set(clp_no_default_path NO_DEFAULT_PATH)
      else()
         set(clp_no_default_path)
      endif()

      find_path(clp_main_inc
         NAMES
            "ClpConfig.h"
         HINTS
            "${CLP_DIR}"
            "$ENV{CLP_DIR}"
            "${COIN_DIR}"
            "$ENV{COIN_DIR}"
         PATH_SUFFIXES
            "coin-or"
            "include/coin-or"
            "dist/include/coin-or"
         VALIDATOR
            clp_validate
         )

      message(DEBUG "With clp_no_default_path=<${clp_no_default_path}>, find_path returned clp_main_inc <${clp_main_inc}>.")
      if(clp_main_inc)
         # The VALIDATOR feature does not exist before 3.25, so we just manually call the validate function.
         if(CMAKE_VERSION VERSION_LESS 3.25)
            message(DEBUG "A-posteriori validating candidate include directory <${clp_main_inc}>.")
            clp_validate(clp_validation_result "${clp_main_inc}")
            if(clp_validation_result)
               break()
            else()
               unset(clp_main_inc CACHE)
            endif()
         else()
            break()
         endif()
      endif()
   endforeach()

   if(clp_main_inc)
      clp_parse_version(clp_version "${clp_main_inc}")
      set(${CMAKE_FIND_PACKAGE_NAME}_VERSION "${clp_version}" PARENT_SCOPE)
      clp_find_main_lib(clp_main_lib "${clp_main_inc}")
      set(clp_list_libs "${clp_main_lib}")
      if(clp_main_lib)

         # C++ library target
         if(NOT TARGET CLP::libclp)
            add_library(CLP::libclp UNKNOWN IMPORTED)
            set_target_properties(CLP::libclp PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${clp_main_inc}"
               IMPORTED_LOCATION "${clp_main_lib}"
               INTERFACE_LINK_LIBRARIES CoinUtils::CoinUtils
               )
         endif()

         # Main interface target
         if(NOT TARGET CLP::CLP)
            # Note that we do not create an ALIAS to support CMake version < 3.15;
            # see https://cmake.org/cmake/help/latest/release/3.15.html#commands
            # Moreover, this allows a consuming package to ALIAS it.
            add_library(CLP::CLP INTERFACE IMPORTED)
            target_link_libraries(CLP::CLP INTERFACE CLP::libclp)
         endif()

         # Deprecated variables refer to a deprecated interface target.
         if(NOT TARGET CLP_LIBRARIES)
            add_library(CLP_LIBRARIES INTERFACE IMPORTED)
            if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
               set_target_properties(CLP_LIBRARIES PROPERTIES
                  DEPRECATION "Do not use ${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS and ${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES anymore, but the CLP::CLP target instead.")
            else()
               message(DEBUG "CMake version ${CMAKE_VERSION} is too low to mark CLP_LIBRARIES target as deprecated!")
            endif()
            target_link_libraries(CLP_LIBRARIES INTERFACE CLP::libclp)
            set(${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES CLP_LIBRARIES PARENT_SCOPE)
            get_target_property(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS CLP::libclp INTERFACE_INCLUDE_DIRECTORIES)
            set(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS "${${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS}" PARENT_SCOPE)
         endif()

         if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME} version ${clp_version}: ${clp_main_inc} and ${clp_list_libs}")
         endif()
         set(${CMAKE_FIND_PACKAGE_NAME}_FOUND TRUE PARENT_SCOPE)
         return()
      else()
         set(clp_fail_reason ": unexpected failure (clp_main_lib=${clp_main_lib})")
      endif()
   else()
      set(clp_fail_reason "")
   endif()

   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE)
      set(clp_fail_version " version ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE}")
   else()
      set(clp_fail_version "")
   endif()
   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find required dependency CLP${clp_fail_version}${clp_fail_reason}. Try specifying suitable CLP_DIR!")
   elseif(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(DEBUG "Could NOT find CLP${clp_fail_version}${clp_fail_reason}.")
   else()
      message(STATUS "Could NOT find optional dependency CLP${clp_fail_version}${clp_fail_reason}. Try specifying suitable CLP_DIR.")
   endif()
   set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE PARENT_SCOPE)

endfunction()

# Invoke main search.
clp_find()

