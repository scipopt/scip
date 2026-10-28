# Finds the MOSEK C library.
#
# == Imported targets ==
#
# This module defines the following IMPORTED targets:
# - MOSEK::MOSEK (interfaces MOSEK::libmosek)
# - MOSEK::libmosek (C library)
#
# == Result variables ==
#
# This module will set the following variables in your project:
# - MOSEK_FOUND (bool indicating presence of MOSEK::libmosek)
# - MOSEK_VERSION (version string)
#
# == Hints ==
#
# MOSEK_DIR (CMake or environment variable) is inspected when searching for the library.

include(FindPackageHandleStandardArgs)
if(NOT COMMAND find_package_check_version)
   function(find_package_check_version version output_var HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE result_var)
      if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION)
         message(WARNING "Accepting version ${version} as compatible with ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION}; CMake version >= 3.19 needed to actually check versions.")
      endif()
      set(${output_var} 1 PARENT_SCOPE)
   endfunction()
endif()

# Parse MOSEK version from an include directory.
function(mosek_parse_version result include_dir)
   file(READ "${include_dir}/mosek.h" header_file)
   string(REGEX MATCH "#define[ \t]+MSK_VERSION_MAJOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_major ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+MSK_VERSION_MINOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_minor ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+MSK_VERSION_REVISION[ \t]+([0-9]+)" _ ${header_file})
   set(version_patch ${CMAKE_MATCH_1})
   if(version_major)
      set(${result} "${version_major}.${version_minor}.${version_patch}" PARENT_SCOPE)
      return()
   else()
      set(${result} NONE PARENT_SCOPE)
      return()
   endif()
endfunction()

# Search for the C library.
function(mosek_find_main_lib result search_dir version)
   message(DEBUG "    mosek_find_main_lib(\"${search_dir}\", \"${version_list}\") with mosek_no_default_path=<${mosek_no_default_path}>")

   # Extract major and minor version.
   string(REPLACE "." ";" version_list "${version}")
   list(GET version_list 0 version_major)
   list(GET version_list 1 version_minor)

   if(MSVC)
      set(lib_name "mosek64_${version_major}_${version_minor}")
   elseif(APPLE)
      set(lib_name "mosek64.${version_major}.${version_minor}.dylib")
   else()
      set(lib_name "mosek64")
   endif()

   find_library(mosek_main_lib
      NAMES
         "${lib_name}"
      PATHS
         "${search_dir}/../bin"
      ${mosek_no_default_path}
      )
   message(DEBUG "    -> ${mosek_main_lib}, stored in ${result}")
   set(${result} "${mosek_main_lib}" PARENT_SCOPE)
endfunction()

# Inspection of a candidate include directory.
function(mosek_validate result include_dir)
   message(DEBUG "  mosek_validate(\"${include_dir}\") with mosek_no_default_path=<${mosek_no_default_path}>")

   # Parse version.
   mosek_parse_version(version "${include_dir}")
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
   mosek_find_main_lib(mosek_main_lib "${include_dir}" "${version}")
   message(DEBUG "  -> mosek_main_lib=${mosek_main_lib}")
   if(NOT mosek_main_lib)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without a library.")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   set(${result} TRUE PARENT_SCOPE)
   return()
endfunction()

# If MOSEK was already found before, check version.
if(${CMAKE_FIND_PACKAGE_NAME}_FOUND
      AND ${CMAKE_FIND_PACKAGE_NAME}_VERSION
      AND TARGET MOSEK::MOSEK
      AND TARGET MOSEK::libmosek
   )
   find_package_check_version(${${CMAKE_FIND_PACKAGE_NAME}_VERSION} mosek_version_result HANDLE_VERSION_RANGE)
   if(mosek_version_result)
      if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
         message(STATUS "Already found suitable ${CMAKE_FIND_PACKAGE_NAME} version \"${${CMAKE_FIND_PACKAGE_NAME}_VERSION}\".")
      endif()
      unset(mosek_version_result)
      return()
   endif()
   unset(mosek_version_result)
endif()

# Main search function whose variables are local.
function(mosek_find)
   if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Searching for ${CMAKE_FIND_PACKAGE_NAME}.")
   endif()

   # Determine path suffix depending on platform.
   if(MSVC)
      set(mosek_platform "win64x86")
   elseif(APPLE)
      set(mosek_platform "osxaarch64")
   elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
      set(mosek_platform "linux64x86")
   else()
      set(mosek_platform "linuxaarch64")
   endif()


   foreach(find_local TRUE;FALSE)
      # Whether to find libraries only at the location derived from the header.
      if(find_local)
         set(mosek_no_default_path NO_DEFAULT_PATH)
      else()
         set(mosek_no_default_path)
      endif()

      find_path(mosek_main_inc
         NAMES
            "mosek.h"
         HINTS
            "${MOSEK_DIR}"
            "$ENV{MOSEK_DIR}"
         PATH_SUFFIXES
            "h"
            "${mosek_platform}/h"
            "platform/${mosek_platform}/h"
            "tools/platform/${mosek_platform}/h/"
         VALIDATOR
            mosek_validate
         )

      message(DEBUG "With mosek_no_default_path=<${mosek_no_default_path}>, find_path returned mosek_main_inc <${mosek_main_inc}>.")
      if(mosek_main_inc)
         # The VALIDATOR feature does not exist before 3.25, so we just manually call the validate function.
         if(CMAKE_VERSION VERSION_LESS 3.25)
            message(DEBUG "A-posteriori validating candidate include directory <${mosek_main_inc}>.")
            mosek_validate(mosek_validation_result "${mosek_main_inc}")
            if(mosek_validation_result)
               break()
            else()
               unset(mosek_main_inc CACHE)
            endif()
         else()
            break()
         endif()
      endif()
   endforeach()

   if(mosek_main_inc)
      mosek_parse_version(mosek_version "${mosek_main_inc}")
      set(${CMAKE_FIND_PACKAGE_NAME}_VERSION "${mosek_version}" PARENT_SCOPE)
      mosek_find_main_lib(mosek_main_lib "${mosek_main_inc}" "${mosek_version}")
      set(mosek_list_libs "${mosek_main_lib}")
      if(mosek_main_lib)

         # C library target
         if(NOT TARGET MOSEK::libmosek)
            add_library(MOSEK::libmosek UNKNOWN IMPORTED)
            set_target_properties(MOSEK::libmosek PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${mosek_main_inc}"
               IMPORTED_LOCATION "${mosek_main_lib}"
               )
         endif()

         # Main interface target
         if(NOT TARGET MOSEK::MOSEK)
            # Note that we do not create an ALIAS to support CMake version < 3.15;
            # see https://cmake.org/cmake/help/latest/release/3.15.html#commands
            # Moreover, this allows a consuming package to ALIAS it.
            add_library(MOSEK::MOSEK INTERFACE IMPORTED)
            target_link_libraries(MOSEK::MOSEK INTERFACE MOSEK::libmosek)
         endif()

         # Deprecated variables refer to a deprecated interface target.
         if(NOT TARGET MOSEK_LIBRARIES)
            add_library(MOSEK_LIBRARIES INTERFACE IMPORTED)
            if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
               set_target_properties(MOSEK_LIBRARIES PROPERTIES
                  DEPRECATION "Do not use ${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS and ${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES anymore, but the MOSEK::MOSEK target instead.")
            else()
               message(DEBUG "CMake version ${CMAKE_VERSION} is too low to mark MOSEK_LIBRARIES target as deprecated!")
            endif()
            target_link_libraries(MOSEK_LIBRARIES INTERFACE MOSEK::libmosek)
            set(${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES MOSEK_LIBRARIES PARENT_SCOPE)
            get_target_property(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS MOSEK::libmosek INTERFACE_INCLUDE_DIRECTORIES)
            set(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS "${${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS}" PARENT_SCOPE)
         endif()

         if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME} version ${mosek_version}: ${mosek_main_inc} and ${mosek_list_libs}")
         endif()
         set(${CMAKE_FIND_PACKAGE_NAME}_FOUND TRUE PARENT_SCOPE)
         return()
      else()
         set(mosek_fail_reason ": unexpected failure (mosek_main_lib=${mosek_main_lib})")
      endif()
   else()
      set(mosek_fail_reason "")
   endif()

   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE)
      set(mosek_fail_version " version ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE}")
   else()
      set(mosek_fail_version "")
   endif()
   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find required dependency MOSEK${mosek_fail_version}${mosek_fail_reason}. Try specifying suitable MOSEK_DIR!")
   elseif(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(DEBUG "Could NOT find MOSEK${mosek_fail_version}${mosek_fail_reason}.")
   else()
      message(STATUS "Could NOT find optional dependency MOSEK${mosek_fail_version}${mosek_fail_reason}. Try specifying suitable MOSEK_DIR.")
   endif()
   set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE PARENT_SCOPE)

endfunction()

# Invoke main search.
mosek_find()

