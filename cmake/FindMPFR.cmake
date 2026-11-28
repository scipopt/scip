# Finds the MPFR C library.
#
# == Imported targets ==
#
# This module defines the following IMPORTED targets:
# - MPFR::MPFR (interfaces MPFR::libmpfr)
# - MPFR::libmpfr (C library)
#
# == Result variables ==
#
# This module will set the following variables in your project:
# - MPFR_FOUND (bool indicating presence of MPFR::libmpfr)
# - MPFR_VERSION (version string)
#
# == Hints ==
#
# MPFR_DIR (CMake or environment variable) is inspected when searching for the library.

include(FindPackageHandleStandardArgs)
if(NOT COMMAND find_package_check_version)
   function(find_package_check_version version output_var HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE result_var)
      if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION)
         message(WARNING "Accepting version ${version} as compatible with ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION}; CMake version >= 3.19 needed to actually check versions.")
      endif()
      set(${output_var} 1 PARENT_SCOPE)
   endfunction()
endif()

# Parse MPFR version from an include directory.
function(mpfr_parse_version result include_dir)
   file(READ "${include_dir}/mpfr.h" header_file)
   string(REGEX MATCH "#define[ \t]+MPFR_VERSION_MAJOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_major ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+MPFR_VERSION_MINOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_minor ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+MPFR_VERSION_PATCHLEVEL[ \t]+([0-9]+)" _ ${header_file})
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
function(mpfr_find_main_lib result search_dir)
   message(DEBUG "    mpfr_find_main_lib(\"${search_dir}\") with mpfr_no_default_path=<${mpfr_no_default_path}>")
   find_library(mpfr_main_lib
      NAMES
         "mpfr"
      PATHS
         "${search_dir}"
         "${search_dir}/../lib"
         "${search_dir}/../lib64"
         "${search_dir}/../lib/${CMAKE_SYSTEM_PROCESSOR}-linux-gnu"
      ${mpfr_no_default_path}
      )
   message(DEBUG "    -> ${mpfr_main_lib}, stored in ${result}")
   set(${result} "${mpfr_main_lib}" PARENT_SCOPE)
endfunction()

# Inspection of a candidate include directory.
function(mpfr_validate result include_dir)
   message(DEBUG "  mpfr_validate(\"${include_dir}\") with mpfr_no_default_path=<${mpfr_no_default_path}>")

   # Parse version.
   mpfr_parse_version(version "${include_dir}")
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
   mpfr_find_main_lib(mpfr_main_lib "${include_dir}")
   message(DEBUG "  -> mpfr_main_lib=${mpfr_main_lib}")
   if(NOT mpfr_main_lib)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without a library.")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   message(DEBUG "    Accepting ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} with library ${mpfr_main_lib} and version ${version}.")
   set(${result} TRUE PARENT_SCOPE)
   return()
endfunction()

# If MPFR was already found before, check version.
if(${CMAKE_FIND_PACKAGE_NAME}_FOUND
      AND ${CMAKE_FIND_PACKAGE_NAME}_VERSION
      AND TARGET MPFR::MPFR
      AND TARGET MPFR::libmpfr
   )
   find_package_check_version(${${CMAKE_FIND_PACKAGE_NAME}_VERSION} mpfr_version_result HANDLE_VERSION_RANGE)
   if(mpfr_version_result)
      if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
         message(STATUS "Already found suitable ${CMAKE_FIND_PACKAGE_NAME} version ${${CMAKE_FIND_PACKAGE_NAME}_VERSION}.")
      endif()
      unset(mpfr_version_result)
      return()
   endif()
   unset(mpfr_version_result)
endif()

# Main search function whose variables are local.
function(mpfr_find)
   if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Searching for ${CMAKE_FIND_PACKAGE_NAME}.")
   endif()

   foreach(find_local TRUE;FALSE)
      # Whether to find libraries only at the location derived from the header.
      if(find_local)
         set(mpfr_no_default_path NO_DEFAULT_PATH)
      else()
         set(mpfr_no_default_path)
      endif()

      find_path(mpfr_main_inc
         NAMES
            "mpfr.h"
         HINTS
            "${MPFR_DIR}"
            "$ENV{MPFR_DIR}"
         PATH_SUFFIXES
            "include"
         VALIDATOR
            mpfr_validate
         )

      message(DEBUG "With mpfr_no_default_path=<${mpfr_no_default_path}>, find_path returned mpfr_main_inc <${mpfr_main_inc}>.")
      if(mpfr_main_inc)
         # The VALIDATOR feature does not exist before 3.25, so we just manually call the validate function.
         if(CMAKE_VERSION VERSION_LESS 3.25)
            message(DEBUG "A-posteriori validating candidate include directory <${mpfr_main_inc}>.")
            mpfr_validate(mpfr_validation_result "${mpfr_main_inc}")
            if(mpfr_validation_result)
               break()
            else()
               unset(mpfr_main_inc CACHE)
            endif()
         else()
            break()
         endif()
      endif()
   endforeach()

   if(mpfr_main_inc)
      mpfr_parse_version(mpfr_version "${mpfr_main_inc}")
      set(${CMAKE_FIND_PACKAGE_NAME}_VERSION "${mpfr_version}" PARENT_SCOPE)
      mpfr_find_main_lib(mpfr_main_lib "${mpfr_main_inc}")
      if(mpfr_main_lib)

         # C library target
         if(NOT TARGET MPFR::libmpfr)
            add_library(MPFR::libmpfr UNKNOWN IMPORTED)
            set_target_properties(MPFR::libmpfr PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${mpfr_main_inc}"
               IMPORTED_LOCATION "${mpfr_main_lib}"
               )
         endif()

         # Main interface target
         if(NOT TARGET MPFR::MPFR)
            # Note that we do not create an ALIAS to support CMake version < 3.15;
            # see https://cmake.org/cmake/help/latest/release/3.15.html#commands
            # Moreover, this allows a consuming package to ALIAS it.
            add_library(MPFR::MPFR INTERFACE IMPORTED)
            target_link_libraries(MPFR::MPFR INTERFACE MPFR::libmpfr)
         endif()

         # Deprecated variables refer to a deprecated interface target.
         if(NOT TARGET MPFR_LIBRARIES)
            add_library(MPFR_LIBRARIES INTERFACE IMPORTED)
            if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
               set_target_properties(MPFR_LIBRARIES PROPERTIES
                  DEPRECATION "Do not use ${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS and ${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES anymore, but the MPFR::MPFR target instead.")
            else()
               message(DEBUG "CMake version ${CMAKE_VERSION} is too low to mark MPFR_LIBRARIES target as deprecated!")
            endif()
            target_link_libraries(MPFR_LIBRARIES INTERFACE MPFR::libmpfr)
            set(${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES MPFR_LIBRARIES PARENT_SCOPE)
            get_target_property(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS MPFR::libmpfr INTERFACE_INCLUDE_DIRECTORIES)
            set(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS "${${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS}" PARENT_SCOPE)
         endif()

         if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME} version ${mpfr_version}: ${mpfr_main_inc} and ${mpfr_main_lib}")
         endif()
         set(${CMAKE_FIND_PACKAGE_NAME}_FOUND TRUE PARENT_SCOPE)
         return()
      else()
         set(mpfr_fail_reason ": unexpected failure (mpfr_main_lib=${mpfr_main_lib})")
      endif()
   else()
      set(mpfr_fail_reason "")
   endif()

   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE)
      set(mpfr_fail_version " version ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE}")
   else()
      set(mpfr_fail_version "")
   endif()
   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find required dependency MPFR${mpfr_fail_version}${mpfr_fail_reason}. Try specifying suitable MPFR_DIR!")
   elseif(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(DEBUG "Could NOT find MPFR${mpfr_fail_version}${mpfr_fail_reason}.")
   else()
      message(STATUS "Could NOT find optional dependency MPFR${mpfr_fail_version}${mpfr_fail_reason}. Try specifying suitable MPFR_DIR.")
   endif()
   set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE PARENT_SCOPE)

endfunction()

# Invoke main search.
mpfr_find()

