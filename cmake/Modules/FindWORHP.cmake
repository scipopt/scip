# Finds the WORHP C library.
#
# == Imported targets ==
#
# This module defines the following IMPORTED targets:
# - WORHP::WORHP (interfaces WORHP::libworhp)
# - WORHP::libworhp (C library)
#
# == Result variables ==
#
# This module will set the following variables in your project:
# - WORHP_FOUND (bool indicating presence of WORHP::libworhp)
# - WORHP_VERSION (version string)
#
# == Hints ==
#
# WORHP_DIR (CMake or environment variable) is inspected when searching for the library.

include(FindPackageHandleStandardArgs)
if(NOT COMMAND find_package_check_version)
   function(find_package_check_version version output_var HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE result_var)
      if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION)
         message(WARNING "Accepting version ${version} as compatible with ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION}; CMake version >= 3.19 needed to actually check versions.")
      endif()
      set(${output_var} 1 PARENT_SCOPE)
   endfunction()
endif()

# Parse WORHP version from an include directory.
function(worhp_parse_version result include_dir)
   file(READ "${include_dir}/worhp_version.h" header_file)
   string(REGEX MATCH "#[ \t]*define[ \t]+WORHP_VERSION[ \t]+\"([.0-9]+)\"" _ ${header_file})
   set(version ${CMAKE_MATCH_1})
   if(version)
      set(${result} "${version}" PARENT_SCOPE)
      return()
   else()
      set(${result} NONE PARENT_SCOPE)
      return()
   endif()
endfunction()

# Search for the C library.
function(worhp_find_main_lib result search_dir)
   message(DEBUG "    worhp_find_main_lib(\"${search_dir}\") with worhp_no_default_path=<${worhp_no_default_path}>")
   find_library(worhp_main_lib
      NAMES
         "worhp"
      PATHS
         "${search_dir}"
         "${search_dir}/../lib"
         "${search_dir}/../../lib"
      ${worhp_no_default_path}
      )
   message(DEBUG "    -> ${worhp_main_lib}, stored in ${result}")
   set(${result} "${worhp_main_lib}" PARENT_SCOPE)
endfunction()

# Inspection of a candidate include directory.
function(worhp_validate result include_dir)
   message(DEBUG "  worhp_validate(\"${include_dir}\") with worhp_no_default_path=<${worhp_no_default_path}>")

   # Parse version.
   worhp_parse_version(version "${include_dir}")
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
   worhp_find_main_lib(worhp_main_lib "${include_dir}")
   message(DEBUG "  -> worhp_main_lib=${worhp_main_lib}")
   if(NOT worhp_main_lib)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without a library.")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   message(DEBUG "    Accepting ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} with library ${worhp_main_lib} and version ${version}.")
   set(${result} TRUE PARENT_SCOPE)
   return()
endfunction()

# If WORHP was already found before, check version.
if(${CMAKE_FIND_PACKAGE_NAME}_FOUND
      AND ${CMAKE_FIND_PACKAGE_NAME}_VERSION
      AND TARGET WORHP::WORHP
      AND TARGET WORHP::libworhp
   )
   find_package_check_version(${${CMAKE_FIND_PACKAGE_NAME}_VERSION} worhp_version_result HANDLE_VERSION_RANGE)
   if(worhp_version_result)
      if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
         message(STATUS "Already found suitable ${CMAKE_FIND_PACKAGE_NAME} version \"${${CMAKE_FIND_PACKAGE_NAME}_VERSION}\".")
      endif()
      unset(worhp_version_result)
      return()
   endif()
   unset(worhp_version_result)
endif()

# Main search function whose variables are local.
function(worhp_find)
   if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Searching for ${CMAKE_FIND_PACKAGE_NAME}.")
   endif()

   foreach(find_local TRUE;FALSE)
      # Whether to find libraries only at the location derived from the header.
      if(find_local)
         set(worhp_no_default_path NO_DEFAULT_PATH)
      else()
         set(worhp_no_default_path)
      endif()

      find_path(worhp_main_inc
         NAMES
            "worhp.h"
         HINTS
            "${WORHP_DIR}"
            "$ENV{WORHP_DIR}"
         PATH_SUFFIXES
            "worhp"
            "include/worhp"
         VALIDATOR
            worhp_validate
         )

      message(DEBUG "With worhp_no_default_path=<${worhp_no_default_path}>, find_path returned worhp_main_inc <${worhp_main_inc}>.")
      if(worhp_main_inc)
         # The VALIDATOR feature does not exist before 3.25, so we just manually call the validate function.
         if(CMAKE_VERSION VERSION_LESS 3.25)
            message(DEBUG "A-posteriori validating candidate include directory <${worhp_main_inc}>.")
            worhp_validate(worhp_validation_result "${worhp_main_inc}")
            if(worhp_validation_result)
               break()
            else()
               unset(worhp_main_inc CACHE)
            endif()
         else()
            break()
         endif()
      endif()
   endforeach()

   if(worhp_main_inc)
      worhp_parse_version(worhp_version "${worhp_main_inc}")
      set(${CMAKE_FIND_PACKAGE_NAME}_VERSION "${worhp_version}" PARENT_SCOPE)
      worhp_find_main_lib(worhp_main_lib "${worhp_main_inc}")
      set(worhp_list_libs "${worhp_main_lib}")
      if(worhp_main_lib)

         # C library target
         if(NOT TARGET WORHP::libworhp)
            add_library(WORHP::libworhp UNKNOWN IMPORTED)
            set_target_properties(WORHP::libworhp PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${worhp_main_inc}"
               IMPORTED_LOCATION "${worhp_main_lib}"
               )
         endif()

         # Main interface target
         if(NOT TARGET WORHP::WORHP)
            # Note that we do not create an ALIAS to support CMake version < 3.15;
            # see https://cmake.org/cmake/help/latest/release/3.15.html#commands
            # Moreover, this allows a consuming package to ALIAS it.
            add_library(WORHP::WORHP INTERFACE IMPORTED)
            target_link_libraries(WORHP::WORHP INTERFACE WORHP::libworhp)
         endif()

         # Deprecated variables refer to a deprecated interface target.
         if(NOT TARGET WORHP_LIBRARIES)
            add_library(WORHP_LIBRARIES INTERFACE IMPORTED)
            if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
               set_target_properties(WORHP_LIBRARIES PROPERTIES
                  DEPRECATION "Do not use ${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS and ${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES anymore, but the WORHP::WORHP target instead.")
            else()
               message(DEBUG "CMake version ${CMAKE_VERSION} is too low to mark WORHP_LIBRARIES target as deprecated!")
            endif()
            target_link_libraries(WORHP_LIBRARIES INTERFACE WORHP::libworhp)
            set(${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES WORHP_LIBRARIES PARENT_SCOPE)
            get_target_property(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS WORHP::libworhp INTERFACE_INCLUDE_DIRECTORIES)
            set(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS "${${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS}" PARENT_SCOPE)
         endif()

         if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME} version ${worhp_version}: ${worhp_main_inc} and ${worhp_list_libs}")
         endif()
         set(${CMAKE_FIND_PACKAGE_NAME}_FOUND TRUE PARENT_SCOPE)
         return()
      else()
         set(worhp_fail_reason ": unexpected failure (worhp_main_lib=${worhp_main_lib})")
      endif()
   else()
      set(worhp_fail_reason "")
   endif()

   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE)
      set(worhp_fail_version " version ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE}")
   else()
      set(worhp_fail_version "")
   endif()
   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find required dependency WORHP${worhp_fail_version}${worhp_fail_reason}. Try specifying suitable WORHP_DIR!")
   elseif(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(DEBUG "Could NOT find WORHP${worhp_fail_version}${worhp_fail_reason}.")
   else()
      message(STATUS "Could NOT find optional dependency WORHP${worhp_fail_version}${worhp_fail_reason}. Try specifying suitable WORHP_DIR.")
   endif()
   set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE PARENT_SCOPE)

endfunction()

# Invoke main search.
worhp_find()

