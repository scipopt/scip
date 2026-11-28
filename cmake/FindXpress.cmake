# Finds the Xpress C library.
#
# == Imported targets ==
#
# This module defines the following IMPORTED targets:
# - Xpress::Xpress (interfaces Xpress::libxpress)
# - Xpress::libxpress (C library)
#
# == Result variables ==
#
# This module will set the following variables in your project:
# - Xpress_FOUND (bool indicating presence of Xpress::libxpress)
# - Xpress_VERSION (version string)
#
# == Hints ==
#
# XPRESS_DIR (CMake or environment variable) is inspected when searching for the library.

include(FindPackageHandleStandardArgs)
if(NOT COMMAND find_package_check_version)
   function(find_package_check_version version output_var HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE result_var)
      if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION)
         message(WARNING "Accepting version ${version} as compatible with ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION}; CMake version >= 3.19 needed to actually check versions.")
      endif()
      set(${output_var} 1 PARENT_SCOPE)
   endfunction()
endif()

# Parse Xpress version from an include directory.
function(xpress_parse_version result include_dir)
   file(READ "${include_dir}/xprs.h" header_file)
   string(REGEX MATCH "#define[ \t]+XPVERSION_MAJOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_major ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+XPVERSION_MINOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_minor ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+XPVERSION_BUILD[ \t]+([0-9]+)" _ ${header_file})
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
function(xpress_find_main_lib result search_dir)
   message(DEBUG "    xpress_find_main_lib(\"${search_dir}\") with xpress_no_default_path=<${xpress_no_default_path}>")
   find_library(xpress_main_lib
      NAMES
         "xprs"
      PATHS
         "${search_dir}"
         "${search_dir}/../lib"
      ${xpress_no_default_path}
      )
   message(DEBUG "    -> ${xpress_main_lib}, stored in ${result}")
   set(${result} "${xpress_main_lib}" PARENT_SCOPE)
endfunction()

# Inspection of a candidate include directory.
function(xpress_validate result include_dir)
   message(DEBUG "  xpress_validate(\"${include_dir}\") with xpress_no_default_path=<${xpress_no_default_path}>")

   # Parse version.
   xpress_parse_version(version "${include_dir}")
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
   xpress_find_main_lib(xpress_main_lib "${include_dir}")
   message(DEBUG "  -> xpress_main_lib=${xpress_main_lib}")
   if(NOT xpress_main_lib)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without a library.")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   message(DEBUG "    Accepting ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} with library ${xpress_main_lib} and version ${version}.")
   set(${result} TRUE PARENT_SCOPE)
   return()
endfunction()

# If Xpress was already found before, check version.
if(${CMAKE_FIND_PACKAGE_NAME}_FOUND
      AND ${CMAKE_FIND_PACKAGE_NAME}_VERSION
      AND TARGET Xpress::Xpress
      AND TARGET Xpress::libxpress
   )
   find_package_check_version(${${CMAKE_FIND_PACKAGE_NAME}_VERSION} xpress_version_result HANDLE_VERSION_RANGE)
   if(xpress_version_result)
      if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
         message(STATUS "Already found suitable ${CMAKE_FIND_PACKAGE_NAME} version \"${${CMAKE_FIND_PACKAGE_NAME}_VERSION}\".")
      endif()
      unset(xpress_version_result)
      return()
   endif()
   unset(xpress_version_result)
endif()

# Main search function whose variables are local.
function(xpress_find)
   if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Searching for ${CMAKE_FIND_PACKAGE_NAME}.")
   endif()

   foreach(find_local TRUE;FALSE)
      # Whether to find libraries only at the location derived from the header.
      if(find_local)
         set(xpress_no_default_path NO_DEFAULT_PATH)
      else()
         set(xpress_no_default_path)
      endif()

      find_path(xpress_main_inc
         NAMES
            "xprs.h"
         HINTS
            "${XPRESS_DIR}"
            "$ENV{XPRESS_DIR}"
         PATH_SUFFIXES
            "include"
         VALIDATOR
            xpress_validate
         )

      message(DEBUG "With xpress_no_default_path=<${xpress_no_default_path}>, find_path returned xpress_main_inc <${xpress_main_inc}>.")
      if(xpress_main_inc)
         # The VALIDATOR feature does not exist before 3.25, so we just manually call the validate function.
         if(CMAKE_VERSION VERSION_LESS 3.25)
            message(DEBUG "A-posteriori validating candidate include directory <${xpress_main_inc}>.")
            xpress_validate(xpress_validation_result "${xpress_main_inc}")
            if(xpress_validation_result)
               break()
            else()
               unset(xpress_main_inc CACHE)
            endif()
         else()
            break()
         endif()
      endif()
   endforeach()

   if(xpress_main_inc)
      xpress_parse_version(xpress_version "${xpress_main_inc}")
      set(${CMAKE_FIND_PACKAGE_NAME}_VERSION "${xpress_version}" PARENT_SCOPE)
      xpress_find_main_lib(xpress_main_lib "${xpress_main_inc}")
      set(xpress_list_libs "${xpress_main_lib}")
      if(xpress_main_lib)

         # C library target
         if(NOT TARGET Xpress::libxpress)
            add_library(Xpress::libxpress UNKNOWN IMPORTED)
            set_target_properties(Xpress::libxpress PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${xpress_main_inc}"
               IMPORTED_LOCATION "${xpress_main_lib}"
               )
         endif()

         # Main interface target
         if(NOT TARGET Xpress::Xpress)
            # Note that we do not create an ALIAS to support CMake version < 3.15;
            # see https://cmake.org/cmake/help/latest/release/3.15.html#commands
            # Moreover, this allows a consuming package to ALIAS it.
            add_library(Xpress::Xpress INTERFACE IMPORTED)
            target_link_libraries(Xpress::Xpress INTERFACE Xpress::libxpress)
         endif()

         # Deprecated variables refer to a deprecated interface target.
         if(NOT TARGET Xpress_LIBRARIES)
            add_library(XPRESS_LIBRARIES INTERFACE IMPORTED)
            if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
               set_target_properties(XPRESS_LIBRARIES PROPERTIES
                  DEPRECATION "Do not use XPRESS_INCLUDE_DIRS and XPRESS_LIBRARIES anymore, but the Xpress::Xpress target instead.")
            else()
               message(DEBUG "CMake version ${CMAKE_VERSION} is too low to mark XPRESS_LIBRARIES target as deprecated!")
            endif()
            target_link_libraries(XPRESS_LIBRARIES INTERFACE Xpress::libxpress)
            set(XPRESS_LIBRARIES XPRESS_LIBRARIES PARENT_SCOPE)
            get_target_property(XPRESS_INCLUDE_DIRS Xpress::libxpress INTERFACE_INCLUDE_DIRECTORIES)
            set(XPRESS_INCLUDE_DIRS "${XPRESS_INCLUDE_DIRS}" PARENT_SCOPE)
         endif()

         if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME} version ${xpress_version}: ${xpress_main_inc} and ${xpress_list_libs}")
         endif()
         set(${CMAKE_FIND_PACKAGE_NAME}_FOUND TRUE PARENT_SCOPE)
         return()
      else()
         set(xpress_fail_reason ": unexpected failure (xpress_main_lib=${xpress_main_lib})")
      endif()
   else()
      set(xpress_fail_reason "")
   endif()

   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE)
      set(xpress_fail_version " version ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE}")
   else()
      set(xpress_fail_version "")
   endif()
   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find required dependency Xpress${xpress_fail_version}${xpress_fail_reason}. Try specifying suitable XPRESS_DIR!")
   elseif(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(DEBUG "Could NOT find Xpress${xpress_fail_version}${xpress_fail_reason}.")
   else()
      message(STATUS "Could NOT find optional dependency Xpress${xpress_fail_version}${xpress_fail_reason}. Try specifying suitable XPRESS_DIR.")
   endif()
   set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE PARENT_SCOPE)

endfunction()

# Invoke main search.
xpress_find()

