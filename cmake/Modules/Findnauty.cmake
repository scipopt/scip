# Finds the nauty C library.
#
# == Imported targets ==
#
# This module defines the following IMPORTED targets:
# - nauty::nauty (interfaces nauty::libnauty)
# - nauty::libnauty (C library)
#
# == Result variables ==
#
# This module will set the following variables in your project:
# - nauty_FOUND (bool indicating presence of nauty::libnauty)
# - nauty_VERSION (version string)
#
# == Hints ==
#
# NAUTY_DIR (CMake or environment variable) is inspected when searching for the library.

include(FindPackageHandleStandardArgs)
if(NOT COMMAND find_package_check_version)
   function(find_package_check_version version output_var HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE result_var)
      if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION)
         message(WARNING "Accepting version ${version} as compatible with ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION}; CMake version >= 3.19 needed to actually check versions.")
      endif()
      set(${output_var} 1 PARENT_SCOPE)
   endfunction()
endif()

# Parse nauty version from an include directory.
function(nauty_parse_version result include_dir)
   file(READ "${include_dir}/nauty/nauty.h" header_file)
   string(REGEX MATCH "[ \t]*#[ \t]*define NAUTYVERSION[ \t]*\"([0-9.]+) \\(.*" _ "${header_file}")
   if(CMAKE_MATCH_1)
      set(${result} "${CMAKE_MATCH_1}" PARENT_SCOPE)
      return()
   else()
      set(${result} NONE PARENT_SCOPE)
      return()
   endif()
endfunction()

# Search for the C library.
function(nauty_find_main_lib result search_dir)
   message(DEBUG "    nauty_find_main_lib(\"${search_dir}\") with nauty_no_default_path=<${nauty_no_default_path}>")
   find_library(nauty_main_lib
      NAMES "nauty"
      PATHS "${search_dir}" "${search_dir}/../lib" "${search_dir}/../lib64" "${search_dir}/../lib/x86_64-linux-gnu"
      ${nauty_no_default_path}
      )
   message(DEBUG "    -> ${nauty_main_lib}, stored in ${result}")
   set(${result} "${nauty_main_lib}" PARENT_SCOPE)
endfunction()

# Inspection of a candidate include directory.
function(nauty_validate result include_dir)
   message(DEBUG "  nauty_validate(\"${include_dir}\") with nauty_no_default_path=<${nauty_no_default_path}>")

   # Parse version.
   nauty_parse_version(version "${include_dir}")
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
   nauty_find_main_lib(nauty_main_lib "${include_dir}")
   message(DEBUG "  -> nauty_main_lib=${nauty_main_lib}")
   if(NOT nauty_main_lib)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without a library.")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   message(DEBUG "    Accepting ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} with library ${nauty_main_lib} and version ${version}.")
   set(${result} TRUE PARENT_SCOPE)
   return()
endfunction()

# If nauty was already found before, check version.
if(${CMAKE_FIND_PACKAGE_NAME}_FOUND
      AND ${CMAKE_FIND_PACKAGE_NAME}_VERSION
      AND TARGET nauty::nauty
      AND TARGET nauty::libnauty
   )
   find_package_check_version(${${CMAKE_FIND_PACKAGE_NAME}_VERSION} nauty_version_result HANDLE_VERSION_RANGE)
   if(nauty_version_result)
      if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
         message(STATUS "Already found suitable ${CMAKE_FIND_PACKAGE_NAME} version \"${${CMAKE_FIND_PACKAGE_NAME}_VERSION}\".")
      endif()
      unset(nauty_version_result)
      return()
   endif()
   unset(nauty_version_result)
endif()

# Main search function whose variables are local.
function(nauty_find)
   if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Searching for ${CMAKE_FIND_PACKAGE_NAME}.")
   endif()

   foreach(find_local TRUE;FALSE)
      # Whether to find libraries only at the location derived from the header.
      if(find_local)
         set(nauty_no_default_path NO_DEFAULT_PATH)
      else()
         set(nauty_no_default_path)
      endif()

      find_path(nauty_main_inc
         NAMES
            "nauty/nauty.h"
         HINTS
            "${NAUTY_DIR}"
            "$ENV{NAUTY_DIR}"
         PATH_SUFFIXES
            "include"
         VALIDATOR
            nauty_validate
         )

      message(DEBUG "With nauty_no_default_path=<${nauty_no_default_path}>, find_path returned nauty_main_inc <${nauty_main_inc}>.")
      if(nauty_main_inc)
         # The VALIDATOR feature does not exist before 3.25, so we just manually call the validate function.
         if(CMAKE_VERSION VERSION_LESS 3.25)
            message(DEBUG "A-posteriori validating candidate include directory <${nauty_main_inc}>.")
            nauty_validate(nauty_validation_result "${nauty_main_inc}")
            if(nauty_validation_result)
               break()
            else()
               unset(nauty_main_inc CACHE)
            endif()
         else()
            break()
         endif()
      endif()
   endforeach()

   if(nauty_main_inc)
      nauty_parse_version(nauty_version "${nauty_main_inc}")
      set(${CMAKE_FIND_PACKAGE_NAME}_VERSION "${nauty_version}" PARENT_SCOPE)
      nauty_find_main_lib(nauty_main_lib "${nauty_main_inc}")
      set(nauty_list_libs "${nauty_main_lib}")
      if(nauty_main_lib)
 
         # C library target
         if(NOT TARGET nauty::libnauty)
            add_library(nauty::libnauty UNKNOWN IMPORTED)
            set_target_properties(nauty::libnauty PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${nauty_main_inc}"
               IMPORTED_LOCATION "${nauty_main_lib}"
               )
         endif()

         # Main interface target
         if(NOT TARGET nauty::nauty)
            # Note that we do not create an ALIAS to support CMake version < 3.15;
            # see https://cmake.org/cmake/help/latest/release/3.15.html#commands
            # Moreover, this allows a consuming package to ALIAS it.
            add_library(nauty::nauty INTERFACE IMPORTED)
            target_link_libraries(nauty::nauty INTERFACE nauty::libnauty)
         endif()

         # Deprecated variables refer to a deprecated interface target.
         if(NOT TARGET nauty_LIBRARIES)
            add_library(nauty_LIBRARIES INTERFACE IMPORTED)
            if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
               set_target_properties(nauty_LIBRARIES PROPERTIES
                  DEPRECATION "Do not use ${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS and ${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES anymore, but the nauty::nauty target instead.")
            else()
               message(DEBUG "CMake version ${CMAKE_VERSION} is too low to mark nauty_LIBRARIES target as deprecated!")
            endif()
            target_link_libraries(nauty_LIBRARIES INTERFACE nauty::libnauty)
            set(${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES nauty_LIBRARIES PARENT_SCOPE)
            get_target_property(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS nauty::libnauty INTERFACE_INCLUDE_DIRECTORIES)
            set(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS "${${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS}" PARENT_SCOPE)
         endif()

         if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME} version ${nauty_version}: ${nauty_main_inc} and ${nauty_list_libs}")
         endif()
         set(${CMAKE_FIND_PACKAGE_NAME}_FOUND TRUE PARENT_SCOPE)
         return()
      else()
         set(nauty_fail_reason ": unexpected failure (nauty_main_lib=${nauty_main_lib})")
      endif()
   else()
      set(nauty_fail_reason "")
   endif()

   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE)
      set(nauty_fail_version " version ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE}")
   else()
      set(nauty_fail_version "")
   endif()
   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find required dependency nauty${nauty_fail_version}${nauty_fail_reason}. Try specifying suitable NAUTY_DIR!")
   elseif(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(DEBUG "Could NOT find nauty${nauty_fail_version}${nauty_fail_reason}.")
   else()
      message(STATUS "Could NOT find optional dependency nauty${nauty_fail_version}${nauty_fail_reason}. Try specifying suitable NAUTY_DIR.")
   endif()
   set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE PARENT_SCOPE)

endfunction()

# Invoke main search.
nauty_find()

