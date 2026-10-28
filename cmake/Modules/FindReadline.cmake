# Finds the Readline C library.
#
# == Imported targets ==
#
# This module defines the following IMPORTED targets:
# - Readline::Readline (interfaces Readline::libreadline)
# - Readline::libreadline (C library)
#
# == Result variables ==
#
# This module will set the following variables in your project:
# - Readline_FOUND (bool indicating presence of Readline::libreadline)
# - Readline_VERSION (version string)
#
# == Hints ==
#
# READLINE_DIR (CMake or environment variable) is inspected when searching for the library.

include(FindPackageHandleStandardArgs)
if(NOT COMMAND find_package_check_version)
   function(find_package_check_version version output_var HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE result_var)
      if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION)
         message(WARNING "Accepting version ${version} as compatible with ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION}; CMake version >= 3.19 needed to actually check versions.")
      endif()
      set(${output_var} 1 PARENT_SCOPE)
   endfunction()
endif()

# Parse Readline version from an include directory.
function(readline_parse_version result include_dir)
   file(READ "${include_dir}/readline/readline.h" header_file)
   string(REGEX MATCH "#define[ \t]+RL_VERSION_MAJOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_major ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+RL_VERSION_MINOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_minor ${CMAKE_MATCH_1})
   if(version_major)
      set(${result} "${version_major}.${version_minor}" PARENT_SCOPE)
      return()
   else()
      set(${result} NONE PARENT_SCOPE)
      return()
   endif()
endfunction()

# Search for the C library.
function(readline_find_main_lib result search_dir)
   message(DEBUG "    readline_find_main_lib(\"${search_dir}\") with readline_no_default_path=<${readline_no_default_path}>")
   find_library(readline_main_lib
      NAMES
         "readline"
      PATHS
         "${search_dir}"
         "${search_dir}/../lib"
         "${search_dir}/../lib64"
         "${search_dir}/../lib/${CMAKE_SYSTEM_PROCESSOR}-linux-gnu"
      ${readline_no_default_path}
      )
   message(DEBUG "    -> ${readline_main_lib}, stored in ${result}")
   set(${result} "${readline_main_lib}" PARENT_SCOPE)
endfunction()

# Inspection of a candidate include directory.
function(readline_validate result include_dir)
   message(DEBUG "  readline_validate(\"${include_dir}\") with readline_no_default_path=<${readline_no_default_path}>")

   # Parse version.
   readline_parse_version(version "${include_dir}")
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
   readline_find_main_lib(readline_main_lib "${include_dir}")
   message(DEBUG "  -> readline_main_lib=${readline_main_lib}")
   if(NOT readline_main_lib)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without a library.")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   message(DEBUG "    Accepting ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} with library ${readline_main_lib} and version ${version}.")

   set(${result} TRUE PARENT_SCOPE)
   return()
endfunction()

# If Readline was already found before, check version.
if(${CMAKE_FIND_PACKAGE_NAME}_FOUND
      AND ${CMAKE_FIND_PACKAGE_NAME}_VERSION
      AND TARGET Readline::Readline
      AND TARGET Readline::libreadline
   )
   find_package_check_version(${${CMAKE_FIND_PACKAGE_NAME}_VERSION} readline_version_result HANDLE_VERSION_RANGE)
   if(readline_version_result)
      if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
         message(STATUS "Already found suitable ${CMAKE_FIND_PACKAGE_NAME} version \"${${CMAKE_FIND_PACKAGE_NAME}_VERSION}\".")
      endif()
      unset(readline_version_result)
      return()
   endif()
   unset(readline_version_result)
endif()

# Main search function whose variables are local.
function(readline_find)
   if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Searching for ${CMAKE_FIND_PACKAGE_NAME}.")
   endif()

   foreach(find_local TRUE;FALSE)
      # Whether to find libraries only at the location derived from the header.
      if(find_local)
         set(readline_no_default_path NO_DEFAULT_PATH)
      else()
         set(readline_no_default_path)
      endif()

      find_path(readline_main_inc
         NAMES
            "readline/readline.h"
         HINTS
            "${READLINE_DIR}"
            "$ENV{READLINE_DIR}"
         PATH_SUFFIXES
            "include"
         VALIDATOR
            readline_validate
         )

      message(DEBUG "With readline_no_default_path=<${readline_no_default_path}>, find_path returned readline_main_inc <${readline_main_inc}>.")
      if(readline_main_inc)
         # The VALIDATOR feature does not exist before 3.25, so we just manually call the validate function.
         if(CMAKE_VERSION VERSION_LESS 3.25)
            message(DEBUG "A-posteriori validating candidate include directory <${readline_main_inc}>.")
            readline_validate(readline_validation_result "${readline_main_inc}")
            if(readline_validation_result)
               break()
            else()
               unset(readline_main_inc CACHE)
            endif()
         else()
            break()
         endif()
      endif()
   endforeach()

   if(readline_main_inc)
      readline_parse_version(readline_version "${readline_main_inc}")
      set(${CMAKE_FIND_PACKAGE_NAME}_VERSION "${readline_version}" PARENT_SCOPE)
      readline_find_main_lib(readline_main_lib "${readline_main_inc}")
      set(readline_list_libs "${readline_main_lib}")
      if(readline_main_lib)

         # C library target
         if(NOT TARGET Readline::libreadline)
            add_library(Readline::libreadline UNKNOWN IMPORTED)
            set_target_properties(Readline::libreadline PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${readline_main_inc}"
               IMPORTED_LOCATION "${readline_main_lib}"
               )
         endif()

         # Main interface target
         if(NOT TARGET Readline::Readline)
            # Note that we do not create an ALIAS to support CMake version < 3.15;
            # see https://cmake.org/cmake/help/latest/release/3.15.html#commands
            # Moreover, this allows a consuming package to ALIAS it.
            add_library(Readline::Readline INTERFACE IMPORTED)
            target_link_libraries(Readline::Readline INTERFACE Readline::libreadline)
         endif()

         # Deprecated variables refer to a deprecated interface target.
         if(NOT TARGET Readline_LIBRARIES)
            add_library(Readline_LIBRARIES INTERFACE IMPORTED)
            if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
               set_target_properties(Readline_LIBRARIES PROPERTIES
                  DEPRECATION "Do not use ${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS and ${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES anymore, but the Readline::Readline target instead.")
            else()
               message(DEBUG "CMake version ${CMAKE_VERSION} is too low to mark Readline_LIBRARIES target as deprecated!")
            endif()
            target_link_libraries(Readline_LIBRARIES INTERFACE Readline::libreadline)
            set(${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES Readline_LIBRARIES PARENT_SCOPE)
            set(${CMAKE_FIND_PACKAGE_NAME}_LIBRARY Readline_LIBRARIES PARENT_SCOPE)
            get_target_property(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS Readline::libreadline INTERFACE_INCLUDE_DIRECTORIES)
            set(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS "${${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS}" PARENT_SCOPE)
         endif()

         if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME} version ${readline_version}: ${readline_main_inc} and ${readline_list_libs}")
         endif()
         set(${CMAKE_FIND_PACKAGE_NAME}_FOUND TRUE PARENT_SCOPE)
         return()
      else()
         set(readline_fail_reason ": unexpected failure (readline_main_lib=${readline_main_lib})")
      endif()
   else()
      set(readline_fail_reason "")
   endif()

   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE)
      set(readline_fail_version " version ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE}")
   else()
      set(readline_fail_version "")
   endif()
   set(readline_fail_message "Could NOT find ${CMAKE_FIND_PACKAGE_NAME}${readline_fail_version}${readline_fail_reason}. Try specifying suitable READLINE_DIR!")
   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find required dependency ${CMAKE_FIND_PACKAGE_NAME}${readline_fail_version}${readline_fail_reason}. Try specifying suitable READLINE_DIR!")
   elseif(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Could NOT find ${CMAKE_FIND_PACKAGE_NAME}${readline_fail_version}${readline_fail_reason}.")
   else()
      message(DEBUG "Could NOT find optional dependency ${CMAKE_FIND_PACKAGE_NAME}${readline_fail_version}${readline_fail_reason}. Try specifying suitable READLINE_DIR!")
   endif()
   set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE PARENT_SCOPE)

endfunction()

# Invoke main search.
readline_find()

