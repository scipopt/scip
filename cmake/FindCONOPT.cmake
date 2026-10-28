# Finds the CONOPT C++ library.
#
# == Imported targets ==
#
# This module defines the following IMPORTED targets:
# - CONOPT::CONOPT (interfaces CONOPT::libconopt)
# - CONOPT::libconopt (C++ library)
#
# == Result variables ==
#
# This module will set the following variables in your project:
# - CONOPT_FOUND (bool indicating presence of CONOPT::libconopt)
# - CONOPT_VERSION (version string)
#
# == Hints ==
#
# CONOPT_DIR (CMake or environment variable) is inspected when searching for the library.
# CONOPT_LICENSE_INT_1 (First license integer)
# CONOPT_LICENSE_INT_2 (Second license integer)
# CONOPT_LICENSE_INT_3 (Third license integer)
# CONOPT_LICENSE_TEXT  (License text, e.g., -DCONOPT_LICENSE_TEXT="Zuse Institute Berlin")

include(FindPackageHandleStandardArgs)
if(NOT COMMAND find_package_check_version)
   function(find_package_check_version version output_var HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE result_var)
      if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION)
         message(WARNING "Accepting version ${version} as compatible with ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION}; CMake version >= 3.19 needed to actually check versions.")
      endif()
      set(${output_var} 1 PARENT_SCOPE)
   endfunction()
endif()

# Parse CONOPT version from an include directory.
function(conopt_parse_version result include_dir)
   file(READ "${include_dir}/conopt.h" header_file)
   string(REGEX MATCH "#define[ \t]+CONOPT_VERSION_MAJOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_major ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+CONOPT_VERSION_MINOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_minor ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+CONOPT_VERSION_PATCH[ \t]+([0-9]+)" _ ${header_file})
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
function(conopt_find_main_lib result search_dir)
   message(DEBUG "    conopt_find_main_lib(\"${search_dir}\") with conopt_no_default_path=<${conopt_no_default_path}>")
   find_library(conopt_main_lib
      NAMES
         "conopt"
         "conopt4"
      PATHS
         "${search_dir}"
         "${search_dir}/../lib"
      ${conopt_no_default_path}
      )
   message(DEBUG "    -> ${conopt_main_lib}, stored in ${result}")
   set(${result} "${conopt_main_lib}" PARENT_SCOPE)
endfunction()

# Inspection of a candidate include directory.
function(conopt_validate result include_dir)
   message(DEBUG "  conopt_validate(\"${include_dir}\") with conopt_no_default_path=<${conopt_no_default_path}>")

   # Parse version.
   conopt_parse_version(version "${include_dir}")
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
   conopt_find_main_lib(conopt_main_lib "${include_dir}")
   message(DEBUG "  -> conopt_main_lib=${conopt_main_lib}")
   if(NOT conopt_main_lib)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without a library.")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   message(DEBUG "    Accepting ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} with library ${conopt_main_lib} and version ${version}.")
   set(${result} TRUE PARENT_SCOPE)
   return()
endfunction()

# If CONOPT was already found before, check version.
if(${CMAKE_FIND_PACKAGE_NAME}_FOUND
      AND ${CMAKE_FIND_PACKAGE_NAME}_VERSION
      AND TARGET CONOPT::CONOPT
      AND TARGET CONOPT::libconopt
   )
   find_package_check_version(${${CMAKE_FIND_PACKAGE_NAME}_VERSION} conopt_version_result HANDLE_VERSION_RANGE)
   if(conopt_version_result)
      if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
         message(STATUS "Already found suitable ${CMAKE_FIND_PACKAGE_NAME} version ${${CMAKE_FIND_PACKAGE_NAME}_VERSION}.")
      endif()
      unset(conopt_version_result)
      return()
   endif()
   unset(conopt_version_result)
endif()

# Main search function whose variables are local.
function(conopt_find)
   if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Searching for ${CMAKE_FIND_PACKAGE_NAME}.")
   endif()

   foreach(find_local TRUE;FALSE)
      # Whether to find libraries only at the location derived from the header.
      if(find_local)
         set(conopt_no_default_path NO_DEFAULT_PATH)
      else()
         set(conopt_no_default_path)
      endif()

      find_path(conopt_main_inc
         NAMES
            "conopt.h"
         HINTS
            "${CONOPT_DIR}"
            "$ENV{CONOPT_DIR}"
         PATH_SUFFIXES
            "include"
         VALIDATOR
            conopt_validate
         )

      message(DEBUG "With conopt_no_default_path=<${conopt_no_default_path}>, find_path returned conopt_main_inc <${conopt_main_inc}>.")
      if(conopt_main_inc)
         # The VALIDATOR feature does not exist before 3.25, so we just manually call the validate function.
         if(CMAKE_VERSION VERSION_LESS 3.25)
            message(DEBUG "A-posteriori validating candidate include directory <${conopt_main_inc}>.")
            conopt_validate(conopt_validation_result "${conopt_main_inc}")
            if(conopt_validation_result)
               break()
            else()
               unset(conopt_main_inc CACHE)
            endif()
         else()
            break()
         endif()
      endif()
   endforeach()

   if(conopt_main_inc)
      conopt_parse_version(conopt_version "${conopt_main_inc}")
      set(${CMAKE_FIND_PACKAGE_NAME}_VERSION "${conopt_version}" PARENT_SCOPE)
      conopt_find_main_lib(conopt_main_lib "${conopt_main_inc}")
      set(conopt_list_libs "${conopt_main_lib}")
      if(conopt_main_lib)

         # C++ library target
         if(NOT TARGET CONOPT::libconopt)
            add_library(CONOPT::libconopt UNKNOWN IMPORTED)
            set_target_properties(CONOPT::libconopt PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${conopt_main_inc}"
               IMPORTED_LOCATION "${conopt_main_lib}"
               )
            set(license_present TRUE)
            foreach(suffix INT_1 INT_2 INT_3 TEXT)
               if(NOT DEFINED CONOPT_LICENSE_${suffix})
                  if(DEFINED ENV{CONOPT_LICENSE_${suffix}})
                     set(CONOPT_LICENSE_${suffix} $ENV{CONOPT_LICENSE_${suffix}})
                  else()
                     set(license_present FALSE)
                  endif()
               endif()
            endforeach()
            if(license_present)
               target_compile_definitions(CONOPT::libconopt
                  INTERFACE
                     "CONOPT_LICENSE_INT_1=${CONOPT_LICENSE_INT_1}"
                     "CONOPT_LICENSE_INT_2=${CONOPT_LICENSE_INT_2}"
                     "CONOPT_LICENSE_INT_3=${CONOPT_LICENSE_INT_3}"
                     "CONOPT_LICENSE_TEXT=\"${CONOPT_LICENSE_TEXT}\""
                  )
            else()
               message(WARNING "No CONOPT license information in (environment) variables CONOPT_LICENSE_{INT_1,INT_2,INT_3,TEXT} found.")
            endif()
         endif()


         # Main interface target
         if(NOT TARGET CONOPT::CONOPT)
            # Note that we do not create an ALIAS to support CMake version < 3.15;
            # see https://cmake.org/cmake/help/latest/release/3.15.html#commands
            # Moreover, this allows a consuming package to ALIAS it.
            add_library(CONOPT::CONOPT INTERFACE IMPORTED)
            target_link_libraries(CONOPT::CONOPT INTERFACE CONOPT::libconopt)
         endif()

         # Deprecated variables refer to a deprecated interface target.
         if(NOT TARGET CONOPT_LIBRARIES)
            add_library(CONOPT_LIBRARIES INTERFACE IMPORTED)
            if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
               set_target_properties(CONOPT_LIBRARIES PROPERTIES
                  DEPRECATION "Do not use ${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS and ${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES anymore, but the CONOPT::CONOPT target instead.")
            else()
               message(DEBUG "CMake version ${CMAKE_VERSION} is too low to mark CONOPT_LIBRARIES target as deprecated!")
            endif()
            target_link_libraries(CONOPT_LIBRARIES INTERFACE CONOPT::libconopt)
            set(${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES CONOPT_LIBRARIES PARENT_SCOPE)
            get_target_property(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS CONOPT::libconopt INTERFACE_INCLUDE_DIRECTORIES)
            set(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS "${${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS}" PARENT_SCOPE)
         endif()

         if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME} version ${conopt_version}: ${conopt_main_inc} and ${conopt_list_libs}")
         endif()
         set(${CMAKE_FIND_PACKAGE_NAME}_FOUND TRUE PARENT_SCOPE)
         return()
      else()
         set(conopt_fail_reason ": unexpected failure (conopt_main_lib=${conopt_main_lib})")
      endif()
   else()
      set(conopt_fail_reason "")
   endif()

   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE)
      set(conopt_fail_version " version ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE}")
   else()
      set(conopt_fail_version "")
   endif()
   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find required dependency CONOPT${conopt_fail_version}${conopt_fail_reason}. Try specifying suitable CONOPT_DIR!")
   elseif(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(DEBUG "Could NOT find CONOPT${conopt_fail_version}${conopt_fail_reason}.")
   else()
      message(STATUS "Could NOT find optional dependency CONOPT${conopt_fail_version}${conopt_fail_reason}. Try specifying suitable CONOPT_DIR.")
   endif()
   set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE PARENT_SCOPE)

endfunction()

# Invoke main search.
conopt_find()

