# Finds the Gurobi C and C++ libraries.
#
# == Imported targets ==
#
# This module defines the following IMPORTED targets:
# - Gurobi::Gurobi (interfaces Gurobi::libgurobi)
# - Gurobi::libgurobi (C library)
# - Gurobi::libgurobixx (C++ library)
#
# == Result variables ==
#
# This module will set the following variables in your project:
# - Gurobi_FOUND (bool indicating presence of Gurobi::libgurobi)
# - Gurobi_cxx_FOUND (bool indicating presence of Gurobi::libgurobixx)
# - Gurobi_VERSION (version string)
#
# == Hints ==
#
# GUROBI_DIR (CMake or environment variable) is inspected when searching for the library.

include(FindPackageHandleStandardArgs)
if(NOT COMMAND find_package_check_version)
   function(find_package_check_version version output_var HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE result_var)
      if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION)
         message(WARNING "Accepting version ${version} as compatible with ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION}; CMake version >= 3.19 needed to actually check versions.")
      endif()
      set(${output_var} 1 PARENT_SCOPE)
   endfunction()
endif()

# Parse Gurobi version from an include directory.
function(gurobi_parse_version result include_dir)
   file(READ "${include_dir}/gurobi_c.h" header_file)
   string(REGEX MATCH "#define[ \t]+GRB_VERSION_MAJOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_major ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+GRB_VERSION_MINOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_minor ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+GRB_VERSION_TECHNICAL[ \t]+([0-9]+)" _ ${header_file})
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
function(gurobi_find_main_lib result search_dir version)
   message(DEBUG "    gurobi_find_main_lib(\"${search_dir}\", \"${version}\") with gurobi_no_default_path=<${gurobi_no_default_path}>")
   string(REPLACE "." ";" version_list "${version}")
   list(GET version_list 0 version_major)
   list(GET version_list 1 version_minor)
   find_library(gurobi_main_lib
      NAMES
         "gurobi${version_major}${version_minor}"
      PATHS
         "${search_dir}/../lib"
      ${gurobi_no_default_path}
      )
   message(DEBUG "    -> ${gurobi_main_lib}, stored in ${result}")
   set(${result} "${gurobi_main_lib}" PARENT_SCOPE)
endfunction()

# Search for the C++ headers.
function(gurobi_find_cxx_inc result search_dir)
   message(DEBUG "    gurobi_find_cxx_inc(\"${search_dir}\")")
   find_path(gurobi_cxx_inc
      NAMES
         "gurobi_c++.h"
      PATHS
         "${search_dir}"
      )
   set(${result} "${gurobi_cxx_inc}" PARENT_SCOPE)
endfunction()

# Search for the C++ library.
function(gurobi_find_cxx_lib result search_dir)
   message(DEBUG "    gurobi_find_cxx_lib(\"${search_dir}\") with gurobi_no_default_path=<${gurobi_no_default_path}>")
   find_library(gurobi_cxx_lib
      NAMES
         "gurobi_c++"
      PATHS
         "${search_dir}/../lib"
      ${gurobi_no_default_path}
      )
   set(${result} "${gurobi_cxx_lib}" PARENT_SCOPE)
endfunction()

# Inspection of a candidate include directory.
function(gurobi_validate result include_dir)
   message(DEBUG "  gurobi_validate(\"${include_dir}\") with gurobi_no_default_path=<${gurobi_no_default_path}>")

   # Parse version.
   gurobi_parse_version(version "${include_dir}")
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
   gurobi_find_main_lib(gurobi_main_lib "${include_dir}" "${version}")
   message(DEBUG "  -> gurobi_main_lib=${gurobi_main_lib}")
   if(NOT gurobi_main_lib)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without a library.")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   if("cxx" IN_LIST ${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS)
      # Check for C++ headers.
      gurobi_find_cxx_inc(gurobi_cxx_inc "${include_dir}")
      message(DEBUG "  -> gurobi_cxx_inc=${gurobi_cxx_inc}")
      if(NOT gurobi_cxx_inc AND ${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED_cxx)
         message(VERBOSE "Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without C++ headers.")
         set(${result} FALSE PARENT_SCOPE)
         return()
      endif()

      # Check for C++ library.
      gurobi_find_cxx_lib(gurobi_cxx_lib "${include_dir}")
      message(DEBUG "  -> gurobi_cxx_lib=${gurobi_cxx_lib}")
      if(NOT gurobi_cxx_lib AND ${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED_cxx)
         message(VERBOSE "Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without a C++ library.")
         set(${result} FALSE PARENT_SCOPE)
         return()
      endif()
      message(DEBUG "    Accepting ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} with libraries ${gurobi_main_lib}, ${gurobi_cxx_lib} and version ${version}.")
   else()
      message(DEBUG "    Accepting ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} with library ${gurobi_main_lib} and version ${version}.")
   endif()

   set(${result} TRUE PARENT_SCOPE)
   return()
endfunction()

# If Gurobi was already found before, check version.
if(${CMAKE_FIND_PACKAGE_NAME}_FOUND
      AND ${CMAKE_FIND_PACKAGE_NAME}_VERSION
      AND TARGET Gurobi::Gurobi
      AND TARGET Gurobi::libgurobi
      AND (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED_cxx OR TARGET Gurobi::libgurobixx)
   )
   find_package_check_version(${${CMAKE_FIND_PACKAGE_NAME}_VERSION} gurobi_version_result HANDLE_VERSION_RANGE)
   if(gurobi_version_result)
      if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
         message(STATUS "Already found suitable ${CMAKE_FIND_PACKAGE_NAME} version ${${CMAKE_FIND_PACKAGE_NAME}_VERSION}.")
      endif()
      unset(gurobi_version_result)
      return()
   endif()
   unset(gurobi_version_result)
endif()

# Main search function whose variables are local.
function(gurobi_find)
   if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Searching for ${CMAKE_FIND_PACKAGE_NAME}.")
   endif()

   foreach(find_local TRUE;FALSE)
      # Whether to find libraries only at the location derived from the header.
      if(find_local)
         set(gurobi_no_default_path NO_DEFAULT_PATH)
      else()
         set(gurobi_no_default_path)
      endif()

      find_path(gurobi_main_inc
         NAMES
            "gurobi_c.h"
         HINTS
            "${GUROBI_DIR}"
            "$ENV{GUROBI_DIR}"
         PATH_SUFFIXES
            "include"
            "linux64/include"
         VALIDATOR
            gurobi_validate
         )

      message(DEBUG "With gurobi_no_default_path=<${gurobi_no_default_path}>, find_path returned gurobi_main_inc <${gurobi_main_inc}>.")
      if(gurobi_main_inc)
         # The VALIDATOR feature does not exist before 3.25, so we just manually call the validate function.
         if(CMAKE_VERSION VERSION_LESS 3.25)
            message(DEBUG "A-posteriori validating candidate include directory <${gurobi_main_inc}>.")
            gurobi_validate(gurobi_validation_result "${gurobi_main_inc}")
            if(gurobi_validation_result)
               break()
            else()
               unset(gurobi_main_inc CACHE)
            endif()
         else()
            break()
         endif()
      endif()
   endforeach()

   if(gurobi_main_inc)
      gurobi_parse_version(gurobi_version "${gurobi_main_inc}")
      set(${CMAKE_FIND_PACKAGE_NAME}_VERSION "${gurobi_version}" PARENT_SCOPE)
      gurobi_find_main_lib(gurobi_main_lib "${gurobi_main_inc}" "${gurobi_version}")
      set(gurobi_list_libs "${gurobi_main_lib}")
      if("cxx" IN_LIST ${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS)
         gurobi_find_cxx_inc(gurobi_cxx_inc "${gurobi_main_inc}")
         gurobi_find_cxx_lib(libgurobixx "${gurobi_main_inc}")
         set(gurobi_list_libs "${gurobi_list_libs},${gurobi_cxx_lib}")
      endif()
      if(gurobi_main_lib
            AND (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED_cxx OR (gurobi_cxx_inc AND gurobi_cxx_lib))
         )
         # C library target
         if(NOT TARGET Gurobi::libgurobi)
            add_library(Gurobi::libgurobi UNKNOWN IMPORTED)
            set_target_properties(Gurobi::libgurobi PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${gurobi_main_inc}"
               IMPORTED_LOCATION "${gurobi_main_lib}"
               )
         endif()

         # Main interface target
         if(NOT TARGET Gurobi::Gurobi)
            # Note that we do not create an ALIAS to support CMake version < 3.15;
            # see https://cmake.org/cmake/help/latest/release/3.15.html#commands
            # Moreover, this allows a consuming package to ALIAS it.
            add_library(Gurobi::Gurobi INTERFACE IMPORTED)
            target_link_libraries(Gurobi::Gurobi INTERFACE Gurobi::libgurobi)
         endif()

         # C++ library target
         if(gurobi_cxx_inc AND gurobi_cxx_lib AND NOT TARGET Gurobi::libgurobixx)
            add_library(Gurobi::libgurobixx UNKNOWN IMPORTED)
            set_target_properties(Gurobi::libgurobixx PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${gurobi_cxx_inc}"
               IMPORTED_LOCATION "${gurobi_cxx_lib}"
               INTERFACE_LINK_LIBRARIES "Gurobi::libgurobi"
               )
         endif()

         # Deprecated variables refer to a deprecated interface target.
         if(NOT TARGET GUROBI_LIBRARIES)
            add_library(GUROBI_LIBRARIES INTERFACE IMPORTED)
            if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
               set_target_properties(GUROBI_LIBRARIES PROPERTIES
                  DEPRECATION "Do not use GUROBI_INCLUDE_DIRS and GUROBI_LIBRARIES anymore, but the Gurobi::Gurobi target instead.")
            else()
               message(DEBUG "CMake version ${CMAKE_VERSION} is too low to mark GUROBI_LIBRARIES target as deprecated!")
            endif()
            target_link_libraries(GUROBI_LIBRARIES INTERFACE Gurobi::libgurobi)
            set(GUROBI_LIBRARIES GUROBI_LIBRARIES PARENT_SCOPE)
            get_target_property(GUROBI_INCLUDE_DIRS Gurobi::libgurobi INTERFACE_INCLUDE_DIRECTORIES)
            set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIRS}" PARENT_SCOPE)
         endif()

         if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME} version ${gurobi_version}: ${gurobi_main_inc} and ${gurobi_list_libs}")
         endif()
         set(${CMAKE_FIND_PACKAGE_NAME}_FOUND TRUE PARENT_SCOPE)
         return()
      else()
         set(gurobi_fail_reason ": unexpected failure (gurobi_main_lib=${gurobi_main_lib}, gurobi_cxx_inc=${gurobi_cxx_inc}, gurobi_cxx_lib=${gurobi_cxx_lib})")
      endif()
   else()
      set(gurobi_fail_reason "")
   endif()

   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE)
      set(gurobi_fail_version " version ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE}")
   else()
      set(gurobi_fail_version "")
   endif()
   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find required dependency Gurobi${gurobi_fail_version}${gurobi_fail_reason}. Try specifying suitable GUROBI_DIR!")
   elseif(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(DEBUG "Could NOT find Gurobi${gurobi_fail_version}${gurobi_fail_reason}.")
   else()
      message(STATUS "Could NOT find optional dependency Gurobi${gurobi_fail_version}${gurobi_fail_reason}. Try specifying suitable GUROBI_DIR.")
   endif()
   set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE PARENT_SCOPE)

endfunction()

# Invoke main search.
gurobi_find()

