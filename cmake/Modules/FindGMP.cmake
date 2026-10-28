# Finds the GMP C and C++ libraries.
#
# == Imported targets ==
#
# This module defines the following IMPORTED targets:
# - GMP::GMP (interfaces GMP::libgmp)
# - GMP::libgmp (C library)
# - GMP::libgmpxx (C++ library)
#
# == Result variables ==
#
# This module will set the following variables in your project:
# - GMP_FOUND (bool indicating presence of GMP::libgmp)
# - GMP_cxx_FOUND (bool indicating presence of GMP::libgmpxx)
# - GMP_VERSION (version string)
#
# == Hints ==
#
# GMP_DIR (CMake or environment variable) is inspected when searching for the library.

include(FindPackageHandleStandardArgs)
if(NOT COMMAND find_package_check_version)
   function(find_package_check_version version output_var HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE result_var)
      if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION)
         message(WARNING "Accepting version ${version} as compatible with ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION}; CMake version >= 3.19 needed to actually check versions.")
      endif()
      set(${output_var} 1 PARENT_SCOPE)
   endfunction()
endif()

# Parse GMP version from an include directory.
function(gmp_parse_version result include_dir)
   message(DEBUG "    gmp_parse_version(\"${include_dir}\")")
   file(READ "${include_dir}/gmp.h" header_file)
   string(REGEX MATCH "#define[ \t]+__GNU_MP_VERSION[ \t]+([0-9]+)" _ ${header_file})
   set(version_major ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+__GNU_MP_VERSION_MINOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_minor ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+__GNU_MP_VERSION_PATCHLEVEL[ \t]+([0-9]+)" _ ${header_file})
   set(version_patch ${CMAKE_MATCH_1})
   message(DEBUG "    -> parsed \"${version_major}.${version_minor}.${version_patch}\".")
   if(version_major)
      set(${result} "${version_major}.${version_minor}.${version_patch}" PARENT_SCOPE)
      return()
   endif()

   # Backup: try to compile.
   file(WRITE "${CMAKE_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/FindGMPVersion.c" "
#include \"gmp.h\"

#define STRINGIFY_IMPL(x) #x
#define STRINGIFY(x) STRINGIFY_IMPL(x)

#pragma message(\"VERSION=<\" STRINGIFY(__GNU_MP_VERSION) \".\" STRINGIFY(__GNU_MP_VERSION_MINOR) \".\" STRINGIFY(__GNU_MP_VERSION_PATCHLEVEL) \">\" )

int main(int argc, char** argv)
{
  return 0;
}
"
      )
   message(DEBUG "    -> wrote \"${CMAKE_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/FindGMPVersion.c\".")
   try_compile(compile_result "${CMAKE_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/" "${CMAKE_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/FindGMPVersion.c"
      COMPILE_DEFINITIONS "-I${include_dir}"
      OUTPUT_VARIABLE compile_output
      )
   message(DEBUG "    -> compilation: ${compile_result}")
   if(compile_result)
      string(REGEX REPLACE ".*VERSION=<([0-9]*)\\.([0-9]*)\\.([0-9]*)>.*" "\\1.\\2.\\3" version "${compile_output}")
      message(DEBUG "    -> version \"${version}\"")
      set(${result} "${version}" PARENT_SCOPE)
   else()
      set(${result} NONE PARENT_SCOPE)
   endif()
endfunction()

# Search for the C library.
function(gmp_find_main_lib result search_dir)
   message(DEBUG "    gmp_find_main_lib(\"${search_dir}\") with gmp_no_default_path=<${gmp_no_default_path}>")
   find_library(gmp_main_lib
      NAMES
         "gmp"
      PATHS
         "${search_dir}"
         "${search_dir}/../lib"
         "${search_dir}/../lib64"
         "${search_dir}/../lib/${CMAKE_SYSTEM_PROCESSOR}-linux-gnu"
      ${gmp_no_default_path}
      )
   message(DEBUG "    -> ${gmp_main_lib}, stored in ${result}")
   set(${result} "${gmp_main_lib}" PARENT_SCOPE)
endfunction()

# Search for the C++ headers.
function(gmp_find_cxx_inc result search_dir)
   message(DEBUG "    gmp_find_cxx_inc(\"${search_dir}\")")
   find_path(gmp_cxx_inc
      NAMES
         "gmpxx.h"
      PATHS
         "${search_dir}"
      )
   set(${result} "${gmp_cxx_inc}" PARENT_SCOPE)
endfunction()

# Search for the C++ library.
function(gmp_find_cxx_lib result search_dir)
   message(DEBUG "    gmp_find_cxx_lib(\"${search_dir}\") with gmp_no_default_path=<${gmp_no_default_path}>")
   find_library(gmp_cxx_lib
      NAMES
         "libgmpxx"
      PATHS
         "${search_dir}"
         "${search_dir}/../lib"
         "${search_dir}/../lib64"
         "${search_dir}/../lib/${CMAKE_SYSTEM_PROCESSOR}-linux-gnu"
      ${gmp_no_default_path}
      )
   set(${result} "${gmp_cxx_lib}" PARENT_SCOPE)
endfunction()

# Inspection of a candidate include directory.
function(gmp_validate result include_dir)
   message(DEBUG "  gmp_validate(\"${include_dir}\") with gmp_no_default_path=<${gmp_no_default_path}>")

   # Parse version.
   gmp_parse_version(version "${include_dir}")

   # Check version.
   if(version)
      find_package_check_version(${version} version_okay HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE reason)
   else()
      set(version_okay FALSE)
      set(reason "unable to determine version.")
   endif()
   if(NOT version_okay)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} and library ${gmp_main_lib} version ${lib_version} ignored: ${reason}")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   # Check for library.
   gmp_find_main_lib(gmp_main_lib "${include_dir}")
   message(DEBUG "  -> gmp_main_lib=${gmp_main_lib}")
   if(NOT gmp_main_lib)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without a library.")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   if("cxx" IN_LIST ${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS)
      # Check for C++ headers.
      gmp_find_cxx_inc(gmp_cxx_inc "${include_dir}")
      if(NOT gmp_cxx_inc AND ${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED_cxx)
         message(VERBOSE "Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without C++ headers.")
         set(${result} FALSE PARENT_SCOPE)
         return()
      endif()

      # Check for C++ library.
      gmp_find_cxx_lib(gmp_cxx_lib "${include_dir}")
      if(NOT gmp_cxx_lib AND ${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED_cxx)
         message(VERBOSE "Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without a C++ library.")
         set(${result} FALSE PARENT_SCOPE)
         return()
      endif()
      message(DEBUG "    Accepting ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} with libraries ${gmp_main_lib}, ${gmp_cxx_lib} and version ${version}.")
   else()
      message(DEBUG "    Accepting ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} with library ${gmp_main_lib} and version ${version}.")
   endif()

   set(${result} TRUE PARENT_SCOPE)
   return()
endfunction()

# If GMP was already found before, check version.
if(${CMAKE_FIND_PACKAGE_NAME}_FOUND
      AND ${CMAKE_FIND_PACKAGE_NAME}_VERSION
      AND TARGET GMP::GMP
      AND TARGET GMP::libgmp
      AND (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED_cxx OR TARGET GMP::libgmpxx)
   )
   find_package_check_version(${${CMAKE_FIND_PACKAGE_NAME}_VERSION} gmp_version_result HANDLE_VERSION_RANGE)
   if(gmp_version_result)
      if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
         message(STATUS "Already found suitable ${CMAKE_FIND_PACKAGE_NAME} version ${${CMAKE_FIND_PACKAGE_NAME}_VERSION}.")
      endif()
      unset(gmp_version_result)
      return()
   endif()
   unset(gmp_version_result)
endif()

# Main search function whose variables are local.
function(gmp_find)
   if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Searching for ${CMAKE_FIND_PACKAGE_NAME}.")
   endif()

   foreach(find_local TRUE;FALSE)
      # Whether to find libraries only at the location derived from the header.
      if(find_local)
         set(gmp_no_default_path NO_DEFAULT_PATH)
      else()
         set(gmp_no_default_path)
      endif()

      find_path(gmp_main_inc
         NAMES
            "gmp.h"
         HINTS
            "${GMP_DIR}"
            "$ENV{GMP_DIR}"
         PATH_SUFFIXES
            "include"
         VALIDATOR
            gmp_validate
         )

      message(DEBUG "With gmp_no_default_path=<${gmp_no_default_path}>, find_path returned gmp_main_inc <${gmp_main_inc}>.")
      if(gmp_main_inc)
         # The VALIDATOR feature does not exist before 3.25, so we just manually call the validate function.
         if(CMAKE_VERSION VERSION_LESS 3.25)
            message(DEBUG "A-posteriori validating candidate include directory <${gmp_main_inc}>.")
            gmp_validate(gmp_validation_result "${gmp_main_inc}")
            if(gmp_validation_result)
               break()
            else()
               unset(gmp_main_inc CACHE)
            endif()
         else()
            break()
         endif()
      endif()
   endforeach()

   if(gmp_main_inc)
      gmp_parse_version(gmp_version "${gmp_main_inc}")
      set(${CMAKE_FIND_PACKAGE_NAME}_VERSION "${gmp_version}" PARENT_SCOPE)
      gmp_find_main_lib(gmp_main_lib "${gmp_main_inc}")
      set(gmp_libs "${gmp_main_lib}")
      if("cxx" IN_LIST ${CMAKE_FIND_PACKAGE_NAME}_FIND_COMPONENTS)
         gmp_find_cxx_inc(gmp_cxx_inc "${gmp_main_inc}")
         gmp_find_cxx_lib(libgmpxx "${gmp_main_inc}")
         set(gmp_libs "${gmp_libs},${gmp_cxx_lib}")
      endif()
      if(gmp_main_lib
            AND (NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED_cxx OR (gmp_cxx_inc AND gmp_cxx_lib))
         )
         # C library target
         if(NOT TARGET GMP::libgmp)
            add_library(GMP::libgmp UNKNOWN IMPORTED)
            set_target_properties(GMP::libgmp PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${gmp_main_inc}"
               IMPORTED_LOCATION "${gmp_main_lib}"
               )
         endif()

         # Main interface target
         if(NOT TARGET GMP::GMP)
            # Note that we do not create an ALIAS to support CMake version < 3.15;
            # see https://cmake.org/cmake/help/latest/release/3.15.html#commands
            # Moreover, this allows a consuming package to ALIAS it.
            add_library(GMP::GMP INTERFACE IMPORTED)
            target_link_libraries(GMP::GMP INTERFACE GMP::libgmp)
         endif()

         # C++ library target
         if(gmp_cxx_inc AND gmp_cxx_lib AND NOT TARGET GMP::libgmpxx)
            add_library(GMP::libgmpxx UNKNOWN IMPORTED)
            set_target_properties(GMP::libgmpxx PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${gmp_cxx_inc}"
               IMPORTED_LOCATION "${gmp_cxx_lib}"
               INTERFACE_LINK_LIBRARIES "GMP::libgmp"
               )
         endif()

         # Deprecated variables refer to a deprecated interface target.
         if(NOT TARGET GMP_LIBRARIES)
            add_library(GMP_LIBRARIES INTERFACE IMPORTED)
            if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
               set_target_properties(GMP_LIBRARIES PROPERTIES
                  DEPRECATION "Do not use ${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS and ${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES anymore, but the GMP::GMP target instead.")
            else()
               message(DEBUG "CMake version ${CMAKE_VERSION} is too low to mark GMP_LIBRARIES target as deprecated!")
            endif()
            target_link_libraries(GMP_LIBRARIES INTERFACE GMP::libgmp)
            set(${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES GMP_LIBRARIES PARENT_SCOPE)
            get_target_property(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS GMP::libgmp INTERFACE_INCLUDE_DIRECTORIES)
            set(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS "${${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS}" PARENT_SCOPE)
         endif()

         if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME} version ${gmp_version}: ${gmp_main_inc} and ${gmp_libs}")
         endif()
         set(${CMAKE_FIND_PACKAGE_NAME}_FOUND TRUE PARENT_SCOPE)
         return()
      else()
         set(gmp_fail_reason ": unexpected failure (gmp_main_lib=${gmp_main_lib}, gmp_cxx_inc=${gmp_cxx_inc}, gmp_cxx_lib=${gmp_cxx_lib})")
      endif()
   else()
      set(gmp_fail_reason "")
   endif()

   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE)
      set(gmp_fail_version " version ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE}")
   else()
      set(gmp_fail_version "")
   endif()
   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find required dependency GMP${gmp_fail_version}${gmp_fail_reason}. Try specifying suitable GMP_DIR!")
   elseif(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(DEBUG "Could NOT find GMP${gmp_fail_version}${gmp_fail_reason}.")
   else()
      message(STATUS "Could NOT find optional dependency GMP${gmp_fail_version}${gmp_fail_reason}. Try specifying suitable GMP_DIR.")
   endif()
   set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE PARENT_SCOPE)

endfunction()

# Invoke main search.
gmp_find()

