# Finds the Ipopt C++ library.
#
# == Imported targets ==
#
# This module defines the following IMPORTED targets:
# - Ipopt::Ipopt (interfaces Ipopt::libipopt)
# - Ipopt::libipopt (C++ library)
#
# == Result variables ==
#
# This module will set the following variables in your project:
# - Ipopt_FOUND (bool indicating presence of Ipopt::libipopt)
# - Ipopt_VERSION (version string)
#
# == Hints ==
#
# IPOPT_DIR (CMake or environment variable) is inspected when searching for the library.

include(FindPackageHandleStandardArgs)
if(NOT COMMAND find_package_check_version)
   function(find_package_check_version version output_var HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE result_var)
      if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION)
         message(WARNING "Accepting version ${version} as compatible with ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION}; CMake version >= 3.19 needed to actually check versions.")
      endif()
      set(${output_var} 1 PARENT_SCOPE)
   endfunction()
endif()

# Parse Ipopt version from an include directory.
function(ipopt_parse_version result include_dir)
   file(READ "${include_dir}/IpoptConfig.h" header_file)
   string(REGEX MATCH "#define[ \t]+IPOPT_VERSION_MAJOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_major ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+IPOPT_VERSION_MINOR[ \t]+([0-9]+)" _ ${header_file})
   set(version_minor ${CMAKE_MATCH_1})
   string(REGEX MATCH "#define[ \t]+IPOPT_VERSION_RELEASE[ \t]+([0-9]+)" _ ${header_file})
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
function(ipopt_find_main_lib result search_dir)
   message(DEBUG "    ipopt_find_main_lib(\"${search_dir}\") with ipopt_no_default_path=<${ipopt_no_default_path}>")
   find_library(ipopt_main_lib
      NAMES
         "ipopt"
         "libipopt"
         "ipopt.dll"
         "ipopt-3"
         "ipopt-0"
         "libipopt-3"
         "libipopt-0"
      PATHS
         "${search_dir}/../../lib"
         "${search_dir}/../lib"
      ${ipopt_no_default_path}
      )
   message(DEBUG "    -> ${ipopt_main_lib}, stored in ${result}")
   set(${result} "${ipopt_main_lib}" PARENT_SCOPE)
endfunction()

# Inspection of a candidate include directory.
function(ipopt_validate result include_dir)
   message(DEBUG "  ipopt_validate(\"${include_dir}\") with ipopt_no_default_path=<${ipopt_no_default_path}>")

   # Parse version.
   ipopt_parse_version(version "${include_dir}")
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
   ipopt_find_main_lib(ipopt_main_lib "${include_dir}")
   message(DEBUG "  -> ipopt_main_lib=${ipopt_main_lib}")
   if(NOT ipopt_main_lib)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without a library.")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   message(DEBUG "    Accepting ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} with library ${ipopt_main_lib} and version ${version}.")
   set(${result} TRUE PARENT_SCOPE)
   return()
endfunction()

# If Ipopt was already found before, check version.
if(${CMAKE_FIND_PACKAGE_NAME}_FOUND
      AND ${CMAKE_FIND_PACKAGE_NAME}_VERSION
      AND TARGET Ipopt::Ipopt
      AND TARGET Ipopt::libipopt
   )
   find_package_check_version(${${CMAKE_FIND_PACKAGE_NAME}_VERSION} ipopt_version_result HANDLE_VERSION_RANGE)
   if(ipopt_version_result)
      if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
         message(STATUS "Already found suitable ${CMAKE_FIND_PACKAGE_NAME} version \"${${CMAKE_FIND_PACKAGE_NAME}_VERSION}\".")
      endif()
      unset(ipopt_version_result)
      return()
   endif()
   unset(ipopt_version_result)
endif()

# Main search function whose variables are local.
function(ipopt_find)
   if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Searching for ${CMAKE_FIND_PACKAGE_NAME}.")
   endif()

   # Indent messages for dependency search.
   list(APPEND CMAKE_MESSAGE_INDENT "  ")

   # First search via PkgConfig.
   if(NOT MSVC)
      find_package(PkgConfig)
      if(PKG_CONFIG_FOUND)
         if(Ipopt_FIND_VERSION_EXACT)
            set(_find_ipopt_target "ipopt=${Ipopt_FIND_VERSION}")
         elseif(Ipopt_FIND_VERSION)
            set(_find_ipopt_target "ipopt>=${Ipopt_FIND_VERSION}")
         else()
            set(_find_ipopt_target "ipopt")
         endif()
         set(ENV{PKG_CONFIG_PATH} "${IPOPT_DIR}/lib/pkgconfig/;$ENV{IPOPT_DIR}/lib/pkgconfig/;$ENV{PKG_CONFIG_PATH}")

         pkg_check_modules(_PC_Ipopt QUIET
            IMPORTED_TARGET ${_find_ipopt_target})

         # Detecting Ipopt via PkgConfig only works properly from version 3.13 on.
         if(_PC_Ipopt_FOUND AND _PC_Ipopt_VERSION VERSION_GREATER_EQUAL 3.13)
            add_library(Ipopt::libipopt INTERFACE IMPORTED)
            target_link_libraries(Ipopt::libipopt INTERFACE PkgConfig::_PC_Ipopt)
            set(ipopt_version ${_PC_Ipopt_VERSION})
         endif()
      endif()
   endif()

   # Otherwise, main search via headers.
   if(NOT TARGET Ipopt::libipopt)
      foreach(find_local TRUE;FALSE)
         # Whether to find libraries only at the location derived from the header.
         if(find_local)
            set(ipopt_no_default_path NO_DEFAULT_PATH)
         else()
            set(ipopt_no_default_path)
         endif()

         find_path(ipopt_main_inc
            NAMES
               "IpoptConfig.h"
            HINTS
               "${IPOPT_DIR}"
               "$ENV{IPOPT_DIR}"
               /opt/homebrew
            PATH_SUFFIXES
               "include/coin"
               "include/coin-or"
            VALIDATOR
               ipopt_validate
            )

         message(DEBUG "With ipopt_no_default_path=<${ipopt_no_default_path}>, find_path returned ipopt_main_inc <${ipopt_main_inc}>.")
         if(ipopt_main_inc)
            # The VALIDATOR feature does not exist before 3.25, so we just manually call the validate function.
            if(CMAKE_VERSION VERSION_LESS 3.25)
               message(DEBUG "A-posteriori validating candidate include directory <${ipopt_main_inc}>.")
               ipopt_validate(ipopt_validation_result "${ipopt_main_inc}")
               if(ipopt_validation_result)
                  break()
               else()
                  unset(ipopt_main_inc CACHE)
               endif()
            else()
               break()
            endif()
         endif()
      endforeach()
   endif()

   if(ipopt_main_inc)
      ipopt_parse_version(ipopt_version "${ipopt_main_inc}")
      set(${CMAKE_FIND_PACKAGE_NAME}_VERSION "${ipopt_version}" PARENT_SCOPE)
      ipopt_find_main_lib(ipopt_main_lib "${ipopt_main_inc}")
      set(ipopt_list_libs "${ipopt_main_lib}")
      if(ipopt_main_lib)

         # C++ library target
         if(NOT TARGET Ipopt::libipopt)
            add_library(Ipopt::libipopt UNKNOWN IMPORTED)
            set_target_properties(Ipopt::libipopt PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${ipopt_main_inc}"
               IMPORTED_LOCATION "${ipopt_main_lib}"
               )
         endif()

         # Check for an Ipopt-specific dependency file.
         find_file(ipopt_dep_file
            NAMES
               "ipopt_addlibs_cpp.txt"
            HINTS
               "${ipopt_include_dir}/../../"
               "/usr"
               "/opt/homebrew"
            PATH_SUFFIXES
               "share/coin/doc/Ipopt"
               "share/doc/coin/Ipopt"
               "share/coin-or/doc/Ipopt"
               "share/doc/coin-or/Ipopt"
            NO_DEFAULT_PATH
            )
         mark_as_advanced(ipopt_dep_file)
   
         if(ipopt_dep_file)
            # Parse the Ipopt-specific dependency file.
            set(ipopt_deps_message "using ${ipopt_dep_file}")
            file(READ "${ipopt_dep_file}" ipopt_dependencies)
            string(REGEX REPLACE "(\n|-lipopt)" "" ipopt_dependencies "${ipopt_dependencies}")
            separate_arguments(ipopt_dependencies)
   
            foreach(ipopt_dependency ${ipopt_dependencies})
               message(DEBUG "  Analyzing extra dependency <${ipopt_dependency}>")
               if(ipopt_dependency MATCHES "-L.*")
                  string(REPLACE "-L" "" ipopt_dependency "${ipopt_dependency}")
                  target_link_directories(Ipopt::libipopt
                     INTERFACE "${ipopt_dependency}")
               elseif(ipopt_dependency MATCHES "-l.*")
                  string(REPLACE "-l" "" ipopt_dependency "${ipopt_dependency}")
                  target_link_libraries(Ipopt::libipopt
                     INTERFACE "${ipopt_dependency}")
               endif()
            endforeach()
            set(ipopt_deps_message " using dependencies from ${ipopt_dep_file}")
         else()
            set(ipopt_deps_message "")
         endif()
      else()
         set(ipopt_fail_reason ": unexpected failure (ipopt_main_lib=${ipopt_main_lib})")
      endif()
   else()
      set(ipopt_fail_reason "")
   endif()

   if(TARGET Ipopt::libipopt)
      # Main interface target
      if(NOT TARGET Ipopt::Ipopt)
         # Note that we do not create an ALIAS to support CMake version < 3.15;
         # see https://cmake.org/cmake/help/latest/release/3.15.html#commands
         # Moreover, this allows a consuming package to ALIAS it.
         add_library(Ipopt::Ipopt INTERFACE IMPORTED)
         target_link_libraries(Ipopt::Ipopt INTERFACE Ipopt::libipopt)
      endif()

      # Deprecated variables refer to a deprecated interface target.
      if(NOT TARGET IPOPT_LIBRARIES)
         add_library(IPOPT_LIBRARIES INTERFACE IMPORTED)
         if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
            set_target_properties(IPOPT_LIBRARIES PROPERTIES
               DEPRECATION "Do not use ${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS and ${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES anymore, but the Ipopt::Ipopt target instead.")
         else()
            message(DEBUG "CMake version ${CMAKE_VERSION} is too low to mark IPOPT_LIBRARIES target as deprecated!")
         endif()
         target_link_libraries(IPOPT_LIBRARIES INTERFACE Ipopt::libipopt)
         set(IPOPT_LIBRARIES IPOPT_LIBRARIES PARENT_SCOPE)
         get_target_property(IPOPT_INCLUDE_DIRS Ipopt::libipopt INTERFACE_INCLUDE_DIRECTORIES)
         set(IPOPT_INCLUDE_DIRS "${IPOPT_INCLUDE_DIRS}" PARENT_SCOPE)
      endif()
   endif()

   # Unindent messages again; equivalent to list(POP_BACK, CMAKE_MESSAGE_INDENT) which requires CMake 3.15.
   list(REVERSE CMAKE_MESSAGE_INDENT)
   list(REMOVE_AT CMAKE_MESSAGE_INDENT 0)
   list(REVERSE CMAKE_MESSAGE_INDENT)

   # Emit success message.
   if(TARGET Ipopt::libipopt)
      if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
         if(ipopt_main_inc)
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME} version ${ipopt_version}: ${ipopt_main_inc} and ${ipopt_list_libs}${ipopt_deps_message}")
         else()
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME} version ${ipopt_version}: pkgconfig ${_PC_Ipopt_PREFIX}")
         endif()
      endif()
      set(${CMAKE_FIND_PACKAGE_NAME}_FOUND TRUE PARENT_SCOPE)
      set(${CMAKE_FIND_PACKAGE_NAME}_VERSION "${ipopt_version}" PARENT_SCOPE)
      return()
   endif()

   # Emit failure message.
   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE)
      set(ipopt_fail_version " version ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE}")
   else()
      set(ipopt_fail_version "")
   endif()
   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find required dependency Ipopt${ipopt_fail_version}${ipopt_fail_reason}. Try specifying suitable IPOPT_DIR!")
   elseif(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(DEBUG "Could NOT find Ipopt${ipopt_fail_version}${ipopt_fail_reason}.")
   else()
      message(STATUS "Could NOT find optional dependency Ipopt${ipopt_fail_version}${ipopt_fail_reason}. Try specifying suitable IPOPT_DIR.")
   endif()
   set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE PARENT_SCOPE)

endfunction()

# Invoke main search.
ipopt_find()

