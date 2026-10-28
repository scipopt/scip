# Finds the CoinUtils C++ library.
#
# == Imported targets ==
#
# This module defines the following IMPORTED targets:
# - CoinUtils::CoinUtils (interfaces CoinUtils::libcoinutils)
# - CoinUtils::libcoinutils (C++ library)
#
# == Result variables ==
#
# This module will set the following variables in your project:
# - COINUTILS_FOUND (bool indicating presence of CoinUtils::libcoinutils)
# - COINUTILS_VERSION (version string)
#
# == Hints ==
#
# COINUTILS_DIR (CMake or environment variable) is inspected when searching for the library.

include(FindPackageHandleStandardArgs)
if(NOT COMMAND find_package_check_version)
   function(find_package_check_version version output_var HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE result_var)
      if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION)
         message(WARNING "Accepting version ${version} as compatible with ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION}; CMake version >= 3.19 needed to actually check versions.")
      endif()
      set(${output_var} 1 PARENT_SCOPE)
   endfunction()
endif()

# Search for the C++ library.
function(coinutils_find_main_lib result search_dir)
   message(DEBUG "    coinutils_find_main_lib(\"${search_dir}\") with coinutils_no_default_path=<${coinutils_no_default_path}>")
   find_library(coinutils_main_lib
      NAMES
         "CoinUtils"
      PATHS
         "${search_dir}/../../lib/"
      ${coinutils_no_default_path}
      )
   message(DEBUG "    -> ${coinutils_main_lib}, stored in ${result}")
   set(${result} "${coinutils_main_lib}" PARENT_SCOPE)
endfunction()

# Inspection of a candidate include directory.
function(coinutils_validate result include_dir)
   message(DEBUG "  coinutils_validate(\"${include_dir}\") with coinutils_no_default_path=<${coinutils_no_default_path}>")

   # Check for library.
   coinutils_find_main_lib(coinutils_main_lib "${include_dir}")
   message(DEBUG "  -> coinutils_main_lib=${coinutils_main_lib}")
   if(NOT coinutils_main_lib)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} without a library.")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   message(DEBUG "    Accepting ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} with library ${coinutils_main_lib}.")
   set(${result} TRUE PARENT_SCOPE)
   return()
endfunction()

# If CoinUtils was already found before, check version.
if(${CMAKE_FIND_PACKAGE_NAME}_FOUND
      AND ${CMAKE_FIND_PACKAGE_NAME}_VERSION
      AND TARGET CoinUtils::CoinUtils
      AND TARGET CoinUtils::libcoinutils
   )
   if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Already found suitable ${CMAKE_FIND_PACKAGE_NAME} version.")
   endif()
   return()
endif()

# Main search function whose variables are local.
function(coinutils_find)
   if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Searching for ${CMAKE_FIND_PACKAGE_NAME}.")
   endif()

   foreach(find_local TRUE;FALSE)
      # Whether to find libraries only at the location derived from the header.
      if(find_local)
         set(coinutils_no_default_path NO_DEFAULT_PATH)
      else()
         set(coinutils_no_default_path)
      endif()

      find_path(coinutils_main_inc
         NAMES
            "CoinUtilsConfig.h"
         HINTS
            "${COINUTILS_DIR}"
            "$ENV{COINUTILS_DIR}"
            "${COIN_DIR}"
            "$ENV{COIN_DIR}"
         PATH_SUFFIXES
            "coin-or"
            "include/coin-or"
            "dist/include/coin-or"
         VALIDATOR
            coinutils_validate
         )

      message(DEBUG "With coinutils_no_default_path=<${coinutils_no_default_path}>, find_path returned coinutils_main_inc <${coinutils_main_inc}>.")
      if(coinutils_main_inc)
         # The VALIDATOR feature does not exist before 3.25, so we just manually call the validate function.
         if(CMAKE_VERSION VERSION_LESS 3.25)
            message(DEBUG "A-posteriori validating candidate include directory <${coinutils_main_inc}>.")
            coinutils_validate(coinutils_validation_result "${coinutils_main_inc}")
            if(coinutils_validation_result)
               break()
            else()
               unset(coinutils_main_inc CACHE)
            endif()
         else()
            break()
         endif()
      endif()
   endforeach()

   if(coinutils_main_inc)
      coinutils_find_main_lib(coinutils_main_lib "${coinutils_main_inc}")
      if(coinutils_main_lib)

         # C++ library target
         if(NOT TARGET CoinUtils::libcoinutils)
            add_library(CoinUtils::libcoinutils UNKNOWN IMPORTED)
            set_target_properties(CoinUtils::libcoinutils PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${coinutils_main_inc}"
               IMPORTED_LOCATION "${coinutils_main_lib}"
               )
         endif()

         # Main interface target
         if(NOT TARGET CoinUtils::CoinUtils)
            # Note that we do not create an ALIAS to support CMake version < 3.15;
            # see https://cmake.org/cmake/help/latest/release/3.15.html#commands
            # Moreover, this allows a consuming package to ALIAS it.
            add_library(CoinUtils::CoinUtils INTERFACE IMPORTED)
            target_link_libraries(CoinUtils::CoinUtils INTERFACE CoinUtils::libcoinutils)
         endif()

         # Deprecated variables refer to a deprecated interface target.
         if(NOT TARGET CoinUtils_LIBRARIES)
            add_library(CoinUtils_LIBRARIES INTERFACE IMPORTED)
            if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
               set_target_properties(CoinUtils_LIBRARIES PROPERTIES
                  DEPRECATION "Do not use ${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS and ${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES anymore, but the CoinUtils::CoinUtils target instead.")
            else()
               message(DEBUG "CMake version ${CMAKE_VERSION} is too low to mark CoinUtils_LIBRARIES target as deprecated!")
            endif()
            target_link_libraries(CoinUtils_LIBRARIES INTERFACE CoinUtils::libcoinutils)
            set(${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES CoinUtils_LIBRARIES PARENT_SCOPE)
            get_target_property(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS CoinUtils::libcoinutils INTERFACE_INCLUDE_DIRECTORIES)
            set(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS "${${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS}" PARENT_SCOPE)
         endif()

         if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME}: ${coinutils_main_inc} and ${coinutils_main_lib}")
         endif()
         set(${CMAKE_FIND_PACKAGE_NAME}_FOUND TRUE PARENT_SCOPE)
         return()
      else()
         set(coinutils_fail_reason ": unexpected failure (coinutils_main_lib=${coinutils_main_lib})")
      endif()
   else()
      set(coinutils_fail_reason "")
   endif()

   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE)
      set(coinutils_fail_version " version ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE}")
   else()
      set(coinutils_fail_version "")
   endif()
   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find required dependency CoinUtils${coinutils_fail_version}${coinutils_fail_reason}. Try specifying suitable COINUTILS_DIR!")
   elseif(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(DEBUG "Could NOT find CoinUtils${coinutils_fail_version}${coinutils_fail_reason}.")
   else()
      message(STATUS "Could NOT find optional dependency CoinUtils${coinutils_fail_version}${coinutils_fail_reason}. Try specifying suitable COINUTILS_DIR.")
   endif()
   set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE PARENT_SCOPE)

endfunction()

# Invoke main search.
coinutils_find()

