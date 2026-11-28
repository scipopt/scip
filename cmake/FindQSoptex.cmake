# Finds the QSoptex C library.
#
# == Imported targets ==
#
# This module defines the following IMPORTED targets:
# - QSoptex::QSoptex (interfaces QSoptex::libqsoptex)
# - QSoptex::libqsoptex (C library)
#
# == Result variables ==
#
# This module will set the following variables in your project:
# - QSoptex_FOUND (bool indicating presence of QSoptex::libqsoptex)
# - QSoptex_VERSION (version string)
#
# == Hints ==
#
# QSOPTEX_DIR (CMake or environment variable) is inspected when searching for the library.

include(FindPackageHandleStandardArgs)
if(NOT COMMAND find_package_check_version)
   function(find_package_check_version version output_var HANDLE_VERSION_RANGE RESULT_MESSAGE_VARIABLE result_var)
      if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION)
         message(WARNING "Accepting version ${version} as compatible with ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION}; CMake version >= 3.19 needed to actually check versions.")
      endif()
      set(${output_var} 1 PARENT_SCOPE)
   endfunction()
endif()

# Search for the C library.
function(qsoptex_find_main_lib result search_dir)
   message(DEBUG "    qsoptex_find_main_lib(\"${search_dir}\") with qsoptex_no_default_path=<${qsoptex_no_default_path}>")
   find_library(qsoptex_main_lib
      "qsopt_ex"
      PATHS "${search_dir}" "${search_dir}/../../lib"
      ${qsoptex_no_default_path}
      )
   message(DEBUG "    -> ${qsoptex_main_lib}, stored in ${result}")
   set(${result} "${qsoptex_main_lib}" PARENT_SCOPE)
endfunction()

# Inspect library for version.
function(qsoptex_query_version version_var include_dir main_lib)
   message(DEBUG "qsoptex_query_version(\"${include_dir}\", \"${main_lib}\")")

   file(WRITE "${CMAKE_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/FindQSoptexVersion.c" "
#include \"QSopt_ex_version.h\"

int main(int argc, char** argv)
{
  QSopt_ex_version();

  return 0;
}
"
      )
   enable_language(C)
   try_run(run_result compile_result "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}" "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FindQSoptexVersion.c"
      COMPILE_DEFINITIONS "-I${include_dir}"
      LINK_LIBRARIES "${main_lib}"
      COMPILE_OUTPUT_VARIABLE compile_output
      RUN_OUTPUT_VARIABLE run_output
      )
   if(compile_result AND NOT run_result)
      string(REGEX REPLACE "Using QSopt_ex ([0-9]*)\\.([0-9]*)\\.([0-9]*)\\.([0-9]*).*" "\\1.\\2.\\3.\\4" version "${run_output}")
      set(${version_var} "${version}" PARENT_SCOPE)
   else()
      message(WARNING "Unable to extract QSoptex version. Compilation output:\n${compile_output}Execution output:\n${run_output}")
      set(${version_var} NONE PARENT_SCOPE)
   endif()
endfunction()


# Inspection of a candidate include directory.
function(qsoptex_validate result include_dir)
   message(DEBUG "  qsoptex_validate(\"${include_dir}\") with qsoptex_no_default_path=<${qsoptex_no_default_path}>")

   # Check for library.
   qsoptex_find_main_lib(qsoptex_main_lib "${include_dir}")
   message(DEBUG "  -> qsoptex_main_lib=${qsoptex_main_lib}")
   if(NOT qsoptex_main_lib)
      message(DEBUG "    Encountered ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} version ${version} without a library.")
      set(${result} FALSE PARENT_SCOPE)
      return()
   endif()

   # Parse version.
   qsoptex_query_version(version "${include_dir}" "${qsoptex_main_lib}")

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

   message(DEBUG "    Accepting ${CMAKE_FIND_PACKAGE_NAME} headers at ${include_dir} with library ${qsoptex_main_lib} and version ${version}.")
   set(${result} TRUE PARENT_SCOPE)
   return()
endfunction()

# If QSoptex was already found before, check version.
if(${CMAKE_FIND_PACKAGE_NAME}_FOUND
      AND ${CMAKE_FIND_PACKAGE_NAME}_VERSION
      AND TARGET QSoptex::QSoptex
      AND TARGET QSoptex::libqsoptex
   )
   find_package_check_version(${${CMAKE_FIND_PACKAGE_NAME}_VERSION} qsoptex_version_result HANDLE_VERSION_RANGE)
   if(qsoptex_version_result)
      if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
         message(STATUS "Already found suitable ${CMAKE_FIND_PACKAGE_NAME} version \"${${CMAKE_FIND_PACKAGE_NAME}_VERSION}\".")
      endif()
      unset(qsoptex_version_result)
      return()
   endif()
   unset(qsoptex_version_result)
endif()

# Main search function whose variables are local.
function(qsoptex_find)
   if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(STATUS "Searching for ${CMAKE_FIND_PACKAGE_NAME}.")
   endif()

   foreach(find_local TRUE;FALSE)
      # Whether to find libraries only at the location derived from the header.
      if(find_local)
         set(qsoptex_no_default_path NO_DEFAULT_PATH)
      else()
         set(qsoptex_no_default_path)
      endif()

      find_path(qsoptex_main_inc
         NAMES
            "QSopt_ex_version.h"
         HINTS
            "${QSOPTEX_DIR}"
            "$ENV{QSOPTEX_DIR}"
         PATH_SUFFIXES
            "qsopt_ex"
            "include/qsopt_ex"
         VALIDATOR
            qsoptex_validate
         )

      message(DEBUG "With qsoptex_no_default_path=<${qsoptex_no_default_path}>, find_path returned qsoptex_main_inc <${qsoptex_main_inc}>.")
      if(qsoptex_main_inc)
         # The VALIDATOR feature does not exist before 3.25, so we just manually call the validate function.
         if(CMAKE_VERSION VERSION_LESS 3.25)
            message(DEBUG "A-posteriori validating candidate include directory <${qsoptex_main_inc}>.")
            qsoptex_validate(qsoptex_validation_result "${qsoptex_main_inc}")
            if(qsoptex_validation_result)
               break()
            else()
               unset(qsoptex_main_inc CACHE)
            endif()
         else()
            break()
         endif()
      endif()
   endforeach()

   if(qsoptex_main_inc)
      qsoptex_find_main_lib(qsoptex_main_lib "${qsoptex_main_inc}")
      qsoptex_query_version(qsoptex_version "${qsoptex_main_inc}" "${qsoptex_main_lib}")
      set(${CMAKE_FIND_PACKAGE_NAME}_VERSION "${qsoptex_version}" PARENT_SCOPE)
      set(qsoptex_list_libs "${qsoptex_main_lib}")
      if(qsoptex_main_lib)

         # C library target
         if(NOT TARGET QSoptex::libqsoptex)
            add_library(QSoptex::libqsoptex UNKNOWN IMPORTED)
            set_target_properties(QSoptex::libqsoptex PROPERTIES
               INTERFACE_INCLUDE_DIRECTORIES "${qsoptex_main_inc}"
               IMPORTED_LOCATION "${qsoptex_main_lib}"
               )
         endif()

         # Main interface target
         if(NOT TARGET QSoptex::QSoptex)
            # Note that we do not create an ALIAS to support CMake version < 3.15;
            # see https://cmake.org/cmake/help/latest/release/3.15.html#commands
            # Moreover, this allows a consuming package to ALIAS it.
            add_library(QSoptex::QSoptex INTERFACE IMPORTED)
            target_link_libraries(QSoptex::QSoptex INTERFACE QSoptex::libqsoptex)
         endif()

         # Deprecated variables refer to a deprecated interface target.
         if(NOT TARGET QSOPTEX_LIBRARIES)
            add_library(QSOPTEX_LIBRARIES INTERFACE IMPORTED)
            if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
               set_target_properties(QSOPTEX_LIBRARIES PROPERTIES
                  DEPRECATION "Do not use QSOPTEX_INCLUDE_DIRS and QSOPTEX_LIBRARIES anymore, but the QSoptex::QSoptex target instead.")
            else()
               message(DEBUG "CMake version ${CMAKE_VERSION} is too low to mark QSOPTEX_LIBRARIES target as deprecated!")
            endif()
            target_link_libraries(QSOPTEX_LIBRARIES INTERFACE QSoptex::libqsoptex)
            set(QSOPTEX_LIBRARIES QSOPTEX_LIBRARIES PARENT_SCOPE)
            get_target_property(QSOPTEX_INCLUDE_DIRS QSoptex::libqsoptex INTERFACE_INCLUDE_DIRECTORIES)
            set(QSOPTEX_INCLUDE_DIRS "${QSOPTEX_INCLUDE_DIRS}" PARENT_SCOPE)
         endif()

         if(NOT ${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
            message(STATUS "Found ${CMAKE_FIND_PACKAGE_NAME} version ${qsoptex_version}: ${qsoptex_main_inc} and ${qsoptex_list_libs}")
         endif()
         set(${CMAKE_FIND_PACKAGE_NAME}_FOUND TRUE PARENT_SCOPE)
         return()
      else()
         set(qsoptex_fail_reason ": unexpected failure (qsoptex_main_lib=${qsoptex_main_lib})")
      endif()
   else()
      set(qsoptex_fail_reason "")
   endif()

   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE)
      set(qsoptex_fail_version " version ${${CMAKE_FIND_PACKAGE_NAME}_FIND_VERSION_COMPLETE}")
   else()
      set(qsoptex_fail_version "")
   endif()
   if(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find required dependency QSoptex${qsoptex_fail_version}${qsoptex_fail_reason}. Try specifying suitable QSOPTEX_DIR!")
   elseif(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      message(DEBUG "Could NOT find QSoptex${qsoptex_fail_version}${qsoptex_fail_reason}.")
   else()
      message(STATUS "Could NOT find optional dependency QSoptex${qsoptex_fail_version}${qsoptex_fail_reason}. Try specifying suitable QSOPTEX_DIR.")
   endif()
   set(${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE PARENT_SCOPE)

endfunction()

# Invoke main search.
qsoptex_find()

