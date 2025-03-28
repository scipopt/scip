# - Try to find QSOPT_ex
# Once done this will define
#  QSOEX_FOUND - System has QSO
#  QSOEX_INCLUDE_DIRS - The QSO include directories
#  QSOEX_LIBRARIES - The libraries needed to use QSO

# check whether environment variable QSO_DIR was set
if(NOT QSOEX_DIR)
   set(QSOEX_DIR_TEST $ENV{QSOEX_DIR})
   if(QSOEX_DIR_TEST)
      set(QSOEX_DIR $ENV{QSOEX_DIR} CACHE PATH "Path to qsopt_ex build directory")
   endif()
endif()

find_path(QSOEX_INCLUDE_DIR
          NAMES QSopt_ex.h
          PATHS "${QSOEX_DIR}/include"
          )

find_library(QSOEX_LIBRARY
             NAMES QSopt_ex
             PATHS "${QSOEX_DIR}/lib"
             )

set(QSOEX_INCLUDE_DIRS "${QSOEX_INCLUDE_DIR}" )
set(QSOEX_LIBRARIES "${QSOEX_LIBRARY}" )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set QSOEX_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(QSOEX  DEFAULT_MSG
                                  QSOEX_LIBRARY QSOEX_INCLUDE_DIR)

mark_as_advanced(QSOEX_INCLUDE_DIR QSOEX_LIBRARY)
