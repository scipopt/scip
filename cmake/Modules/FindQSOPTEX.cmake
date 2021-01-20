# - Try to find QSOPT_ex
# Once done this will define
#  QSOPTEX_FOUND - System has QSO
#  QSOPTEX_INCLUDE_DIRS - The QSO include directories
#  QSOPTEX_LIBRARIES - The libraries needed to use QSO

# check whether environment variable QSO_DIR was set
if(NOT QSOPTEX_DIR)
   set(QSOPTEX_DIR_TEST $ENV{QSOPTEX_DIR})
   if(QSOPTEX_DIR_TEST)
      set(QSOPTEX_DIR $ENV{QSOPTEX_DIR} CACHE PATH "Path to qsopt_ex build directory")
   endif()
endif()

find_path(QSOPTEX_INCLUDE_DIR
          NAMES QSopt_ex.h
          PATHS "${QSOPTEX_DIR}/include"
          )

find_library(QSOPTEX_LIBRARY
             NAMES QSopt_ex
             PATHS "${QSOPTEX_DIR}/lib"
             )

set(QSOPTEX_INCLUDE_DIRS "${QSOPTEX_INCLUDE_DIR}" )
set(QSOPTEX_LIBRARIES "${QSOPTEX_LIBRARY}" )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set QSOPTEX_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(QSOPTEX  DEFAULT_MSG
                                  QSOPTEX_LIBRARY QSOPTEX_INCLUDE_DIR)

mark_as_advanced(QSOPTEX_INCLUDE_DIR QSOPTEX_LIBRARY)