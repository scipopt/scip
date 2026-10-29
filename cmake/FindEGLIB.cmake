# - Try to find EGlib
# Once done this will define
#  EGLIB_FOUND - System has QSO
#  EGLIB_INCLUDE_DIRS - The QSO include directories
#  EGLIB_LIBRARIES - The libraries needed to use QSO

# check whether environment variable QSO_DIR was set
if(NOT EGLIB_DIR)
   set(EGLIB_DIR_TEST $ENV{EGLIB_DIR})
   if(EGLIB_DIR_TEST)
      set(EGLIB_DIR $ENV{EGLIB_DIR} CACHE PATH "Path to qsopt_ex build directory")
   endif()
endif()

find_path(EGLIB_INCLUDE_DIR
          NAMES EGlib.h
          PATHS "${EGLIB_DIR}/include"
          )

find_library(EGLIB_LIBRARY
             NAMES EGlib
             PATHS "${EGLIB_DIR}/lib"
             )

set(EGLIB_INCLUDE_DIRS "${EGLIB_INCLUDE_DIR}" )
set(EGLIB_LIBRARIES "${EGLIB_LIBRARY}" )


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set EGLIB_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(EGLIB  DEFAULT_MSG
                                  EGLIB_LIBRARY EGLIB_INCLUDE_DIR)

mark_as_advanced(EGLIB_INCLUDE_DIR EGLIB_LIBRARY)