# - Try to find QSO
# Once done this will define
#  QSO_FOUND - System has QSO
#  QSO_INCLUDE_DIRS - The QSO include directories
#  QSO_LIBRARIES - The libraries needed to use QSO

# check whether environment variable QSO_DIR was set
if(NOT QSO_DIR)
   set(QSO_DIR_TEST $ENV{QSO_DIR})
   if(QSO_DIR_TEST)
      set(QSO_DIR $ENV{QSO_DIR} CACHE PATH "Path to qsopt build directory")
   endif()
endif()

find_path(QSO_INCLUDE_DIR
          NAMES qsopt.h
          PATHS "${QSO_DIR}"
          )

find_library( QSO_LIBRARY
              NAMES qsopt
              PATHS "${QSO_DIR}"
              )

set(QSO_INCLUDE_DIRS "${QSO_INCLUDE_DIR}" )
set(QSO_LIBRARIES "${QSO_LIBRARY}" )


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set QSO_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(QSO  DEFAULT_MSG
                                  QSO_LIBRARY QSO_INCLUDE_DIR)

mark_as_advanced(QSO_INCLUDE_DIR QSO_LIBRARY)

