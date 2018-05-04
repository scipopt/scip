# - Try to find CLP
# Once done this will define
#  CLP_FOUND - System has CLP
#  CLP_INCLUDE_DIRS - The CLP include directories
#  CLP_LIBRARIES - The libraries needed to use CLP

# check whether environment variable CLP_DIR was set
if(NOT CLP_DIR)
   set(CLP_DIR_TEST $ENV{CLP_DIR})
   if(CLP_DIR_TEST)
      set(CLP_DIR $ENV{CLP_DIR} CACHE PATH "Path to clp build directory")
   endif()
endif()

find_path(CLP_INCLUDE_DIR
          NAMES ClpConfig.h
          PATHS "${CLP_DIR}/include/coin"
                "${CBC_DIR}/include/coin"
                 "/usr/include/coin"
                 "C:\\libs\\clp\\include"
                 "C:\\libs\\cbc\\include"
          )

find_library( COIN_LIBRARY
              NAMES CoinUtils
              PATHS "${CLP_DIR}/lib"
                    "${CBC_DIR}/lib"
                    "/usr/lib"
                    "/usr/lib/coin"
                    "C:\\libs\\clp\\lib"
                    "C:\\libs\\cbc\\lib"
              )

find_library( CLP_LIBRARY
              NAMES Clp
              PATHS "${CLP_DIR}/lib"
                    "${CBC_DIR}/lib"
                    "/usr/lib"
                    "/usr/lib/coin"
                    "C:\\libs\\clp\\lib"
                    "C:\\libs\\cbc\\lib"
              )

set(CLP_INCLUDE_DIRS "${CLP_INCLUDE_DIR}" )
set(CLP_LIBRARIES "${CLP_LIBRARY}" "${COIN_LIBRARY}")


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CLP_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CLP  DEFAULT_MSG
                                  CLP_LIBRARY CLP_INCLUDE_DIR)

mark_as_advanced(CLP_INCLUDE_DIR CLP_LIBRARY)

