find_path(CPLEX_INCLUDE_DIRS
    NAMES cplex.h
    HINTS ${CPLEX_DIR} $ENV{CPLEX_DIR}
    PATH_SUFFIXES include/ilcplex include)

if(MSVC)
   string(REGEX REPLACE "/VC/bin/.*" "" VISUAL_STUDIO_PATH ${CMAKE_CXX_COMPILER})
   string(REGEX MATCH "Studio/[0-9]+/" CPLEX_WIN_VS_VERSION ${VISUAL_STUDIO_PATH})
   string(REGEX REPLACE "Studio/" "" CPLEX_WIN_VS_VERSION ${CPLEX_WIN_VS_VERSION})
   string(REGEX REPLACE "/" "" CPLEX_WIN_VS_VERSION ${CPLEX_WIN_VS_VERSION})

   if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
      if(MT)
         set(CPLEX_WIN_RUNTIME mtd)
      else(MT)
         set(CPLEX_WIN_RUNTIME mta)
      endif(MT)
   else()
      if(MT)
         set(CPLEX_WIN_RUNTIME mdd)
      else(MT)
         set(CPLEX_WIN_RUNTIME mda)
      endif(MT)
   endif()

   find_library(CPLEX_LIBRARY
      NAMES cplex2212 cplex2211 cplex2210 cplex2010 cplex
      HINTS ${CPLEX_DIR} $ENV{CPLEX_DIR}
      PATH_SUFFIXES lib/x64_windows_vs${CPLEX_WIN_VS_VERSION}/stat_${CPLEX_WIN_RUNTIME})

else(MSVC)

   find_library(CPLEX_LIBRARY
      NAMES cplex
      HINTS ${CPLEX_DIR} $ENV{CPLEX_DIR}
      PATH_SUFFIXES lib/x86-64_linux/static_pic
                    lib/x86-64_osx/static_pic
                    lib/arm64_linux/static_pic
                    lib/arm64_osx/static_pic
                    lib)

endif(MSVC)

# todo properly check when pthread is necessary
set(CPLEX_LIBRARIES ${CPLEX_LIBRARY} pthread ${CMAKE_DL_LIBS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CPLEX DEFAULT_MSG CPLEX_INCLUDE_DIRS CPLEX_LIBRARIES)
