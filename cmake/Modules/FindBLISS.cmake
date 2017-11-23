include(FindPackageHandleStandardArgs)

# check whether environment variable BLISS_DIR was set
if(NOT BLISS_DIR)
   set(BLISS_DIR_TEST $ENV{BLISS_DIR})
   if(BLISS_DIR_TEST)
      set(BLISS_DIR $ENV{BLISS_DIR} CACHE PATH "Path to bliss build directory")
   endif()
endif()

# if the bliss directory is specified, first try to use exactly that bliss
if(BLISS_DIR)
   # look for the includes with subdirectory bliss
   find_path(BLISS_INCLUDE_DIR
       NAMES bliss/graph.hh
       PATHS ${BLISS_DIR}
       PATH_SUFFIXES include
       NO_DEFAULT_PATH
       )

   # if not found look for the includes without bliss subdirectory
   if(NOT BLISS_INCLUDE_DIR)
      find_path(BLISS_INCLUDE_DIR
          NAMES graph.hh
          PATHS ${BLISS_DIR}
          PATH_SUFFIXES include
          NO_DEFAULT_PATH
          )

      # if we found the headers there we copy the folder to a bliss folder in the binary dir and use that as include
      if(BLISS_INCLUDE_DIR)
         set(COPY_BLISS_HEADERS TRUE)
      endif()
   endif()

   # look for the library in the bliss directory
   find_library(BLISS_LIBRARY
       NAMES bliss
       PATHS ${BLISS_DIR}
       PATH_SUFFIXES lib
       NO_DEFAULT_PATH
       )

   # set variables and call handle standard args
   set(BLISS_LIBRARIES ${BLISS_LIBRARY})
   set(BLISS_INCLUDE_DIRS ${BLISS_INCLUDE_DIR})

   find_package_handle_standard_args(BLISS DEFAULT_MSG BLISS_INCLUDE_DIRS BLISS_LIBRARIES)

   if(BLISS_FOUND AND COPY_BLISS_HEADERS)
      file(GLOB bliss_headers ${BLISS_INCLUDE_DIR}/*.hh)
      file(COPY ${bliss_headers} DESTINATION ${CMAKE_BINARY_DIR}/bliss)
      set(BLISS_INCLUDE_DIR ${CMAKE_BINARY_DIR} CACHE PATH "Include path for bliss headers" FORCE)
      set(BLISS_INCLUDE_DIRS ${BLISS_INCLUDE_DIR})
   endif()
endif()

# if bliss is not already found by the code above we look for it including system directories
if(NOT BLISS_FOUND)
   find_path(BLISS_INCLUDE_DIR
       NAMES bliss/graph.hh
       PATH_SUFFIXES include)

   find_library(BLISS_LIBRARY
       NAMES bliss
       PATH_SUFFIXES lib)

   set(BLISS_LIBRARIES ${BLISS_LIBRARY})
   set(BLISS_INCLUDE_DIRS ${BLISS_INCLUDE_DIR})

   find_package_handle_standard_args(BLISS DEFAULT_MSG BLISS_INCLUDE_DIRS BLISS_LIBRARIES)
endif()

if(BLISS_FOUND AND ${BLISS_LIBRARY} MATCHES "bliss\\.so.*$")
   # check whether the bliss library with built with gmp, because then we need to add flags for the build
   include(GetPrerequisites)
   get_prerequisites(${BLISS_LIBRARY} BLISS_PREREQUISITES 0 0 "" "")
   foreach(prerequisit ${BLISS_PREREQUISITES})
      if(prerequisit MATCHES "gmp")
         set(BLISS_CXX_FLAGS "-DBLISS_USE_GMP" CACHE STRING "Extra CXX flags required for bliss")
      endif()
   endforeach()
endif()
