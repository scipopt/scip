cmake_minimum_required(VERSION 3.9)

if(${CMAKE_VERSION} VERSION_LESS "3.12")
  project(
    Bliss
    VERSION 0.77
    LANGUAGES CXX C
    DESCRIPTION "A Tool for Computing Automorphism Groups and Canonical Labelings of Graphs. http://www.tcs.hut.fi/Software/bliss/."
  )
else()
  project(
    Bliss
    VERSION 0.77
    LANGUAGES CXX C
    DESCRIPTION "A Tool for Computing Automorphism Groups and Canonical Labelings of Graphs"
    HOMEPAGE_URL "http://www.tcs.hut.fi/Software/bliss/"
  )
endif()

option(USE_GMP "Use GNU Multiple Precision Arithmetic library" OFF)


set(
  BLISS_SOURCE_FILES
  src/abstractgraph.cc
  src/bliss_C.cc
  src/defs.cc
  src/digraph.cc
  src/graph.cc
  src/orbit.cc
  src/partition.cc
  src/uintseqhash.cc
  src/utils.cc
)

# In MSVC, the executable and libraries should go to the top directory
if(MSVC)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})
endif(MSVC)

# Add the shared library
add_library(libbliss ${BLISS_SOURCE_FILES})

set_target_properties(libbliss PROPERTIES OUTPUT_NAME bliss )
add_library(Bliss::libbliss ALIAS libbliss)

target_include_directories(
  libbliss
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

set_target_properties(libbliss PROPERTIES
  # All code ending in a shared library should be made PIC
  POSITION_INDEPENDENT_CODE ON
  # Compiling with hidden visibility
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
)

add_executable(bliss src/bliss.cc)
add_executable(Bliss::bliss ALIAS bliss)
target_link_libraries(bliss PRIVATE libbliss)
if(USE_GMP)
  target_link_libraries(bliss ${GMP_LIBRARIES})
endif(USE_GMP)

set_target_properties(libbliss bliss PROPERTIES COMPILE_OPTIONS "")

# specify the C++ standard
set_target_properties(libbliss bliss
    PROPERTIES CXX_STANDARD 11
               CXX_STANDARD_REQUIRED True)

if (MSVC)
    # /permissive- for standard C++ with "and" and "or" logical operators
    target_compile_options(libbliss PUBLIC /permissive-)
    target_compile_options(bliss PUBLIC /permissive-)
else()
    # Warnings, optimization, no assertions
    target_compile_options(libbliss PRIVATE -O3 -DNDEBUG) #-Wextra -Werror ### disable as part of scip build
    target_compile_options(bliss PRIVATE -O3 -DNDEBUG) #-Wextra -Werror ### disable as part of scip build
endif()

if(USE_GMP)
  find_path(GMP_INCLUDE_DIR NAMES gmp.h)
  find_library(GMP_LIBRARIES NAMES gmp libgmp mpir REQUIRED)
  if (MSVC)
      target_compile_options(libbliss PUBLIC /DBLISS_USE_GMP /I${GMP_INCLUDE_DIR})
      target_compile_options(bliss PUBLIC /DBLISS_USE_GMP /I${GMP_INCLUDE_DIR})
  else()
      target_compile_options(libbliss PUBLIC -DBLISS_USE_GMP -I${GMP_INCLUDE_DIR})
      target_compile_options(bliss PUBLIC -DBLISS_USE_GMP -I${GMP_INCLUDE_DIR})
  endif()
endif(USE_GMP)

# Set the default build type to the given value if no build type was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to Release as none was specified")
  set(
    CMAKE_BUILD_TYPE Release
    CACHE STRING "Choose the type of build" FORCE
  )
  # Set the possible values of build type for cmake-gui, ccmake
  set_property(
    CACHE CMAKE_BUILD_TYPE
    PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo"
  )
endif()


install(
  TARGETS libbliss bliss
  EXPORT Bliss
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

if(${CMAKE_VERSION} VERSION_LESS "3.14")
  install(
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/bliss"
    DESTINATION include
    COMPONENT lib
  )
else()
  install(
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/bliss"
    TYPE INCLUDE
    COMPONENT lib
  )
endif()

install(EXPORT Bliss
  DESTINATION lib/cmake/Bliss
  NAMESPACE Bliss::
  FILE BlissConfig.cmake
)
