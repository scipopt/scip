cmake_minimum_required(VERSION 3.10)

if(${CMAKE_VERSION} VERSION_LESS "3.12")
  project(
    Bliss
    VERSION 0.73.3
    LANGUAGES CXX C
    DESCRIPTION "A Tool for Computing Automorphism Groups and Canonical Labelings of Graphs. http://www.tcs.hut.fi/Software/bliss/."
  )
else()
  project(
    Bliss
    VERSION 0.73.3
    LANGUAGES CXX C
    DESCRIPTION "A Tool for Computing Automorphism Groups and Canonical Labelings of Graphs"
    HOMEPAGE_URL "http://www.tcs.hut.fi/Software/bliss/"
  )
endif()

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


set(
  BLISS_SOURCE_FILES
  src/bliss_C.cc
  src/defs.cc
  src/graph.cc
  src/heap.cc
  src/orbit.cc
  src/partition.cc
  src/timer.cc
  src/uintseqhash.cc
  src/utils.cc
)
if(WIN32)
  list(REMOVE_ITEM BLISS_SOURCE_FILES src/timer.cc)
endif()

add_library(libbliss STATIC ${BLISS_SOURCE_FILES})

set_target_properties(libbliss PROPERTIES OUTPUT_NAME bliss)
add_library(Bliss::libbliss ALIAS libbliss)

target_include_directories(
  libbliss
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
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
