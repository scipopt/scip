cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME PolySCIP)
project(${PROJECT_NAME})

# set variables
set(EXECUTABLE_NAME "polyscip")
set(POLYSCIP_VERSION_MAJOR 2)
set(POLYSCIP_VERSION_MINOR 0)
set(POLYSCIP_MAX_NUMBER_OBJS "30" CACHE STRING "Maximum number of objectives for PolySCIP instances")
set(SRC src)
set(DOC doc)
set(CMAKE cmake)
set(CMAKE_BUILD_TYPE Debug)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/${CMAKE})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Werror -pedantic -Wno-unused-variable")


# configure a header file to pass CMake settings to source code
configure_file("${PROJECT_SOURCE_DIR}/${CMAKE}/PolySCIPConfig.h.in" "${PROJECT_SOURCE_DIR}/${SRC}/PolySCIPConfig.h" ESCAPE_QUOTES)

include(FindPackageHandleStandardArgs)
find_package(SCIP REQUIRED)
#find_package(SOPLEX REQUIRED)
find_package(ZLIB)
find_package(GMP)
find_package(Readline)


include_directories(${PROJECT_SOURCE_DIR}/${SRC})
include_directories(${SCIP_INC})

add_subdirectory(${SRC})
add_subdirectory(${DOC})

