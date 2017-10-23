cmake_minimum_required(VERSION 3.0)

project(MIPSolver)

find_package(SCIP REQUIRED)
include_directories(${SCIP_INCLUDE_DIRS})

add_executable(scipmip
   src/cppmain.cpp)

target_link_libraries(scipmip ${SCIP_LIBRARIES})

if( TARGET examples )
    add_dependencies( examples scipmip )
endif()
