cmake_minimum_required(VERSION 3.9)

project(MIPSolver)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(scipmip
   src/cppmain.cpp)

target_link_libraries(scipmip ${SCIP_LIBRARIES})

if( TARGET examples )
    add_dependencies( examples scipmip )
endif()
