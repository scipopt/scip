cmake_minimum_required(VERSION 3.9)

project(Sudoku)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(sudoku
   src/sudoku_main.cpp
   src/sudoku_utils.h)

target_link_libraries(sudoku ${SCIP_LIBRARIES})

if( TARGET examples )
    add_dependencies( examples sudoku )
endif()

#add_subdirectory(check)
#enable_testing()
