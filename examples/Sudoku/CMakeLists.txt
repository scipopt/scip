cmake_minimum_required(VERSION 3.3)

project(Sudoku)

find_package(SCIP REQUIRED)
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
