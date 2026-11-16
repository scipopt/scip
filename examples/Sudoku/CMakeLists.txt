cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
   cmake_policy(SET CMP0069 NEW)
endif()

project(Sudoku)

# Find SCIP unless we already have it because we build it.
if(NOT TARGET SCIP::SCIP)
  find_package(SCIP
     REQUIRED
     CONFIG
     HINTS "${SCIP_DIR}"
     )
endif()

add_executable(sudoku
   src/sudoku_main.cpp
   src/sudoku_utils.h)
target_link_libraries(sudoku SCIP::SCIP)

if(SCIP_HAS_LTO)
   set_target_properties(sudoku PROPERTIES
      INTERPROCEDURAL_OPTIMIZATION 1)
endif()

if(TARGET examples)
    add_dependencies(examples sudoku)
endif()

