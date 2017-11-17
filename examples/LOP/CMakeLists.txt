cmake_minimum_required(VERSION 3.0)

project(LOP)
include(CTest)
find_package(SCIP REQUIRED)
include_directories(${SCIP_INCLUDE_DIRS})

add_executable(lop
   src/cmain.c
   src/cons_lop.c
   src/reader_lop.c)

target_link_libraries(lop ${SCIP_LIBRARIES})

if( TARGET examples )
    add_dependencies( examples lop )
endif()

add_subdirectory(check)
enable_testing()
