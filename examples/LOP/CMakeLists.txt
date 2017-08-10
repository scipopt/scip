project(LOP)
include(CTest)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src)

add_executable(lop
   src/cmain.c
   src/cons_lop.c
   src/reader_lop.c)

target_link_libraries(lop libscip)
add_dependencies(examples lop)

add_subdirectory(check)
enable_testing()
