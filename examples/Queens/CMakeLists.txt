project(Queens)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src)

add_executable(queens
   src/queens.cpp
   src/queens_main.cpp)

target_link_libraries(queens libscip)
add_dependencies(examples queens)
