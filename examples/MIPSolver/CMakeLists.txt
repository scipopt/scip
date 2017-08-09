project(MIPSolver)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src)

add_executable(scipmip
   src/cppmain.cpp)

target_link_libraries(scipmip libscip m)
add_dependencies(examples scipmip)
