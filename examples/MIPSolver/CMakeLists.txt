project(MIPSolver)

find_package(SCIP REQUIRED)
include_directories(${SCIP_INCLUDE_DIRS})

add_executable(scipmip
   src/cppmain.cpp)

target_link_libraries(scipmip ${SCIP_LIBRARIES})
add_dependencies(examples scipmip)
