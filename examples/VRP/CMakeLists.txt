project(VRP)

find_package(SCIP REQUIRED)
include_directories(${SCIP_INCLUDE_DIRS})

add_executable(vrp
   src/main_vrp.cpp
   src/pricer_vrp.cpp)

target_link_libraries(vrp ${SCIP_LIBRARIES})
add_dependencies(examples vrp)


add_subdirectory(check)
enable_testing()
