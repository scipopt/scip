project(VRP)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src)

add_executable(vrp
   src/main_vrp.cpp
   src/pricer_vrp.cpp)

target_link_libraries(vrp libscip)
add_dependencies(examples vrp)


add_subdirectory(check)
enable_testing()
