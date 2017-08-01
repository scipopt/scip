project(CallableLibrary)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src)

add_executable(circle
   src/circle.c)

add_executable(gastrans
   src/gastrans.c)

add_executable(string
   src/string.c)

target_link_libraries(circle libscip)
target_link_libraries(gastrans libscip m)
target_link_libraries(string libscip)

add_dependencies(examples circle gastrans string)
