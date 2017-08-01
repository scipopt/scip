project(GMI)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src)

add_executable(gmi
   src/cmain.c
   src/sepa_gmi.c)

target_link_libraries(gmi libscip m)
add_dependencies(examples gmi)
