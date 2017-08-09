project(Relaxator)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src)

add_executable(relaxator
   src/cmain.c
   src/relax_lp.c
   src/relax_nlp.c)

target_link_libraries(relaxator libscip)
add_dependencies(examples relaxator)

add_subdirectory(check)
enable_testing()
