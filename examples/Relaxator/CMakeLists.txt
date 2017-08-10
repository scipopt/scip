project(Relaxator)

find_package(SCIP REQUIRED)
include_directories(${SCIP_INCLUDE_DIRS})

add_executable(relaxator
   src/cmain.c
   src/relax_lp.c
   src/relax_nlp.c)

target_link_libraries(relaxator ${SCIP_LIBRARIES})
add_dependencies(examples relaxator)

add_subdirectory(check)
enable_testing()
