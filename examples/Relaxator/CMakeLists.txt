cmake_minimum_required(VERSION 3.0)

project(Relaxator)

find_package(SCIP REQUIRED)
include_directories(${SCIP_INCLUDE_DIRS})

add_executable(relaxator
   src/cmain.c
   src/relax_lp.c
   src/relax_nlp.c)

target_link_libraries(relaxator ${SCIP_LIBRARIES})

if( TARGET examples )
    add_dependencies( examples relaxator )
endif()

add_subdirectory(check)
enable_testing()
