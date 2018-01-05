cmake_minimum_required(VERSION 3.0)

project(Rinpacking)

find_package(SCIP REQUIRED)
include_directories(${SCIP_INCLUDE_DIRS})

add_executable(rinpacking
   src/cmain.c
   src/reader_rpa.c
   src/probdata_rpa.c
   src/pattern.c)

target_link_libraries(rinpacking ${SCIP_LIBRARIES})
if( TARGET examples )
    add_dependencies( examples rinpacking )
endif()

#
# add check subdirectory for tests
#
#add_subdirectory(check)

enable_testing()
