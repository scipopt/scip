cmake_minimum_required(VERSION 3.0)

project(CAP)

find_package(SCIP REQUIRED)
include_directories(${SCIP_INCLUDE_DIRS})

add_executable(cap
   src/cmain.c
   src/probdata_cap.c
   src/reader_cap.c )

target_link_libraries(cap ${SCIP_LIBRARIES})
target_link_libraries(cap m)
if( TARGET examples )
    add_dependencies( examples cap )
endif()

#
# add check subdirectory for tests
#
add_subdirectory(check)

enable_testing()
