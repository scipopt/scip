cmake_minimum_required(VERSION 3.0)

project(MinIISC)

find_package(SCIP REQUIRED)

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(miniisc
    src/benders.c
    src/miniisc.c
    src/readargs.c
)

target_link_libraries(miniisc ${SCIP_LIBRARIES} m)

if( TARGET applications )
    add_dependencies( applications miniisc )
endif()

add_subdirectory(check)

enable_testing()
