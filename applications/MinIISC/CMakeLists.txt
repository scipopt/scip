project(MinIISC)

find_package(SCIP REQUIRED)

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(miniisc 
    src/benders.c
    src/miniisc.c
    src/readargs.c
)

target_link_libraries(miniisc ${SCIP_LIBRARIES} m)
                                                 
add_dependencies(applications miniisc)

add_subdirectory(check)

enable_testing()
