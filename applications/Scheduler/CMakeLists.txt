cmake_minimum_required(VERSION 3.0)

project(Scheduler)

find_package(SCIP REQUIRED)

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(scheduler
    src/main.cpp
    src/cons_optcumulative.c
    src/heur_listscheduling.c
    src/heur_optcumulative.c
    src/reader_cmin.c
    src/reader_rcp.c
    src/reader_sch.c
    src/reader_sm.c
)

target_link_libraries(scheduler ${SCIP_LIBRARIES})

if( TARGET applications )
    add_dependencies( applications scheduler )
endif()

add_subdirectory(check)

enable_testing()
